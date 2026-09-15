# CP-UMAT Review Request: Multi-Element Convergence Problem

Дата: 2026-08-12
Автор: [project maintainer]
Статус: просьба о ревью

## 1. Проблема (одним абзацем)

UMAT кристаллической пластичности (неявный Newton-Raphson по напряжениям, 6×6) **сходится на одноэлементных тестах и на многоэлементной сетке с одинаковой ориентацией всех зёрен**, но **падает (TOO MANY ATTEMPTS) на многоэлементной сетке, где у каждого элемента своя случайная начальная ориентация** (материал-per-element, каждый элемент — отдельный `*MATERIAL` с одинаковыми 19 PROPS, отличающимися только кватернионом PROPS[10-13]).

Abaqus режет шаг от 0.005 до 3e-4 и всё равно не сходится за 5 попыток. Ошибка: `***ERROR: TOO MANY ATTEMPTS MADE FOR THIS INCREMENT` в глобальном Ньютоне (equilibrium iterations).

## 2. Что работает (контрольные точки)

| Тест | Ориентации | Результат |
|------|-----------|-----------|
| Single element, `tension.inp` | identity (1,0,0,0) | ✅ S33=49.78 МПа, 200 inc, 0 cutbacks |
| Single element, `tension_random.inp` | один случайный кватернион | ✅ S33=25.50 МПа (анизотропия видна), но Abaqus режет шаг |
| 8 элементов, `tension_8el.inp` | все identity | ✅ S33=49.78 МПа равномерно |
| 8 элементов, `tension_8el_random.inp` (material-per-element) | все identity | ✅ S33=49.78 МПа равномерно |
| 8 элементов, `tension_8el_random.inp` (material-per-element) | случайные | ❌ TOO MANY ATTEMPTS на инкременте 1 |

Наблюдение: один элемент с НЕidentity кватернионом тоже требует резания шага (200 инкрементов не хватает до ε=0.1, останавливается на 0.75). Это намекает, что проблема — в тангенсе (DDSDDE) для повёрнутой ориентации, а не в архитектуре.

## 3. Архитектура UMAT

Один `Grain*` (singleton `static CPUMAT`), переиспользуется между всеми точками интегрирования. Состояние каждой точки (stress, tc, shifts, ориентация) сериализуется через STATEV (70 значений):

```
STATEV[0-5]    = stress_v (Voigt, кристальный фрейм)
STATEV[6-9]    = кватернион ориентации (w,x,y,z)
STATEV[10-33]  = tc[24] — критические напряжения
STATEV[34-57]  = shifts[24] — накопленные сдвиги
STATEV[67]     = sim_time
STATEV[68]     = acc_def
STATEV[69]     = guard 123456.0 (инициализировано?)
```

На каждом вызове `umat_()`:
1. Если guard не выставлен (первый заход в точку): читаем кватернион из PROPS[10-13], `set_orientation`, `reset_state`, сохраняем в STATEV.
2. Иначе: `load_state(STATEV)` → восстанавливаем состояние.
3. `step()`: `DSTRAN → A_to_V → voigtToTensor → L = D/dt → Grain::Simulate_Step` → `get_stress_lab` → `STRESS`, `get_tangent_lab_into` → `DDSDDE`.
4. `save_state(STATEV)`.

Grain вычисляет консистентный тангенс `J⁻¹(P − c·gᵀ)` на сошедшемся состоянии (crystal frame), затем поворачивается в lab frame через `T_voigt` (`get_tangent_lab_into`).

## 4. Подозрение №1: кэш `T_voigt` при смене ориентации

В `Grain`:

```cpp
mutable Matrix6d T_voigt_buf;
mutable bool T_voigt_valid = false;
void computeVoigtRotationInPlace(Matrix6d& T) const;  // T от R_buf

void Grain::get_elasticity_lab_into(Matrix6d& result) const {
    if (!T_voigt_valid) {
        computeVoigtRotationInPlace(T_voigt_buf);
        T_voigt_valid = true;
    }
    Buf_P.noalias() = T_voigt_buf * params->P_v;
    result.noalias() = Buf_P * T_voigt_buf.transpose();
}
```

`T_voigt_valid` инвалидируется в `set_orientation` и `updateOrientation`. Но между точками интеграции объект `Grain` переиспользуется. Если для точки A (identity) `T_voigt` вычислен и закэширован, а затем точка B (случайная ориентация) вызывает `load_state` → `set_orientation(q_B)` → должно инвалидировать кэш. **Проверить: не забыта ли инвалидация где-то в пути, который пропускает `set_orientation` (например, `reset_state`).**

`reset_state()` сейчас НЕ трогает `orient_q` и `T_voigt_valid`. Flow первого захода: `set_orientation(q)` (инвалидирует) → `reset_state()` (не трогает). Это корректно. Но на ПОВТОРНОМ заходе в новую точку (guard уже выставлен, но это другая IP) — `load_state` вызывает `set_orientation` из STATEV. Корректно.

**НО:** возможно проблема в том, что на ПЕРВОМ инкременте Abaqus вызывает UMAT для всех 8 IP с разными PROPS (разные кватернионы), а `static bool ready` уже выставлен после первой точки. `init_params` вызывается один раз с PROPS первой точки — это ок (материал одинаков). Ориентация читается из PROPS[10-13] КАЖДОЙ точки на её первом заходе. Ок.

## 5. Подозрение №2: упругий тангенс на первом шаге не согласован

На первом инкременте, ДО первого `Simulate_Step`, `tangent_valid == false`. `get_tangent_lab_into` под `ABAQUS_UMAT` возвращает **упругую** матрицу (`get_elasticity_lab_into`) вместо тангенса. Для identity-ориентации упругая матрица в lab frame = cubic stiffness. Для случайной — повёрнутая cubic stiffness.

Возможная проблема: Abaqus на первой итерации использует DDSDDE для сборки глобальной матрицы жёсткости и решает систему. Если у некоторых элементов упругая матрица сильно отличается (жёсткая/мягкая), глобальная система плохо обусловлена → первый шаг плохой → итерации расходятся.

Но identity-8el работает, значит сборка упругой матрицы не является катастрофой сама по себе. Разница — разброс жёсткостей между элементами.

## 6. Подозрение №3: `A_to_V` ×2 на сдвигах + `P_v(3,3) = c44`

В `cp_umat_types.hpp`:

```cpp
A_to_V:  Voigt[3]=2ε23 ← Abaqus[5]=ε23 (×2)   // для DSTRAN
```

Коллега ранее указал: `DSTRAN` в Abaqus для сдвигов — УЖЕ инженерная деформация (γ=2ε), та же, что наш Voigt. Значит `A_to_V` должна быть чистой перестановкой БЕЗ ×2. Сейчас двойка компенсируется тем, что `P_v(3,3)=c44` вместо `2*c44` — но только в упругом режиме. Для повёрнутых ориентаций и пластики двойки НЕ сокращаются в `computeDeviatoricDu()` (du завышен для сдвиговых компонент → rate_factor завышен → неправильная пластическая кинетика).

Это влияет на КОРРЕКТНОСТЬ отклика, но вряд ли на СХОДИМОСТЬ. Отмечено для исправления, но не является первопричиной краха.

## 7. Что уже пробовали (не помогло)

1. **`UNSYMM` флаг** в `*USER MATERIAL` — не помог (глобальный Ньютон всё равно не сходится).
2. **Уменьшение шага** (0.005 → 0.0003) — не помогло, Abaqus режет шаг и всё равно 5 попыток не хватает.
3. **Material-per-element vs единый материал с таблицей кватернионов в PROPS** — первый работает (парсинг), второй падает на парсере из-за ограничения "8 значений на строку" (решаемо, но не суть).

## 8. Просьба к ревьюеру

Посмотреть на:

1. **`cp_umat/cp_umat.cpp`** — логика `umat_()` (инициализация, load/save, step). Особенно: не теряется ли состояние между точками интеграции, корректна ли изоляция при переиспользуемом `Grain*`.
2. **`cp_model/Grain.hpp`** — `solve_implicit`, `assemble_J`, `evaluate_Dp`, `get_tangent_lab_into`, `get_elasticity_lab_into`, `reset_state`, кэш `T_voigt`. Ищем, почему для повёрнутой ориентации тангенс/стресс рассинхронизируются с глобальным Ньютоном.
3. **`cp_umat/cp_umat_types.hpp`** — `A_to_V`, `V_to_A` (п.6).
4. **`cp_umat_validation/gen_per_material.py`** — генератор `.inp` (корректность формата).

Вопросы:
- Есть ли ошибка в консистентном тангенсе для повёрнутых ориентаций (что-то в повороте `T_voigt` при неidentity-кватернионе)?
- Не является ли проблема в том, что `get_stress_lab()` и `get_tangent_lab_into()` используют РАЗНЫЕ состояния `R_buf`/`T_voigt` после шага (гонка между обновлением ориентации и кэшем)?
- Правильно ли вообще формулировать DDSDDE для вязкопластической модели с явной зависимостью γ̇ от D_eq (member `rate_factor_impl`)?

## 9. Файлы для ревью

Основные:
- `cp_umat/cp_umat.cpp`
- `cp_umat/cp_umat_interface.hpp`
- `cp_umat/cp_umat_types.hpp`
- `cp_model/Grain.hpp`
- `cp_model/MaterialParams.hpp`
- `cp_model/tensor_utils.hpp`

Валидация:
- `cp_umat_validation/tension.inp` (рабочий, identity)
- `cp_umat_validation/tension_random.inp` (рабочий, 1 элемент, случайный)
- `cp_umat_validation/tension_8el.inp` (рабочий, 8 элементов, identity)
- `cp_umat_validation/tension_8el_random.inp` (НЕ работает — генерится `gen_per_material.py`)
- `cp_umat_validation/gen_per_material.py` (генератор)

Документация:
- `cp_umat/README.md` (архитектура, PROPS/STATEV)
- `cp_umat/DEBUG_CONTEXT.md` (история отладки)

Зависимость: Eigen 3.4 (внешний, `D:\eigen-3.4.0`), MSVC 2022, Abaqus 2022.