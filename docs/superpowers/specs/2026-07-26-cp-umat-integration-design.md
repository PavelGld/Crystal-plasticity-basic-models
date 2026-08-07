# Дизайн: CP-UMAT – Abaqus Standard Integration

Дата: 2026-08-07
Ветка: `implicit-integration` (поверх неявной схемы)
Статус: согласовано, ожидает реализации

## 1. Контекст

Неявная схема интегрирования (Newton–Raphson) с консистентным тангенсом реализована в `implicit-integration`. Следующий шаг — встраивание в Abaqus/Standard как UMAT. Первый этап: **одно зерно на точку интегрирования** (RVE-per-IP — следующим этапом).

Референс: работающий EMSC UMAT в `UMAT/` (нейросетевая модель) — шаблон для сборки MSVC → `.obj`, extern "C" входа, STATEV/PROPS/DDSDDE.

## 2. Решения по скоупу

- **Одно зерно на IP** (v1). RVE-per-IP — следующим этапом (RVE.hpp готовится header-only в cp_model/).
- **Header-only рефакторинг** — модель переносится в `cp_model/`, `.h`+`.cpp` сливаются в `.hpp`.
- **Интерфейсный слой** — конвертация Abaqus↔Voigt только в `cp_umat/`, модель не знает про Abaqus.
- **NLGEOM=NO** (малые деформации). Большие — следующим этапом (нужен DFGRD + DROT).
- **FCC {111}⟨110⟩** 12 систем зашиты в код (не в PROPS). General — следующим этапом.
- **Нет файлового I/O** (параметры только из PROPS), нет cout/cin.
- **Сборка** — чистый MSVC (как EMSC-референс), `/EHsc /O2`, один `.cpp` → один `.obj`.

## 3. Структура проекта после рефакторинга

```
EigenProject/
├── cp_model/                         ← библиотека модели (header-only)
│   ├── tensor_utils.hpp              (переименован из tensor_utils.h)
│   ├── MaterialParams.hpp            (слияние MaterialParams.h + MaterialParams.cpp)
│   ├── Grain.hpp                     (слияние Grain.h + Grain.cpp)
│   └── RVE.hpp                       (слияние RVE.h + RVE.cpp)
│
├── cp_umat/                          ← UMAT-специфичный код
│   ├── build_cp_umat.bat            (MSVC сборка → cp_umat.obj)
│   ├── cp_umat.cpp                  (extern "C" UMAT + UMAT_)
│   ├── cp_umat_interface.hpp        (класс CPUMAT)
│   └── cp_umat_types.hpp            (Abaqus↔Voigt матрицы 6×6)
│
├── cp_umat_validation/               ← .inp для Abaqus
│   ├── tension.inp                   (одноосное растяжение)
│   ├── compression.inp               (одноосное сжатие)
│   ├── shear.inp                     (простой сдвиг)
│   └── cyclic.inp                    (циклическое нагружение)
│
├── tests.cpp                         (обновлённые #include на cp_model/)
├── EigenProject.cpp                  (обновлённые #include)
└── удалены из корня:
    tensor_utils.h, MaterialParams.h/.cpp, Grain.h/.cpp, RVE.h/.cpp
```

Сборка тестов: `cl ... tests.cpp` (один файл, модель header-only).

## 4. Формат CPUMAT (cp_umat_interface.hpp)

### Матрицы перехода (cp_umat_types.hpp)

Abaqus NTENS=6: [σ₁₁,σ₂₂,σ₃₃, σ₁₂,σ₁₃,σ₂₃] (тензорные, без ×2).
Мой Voigt: [σ₀₀,σ₁₁,σ₂₂, 2σ₂₃,2σ₁₃,2σ₁₂].

```
A_to_V: Abaqus Δε → Voigt Δε (6×6, пермутация + ×2 на сдвигах)
V_to_A: Voigt Δσ → Abaqus Δσ (6×6, пермутация + /2 на сдвигах)
DDSDDE_abaqus = V_to_A · T · A_to_V   (T = тангенс в Фойгте, КСК)
```

### PROPS (14-19 чисел)

| PROPS[i] | Значение |
|----------|----------|
| 0 | NPROPS (проверка) |
| 1-3 | c11, c12, c44 (Па) |
| 4 | hard_mode (0,1,2) |
| 5 | rot_mode (0,1,2) |
| 6 | integ_mode (должен быть 1) |
| 7 | tc0 (Па) |
| 8 | strain_rate_infl (γ̇₀) |
| 9 | n_power |
| 10-13 | qw, qx, qy, qz (начальный кватернион) |
| +5 если hard_mode=2 | A, q, h0, ts, alfa |
| +2 если hard_mode=1 | a, b |
| +0 если hard_mode=0 | — |

### STATEV (70 чисел)

| STATEV[i] | Содержимое |
|-----------|-----------|
| 0 | флаг инициализации (123456.0 = инициализировано) |
| 1-24 | tc[0..23] |
| 25-48 | shifts[0..23] |
| 49-52 | orient_q (w,x,y,z) |
| 53-58 | stress_v (Voigt, кристальный фрейм) |
| 59-67 | F (3×3, row-major) |
| 68 | acc_def |
| 69 | sim_time |

### Алгоритм CPUMAT::step

1. Преобразовать DSTRAN (Abaqus) → D_v (Voigt) через A_to_V
2. Построить L_3×3 = voigtToTensor(D_v / dt) (чистая деформация, без спина)
3. Вызвать Grain::Simulate_Step(dt, L)
4. Прочитать σ_lab = Grain::get_stress_lab()
5. Преобразовать σ_lab → Abaqus STRESS через V_to_A
6. Построить DDSDDE = V_to_A · T_crystal · A_to_V (через get_tangent_crystal_into)
7. При несходимости Ньютона → PNEWDT = 0.5
8. Сохранить состояние в STATEV

## 5. План задач

### Task A: Header-only рефакторинг (cp_model/)
- Переименовать/слить файлы
- Обновить includes в tests.cpp и EigenProject.cpp
- Проверить: сборка тестов, 16/16 зелёные

### Task B: Abaqus-интерфейс (cp_umat/)
- cp_umat_types.hpp (матрицы)
- cp_umat_interface.hpp (CPUMAT)
- cp_umat.cpp (extern "C" UMAT)
- build_cp_umat.bat
- Проверить: компиляция cp_umat.obj

### Task C: Валидация и документация
- 4 .inp файла в cp_umat_validation/
- README.md в cp_umat/
- Инструкция по запуску в Abaqus

### Гейты приёмки
- [ ] tests.exe — 16/16 после рефакторинга
- [ ] build_cp_umat.bat компилирует cp_umat.obj
- [ ] tension.inp запускается в Abaqus без ошибок
