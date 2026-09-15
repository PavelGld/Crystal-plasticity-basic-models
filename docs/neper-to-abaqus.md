# Интеграция сеток Neper в Abaqus (CPFEM)

Пошаговая инструкция: как из структуры Neper получить готовый Abaqus-дек с
кристаллической пластичностью и запустить расчёт.

---

## 0. TL;DR

```bash
# 1. собрать UMAT (один раз)
cd cp_umat && build_cp_umat.bat && cd ..

# 2. сгенерировать дек из сетки Neper
python scripts/neper2abaqus.py meshes/neper_200gr_11135el out.inp \
       1.0 0.0001 0.001 mesh_c free

# 3. запустить
abaqus job=out user=cp_umat/cp_umat.obj cpus=4 mp_mode=threads \
       standard_parallel=all int
```

---

## 1. Что выгружает Neper

Команда Neper с `-M` создаёт набор файлов (пример — `meshes/neper_200gr_11135el/`):

| Файл | Содержимое | Нужен? |
|---|---|---|
| **`mesh.inp`** | Сетка в формате Abaqus: узлы, элементы, наборы `poly1…polyN` (по одному на зерно) | **ДА — основной** |
| **`mesh.ori`** | Ориентации зёрен: углы Бунге (секция `$ElsetOrientations`) | **ДА** |
| `mesh.msh` | Сетка в формате Gmsh | ❌ не использовать (см. §2) |
| `mesh.stelset` | По зёрнам: `id, число_элементов, объём, φ1, Φ, φ2` | для сверки |
| `mesh.stelt` | По элементам: `id_элемента, id_зерна` | для сверки |

> ⚠️ **Главная ловушка:** порядок узлов у тетраэдров `C3D10` в `mesh.msh` и
> `mesh.inp` **разный** (последние два узла n9/n10 переставлены). Если собрать
> дек из `.msh`, все элементы получат отрицательный объём и Abaqus упадёт на
> этапе Flattener с `Assertion Failure: i >= 0 && i < nrows`.
> **Всёгда берите `mesh.inp`.**

Проверить, что сетка согласована (номер зерна в `stelset` совпадает с углами):

```bash
# число элементов должно дать общее количество
awk '{s+=$2} END {print s}' meshes/neper_200gr_11135el/mesh_c.stelset
```

---

## 2. Генератор дека: `scripts/neper2abaqus.py`

### Синтаксис

```bash
python scripts/neper2abaqus.py <mesh_dir> <out.inp> \
       [strain] [dt_start] [dt_max] [mesh_basename] [bc_mode]
```

| Аргумент | По умолчанию | Смысл |
|---|---|---|
| `mesh_dir` | — | папка с `mesh.inp` + `mesh.ori` |
| `out.inp` | — | выходной дек |
| `strain` | `0.3` | конечная инженерная деформация (1.0 = 100%) |
| `dt_start` | `0.001` | стартовый шаг деформации за инкремент |
| `dt_max` | `10·dt_start` | максимальный шаг (Abaqus растёт сам) |
| `mesh_basename` | `mesh` | префикс файлов (`mesh` → `mesh.inp`; `mesh_c` → `mesh_c.inp`) |
| `bc_mode` | `sym` | граничные условия: `sym` / `free` / `clamp` |

### Что делает генератор

1. Читает `mesh.inp` → узлы, элементы, наборы зёрен `polyN`.
2. Читает `mesh.ori` → углы Бунге `(φ1, Φ, φ2)` → кватернион.
3. **Проверяет**: все ссылки на узлы/элементы существуют, знак объёма каждого
   тетраэдра. Печатает:
   ```
   nodes=2213 elements=11135 poly-elsets=200 orientations=200
   missing node refs=0, missing elem refs=0
   negative/zero-volume corner tets: 0 / 11135
   ```
   Если `negative/zero-volume` > 0 — сетка или её источник неверны (§1).
4. Пишет **плоский** дек (без `*Part`/`*Instance`):
   - `*NODE`, `*ELEMENT, TYPE=C3D4|C3D8|C3D10` (тип определяется по числу узлов)
   - на каждое зерно: `*ELSET polyN` + `*SOLID SECTION` + `*MATERIAL`
   - `*NSET` для граней
   - `*STEP` с ГУ

### Формат материала (19 констант)

```
*MATERIAL, NAME=GRAIN_1
*USER MATERIAL, CONSTANTS=19, UNSYMM
0, 106750, 60410, 28340, 2, 2, 1, 6
0.588235, 50, <qw>, <qx>, <qy>, <qz>, 1, 2
115, 34, 2.25
*DEPVAR
70
```

| Строка | PROPS | Смысл |
|---|---|---|
| 1 | `[0..7]` | filler, c11, c12, c44, hard_mode(=2), rot_mode(=2), integ_mode(=1), tc0 |
| 2 | `[8..15]` | strain_rate_infl, **n_power(=50)**, **кватернион [10..13]**, A(=1), q(=2) |
| 3 | `[16..18]` | **h0(=115), ts(=34), alfa(=2.25)** — Kocks-Mecking |

Кватернион — поворот КСК→ЛСК, из углов Бунге. Хардкод параметров — в
`neper2abaqus.py` (функция `main`, блок генерации `*MATERIAL`).

---

## 3. Граничные условия (`bc_mode`)

| Режим | Схема | Когда использовать |
|---|---|---|
| **`sym`** | 3 грани на плоскостях симметрии (1/8 модели), 4-я растягивается | быстрый расчёт, если нужна только текстура |
| **`free`** | 1 грань `ux=0`, 2 узла убирают твёрдое тело, 4 боковые грани **свободны** | честное одноосное растяжение RVE |
| **`clamp`** | 1 грань закреплена **полностью** (ux=uy=uz=0) | имитация жёсткого захвата; даёт неоднородное поле |

**`sym`** (1/8 модели):
```
LEFT   (x=0)  ux = 0
BOTTOM (y=0)  uy = 0
BACK   (z=0)  uz = 0
RIGHT  (x=L)  ux = strain·L
```

**`free`** (честное растяжение):
```
LEFT   (x=0)  ux = 0
узел-угол     uy = uz = 0     ← убираем y,z-трансляцию
узел y=max    uz = 0          ← убираем поворот вокруг x
RIGHT  (x=L)  ux = strain·L
(грани y=0, y=L, z=0, z=L — свободны)
```

**`clamp`** (жёсткий захват):
```
LEFT   (x=0)  ux = uy = uz = 0
RIGHT  (x=L)  ux = strain·L
```

Что выбрать: если нужна **эволюция текстуры** — берите `free` (или `sym` для
скорости). Если хотите увидеть **неоднородность** (концентрацию у зажима,
градиент поперечной деформации) — `clamp`. Подробнее — в
[`examples/04_neper_200grains_free/README.md`](../examples/04_neper_200grains_free/README.md).

---

## 4. Настройки шага (важно)

Генератор пишет один шаг на всю деформацию:

```
*STEP, NLGEOM=YES, INC=200000, UNSYMM=YES
*STATIC
<dt_start>, <strain>, 1e-9, <dt_max>
*CONTROLS, PARAMETERS=TIME INCREMENTATION
8, 10, 9, 16, 10, 4, 12, 5
*CONTROLS, PARAMETERS=LINE SEARCH
5, 1.0, 0.0001, 0.25, 0.10
```

Три вещи, без которых расчёт либо падает, либо встаёт:

1. **`UNSYMM=YES`** — консистентный тангенс несимметричен.
2. **`*CONTROLS, PARAMETERS=LINE SEARCH`** — в Abaqus line search **выключен**
   по умолчанию. Без него первый инкремент не сходится.
3. **Один `*STEP` на всю деформацию.** Каждый новый `*STEP` Abaqus начинает без
   экстраполяции и с пробным вызовом UMAT при `DSTRAN=0` (упругий тангенс) —
   на пластически нагруженной модели это срезает `dt` и всё замедляет.

Плюс: UMAT **не должен** безусловно ставить `PNEWDT=1.0` — Abaqus берёт минимум
PNEWDT по всем вызовам как множитель нового `dt`, так что «1.0» запрещает рост
шага (см. `cp_umat/cp_umat.cpp`, функция `CPUMAT::step`).

---

## 5. Сборка и запуск

### UMAT (один раз)
```bash
cd cp_umat
build_cp_umat.bat          # -> cp_umat/cp_umat.obj
```
Требуется MSVC 2022 и Eigen (путь `D:\eigen-3.4.0` в `build_cp_umat.bat`).

### Запуск
```bash
abaqus job=<jobname> user=cp_umat/cp_umat.obj cpus=4 mp_mode=threads \
       standard_parallel=all int
```
> Для линковки Abaqus нужен Intel oneAPI в окружении (`setvars.bat`).
> Готовый `run.bat` в каждом примере это делает сам.

### Запуск через `run.bat` примера
```bash
cd examples/03_neper_200grains_symmetry
run.bat
```

---

## 6. Извлечение результатов

Все скрипты — в `scripts/`. Запуск через `abaqus python` (кроме графики).

| Задача | Команда |
|---|---|
| Напряжения S33 по элементам | `abaqus python scripts/extract_s33.py <job> s33.txt` |
| Ориентации решётки (текстура) | `abaqus python scripts/extract_quats.py <job> quats.txt` |
| Перемещения (для анализа ГУ) | `abaqus python scripts/extract_u.py <job> u.txt` |
| Контроль объёма (непожимаемость) | `abaqus python scripts/check_volume.py <job>` |
| Полюсные фигуры | `python scripts/polefigures.py quats.txt PF.png "title" 0` |

Последний аргумент `polefigures.py` — ось проекции: `0`=X (растяжение),
`1`=Y, `2`=Z.

Пример полного цикла для примера 03:
```bash
cd examples/03_neper_200grains_symmetry
run.bat
abaqus python ../../scripts/extract_quats.py model quats.txt
python ../../scripts/polefigures.py quats.txt PF.png "200 grains, 100%" 0
```

### Что смотреть в ODB

| Поле | Смысл |
|---|---|
| `S` | напряжения (6 компонент) |
| `LE` | логарифмическая деформация (`tr(LE) = ln J`) |
| `SDV7…10` | кватернион ориентации решётки (текстура) |
| `SDV11…34` | критическое напряжение `tc` по системам |
| `SDV35…58` | накопленные сдвиги |
| `U` | перемещения |

> **Внимание:** Abaqus может писать SDV как **отдельные поля** `SDV1…SDV70`,
> а не одним `SDV`. Скрипты это учитывают.

Интерактивно: `abaqus cae database=<job>.odb`.

---

## 7. Производительность: как выбрать сетку

Стоимость ∝ числу элементов × числу точек интегрирования (`C3D4` — 1 IP,
`C3D8` — 8, `C3D10` — 4) и сильно зависит от размера задачи для решателя.

Ориентир (200 зёрен, n_power=50, `NLGEOM=YES`, cpus=4, до 100% деформации):

| Сетка | Элементов | Элементов/зерно | Время до 100% |
|---|---|---|---|
| C3D10 мелкая | 46 499 | 232 | ~50 ч |
| C3D4 грубая | 11 135 | 56 | **~8 мин** |

Для **эволюции текстуры** плотная сетка избыточна — важна статистика по
зёрнам, а не внутризёренные градиенты. ~50 элементов/зерно достаточно.
Сгенерировать грубее в Neper: уменьшить `-rcl` (например `-rcl 0.7`).

---

## 8. Чек-лист

- [ ] UMAT собран (`cp_umat/cp_umat.obj` существует)
- [ ] Есть `mesh.inp` и `mesh.ori` из Neper (не `.msh`!)
- [ ] `neper2abaqus.py` напечатал `negative/zero-volume corner tets: 0 / N`
- [ ] В деке: `UNSYMM=YES`, `*CONTROLS ... LINE SEARCH`, **один** `*STEP`
- [ ] UMAT не пиннит `PNEWDT=1.0`
- [ ] Выбран подходящий `bc_mode`
- [ ] Abaqus запущен с Intel oneAPI в окружении
- [ ] `.sta` заканчивается на `THE ANALYSIS HAS COMPLETED SUCCESSFULLY`
- [ ] Текстура извлечена и построены полюсные фигуры

---

## 9. Типичные ошибки

| Симптом | Причина | Лечение |
|---|---|---|
| `Assertion Failure: i >= 0 && i < nrows`, Flattener | дек собран из `mesh.msh` | использовать `mesh.inp` |
| `The volume of N elements is zero, small, or negative` | неверный порядок узлов | использовать `mesh.inp` |
| `8 NEGATIVE EIGENVALUES`, невязка не падает | нет line search | добавить `*CONTROLS, PARAMETERS=LINE SEARCH` |
| `dt` навсегда равен стартовому | UMAT ставит `PNEWDT=1.0` | убрать (см. §4) |
| Первый инкремент нового шага падает | Abaqus стартует шаг без экстраполяции | **один** `*STEP` на всю деформацию |
| `USER SUBROUTINE ... STATEV DIMENSIONED TO ZERO` | нет `*DEPVAR` | добавить `*DEPVAR` / `NSTATV=70` |
| Дек не линкуется | нет Intel oneAPI | `setvars.bat` перед `abaqus` |
