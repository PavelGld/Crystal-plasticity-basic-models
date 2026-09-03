# Crystal Plasticity Grain Model

C++ реализация модели кристаллической пластичности для моделирования деформации и эволюции текстуры поликристаллических материалов.

## Возможности

- Упруговязкопластическая модель для монокристалла (закон Хатчинсона)
- Поддержка произвольного числа систем скольжения
- Три режима упрочнения:
  - **Без упрочнения** (hard_mode=0)
  - **Изотропное** (hard_mode=1): `tc = tc0 + a · (Σγ)^b`
  - **Kocks-Mecking** (hard_mode=2): матрица взаимодействия H
- Два режима ротации кристаллической решётки:
  - **Материальный поворот** (Mandel spin)
  - **Стеснённый поворот Тейлора** (Taylor constrained rotation)
- Представление ориентации через кватернионы
- **RVE** (Representative Volume Element) — осреднение по множеству зёрен
- **Жёсткое и мягкое нагружение** (заданный L или заданный dS/dt)
- **Параллелизация** через OpenMP
- **Неявная схема интегрирования** (Newton-Raphson) с консистентным тангенсом — для больших шагов и встраивания в Abaqus/Standard (UMAT)
- **Multi-element UMAT** — thread-safe реализация с `__declspec(thread)`, валидирована на 8 и 343 элементах
- **Параллельное выполнение в Abaqus** — ускорение ~3.7x на 4 ядрах (343 элемента: 48 мин vs ~3 часа)
- Валидация текстур через MTEX pole figures

## Валидация текстур

Модель прошла валидацию для FCC материала (Cu, Al) — 12 систем скольжения {111}<110>. Полюсные фигуры на 100% деформации совпали с литературными данными:

### Одноосное сжатие (100%)
Текстура `<110>||ND` — классическая для FCC compression:

![Compression {111}](rotation_validation_test_results/PF_uniaxial_compression_100pct.png)

### Одноосное растяжение (100%)
Текстура `<111>||ND` — классическая для FCC tension:

![Tension {111}](rotation_validation_test_results/PF_uniaxial_tension_100pct.png)

### Плоское растяжение (100%)
Brass-type текстура `{110}<112>`:

![Plane Strain {111}](rotation_validation_test_results/PF_plane_strain_100pct.png)

### Простой сдвиг (100%)

![Simple Shear {111}](rotation_validation_test_results/PF_simple_shear_100pct.png)

## Архитектура

```
┌─────────────────────────────────────────────────────────────┐
│                        MaterialParams                        │
│  (immutable, shared_ptr между зёрнами одного материала)      │
├─────────────────────────────────────────────────────────────┤
│  M[i], Wpm[i], Ms_v[i]  — системы скольжения                 │
│  P_v                    — матрица упругости (6×6)           │
│  tc0, strain_rate_infl  — параметры пластичности             │
│  hard_mode, H, A, q...  — параметры упрочнения               │
│  iso_a, iso_b           — изотропное упрочнение (hard_mode=1)│
│  rot_mode               — режим ротации                     │
└─────────────────────────────────────────────────────────────┘
                              │
                              │ shared_ptr
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                         Grain                                │
├─────────────────────────────────────────────────────────────┤
│  orient_q    — кватернион ориентации (КСК → ЛСК)             │
│  Stress_v    — напряжения в КСК (Фойгт)                      │
│  tc[i]       — критические напряжения (эволюция)             │
│  sp_shifts   — скорости сдвига γ̇^α                          │
└─────────────────────────────────────────────────────────────┘
                              │
                              │ множество зёрен (OpenMP)
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                          RVE                                 │
├─────────────────────────────────────────────────────────────┤
│  stress_av   — средние напряжения                            │
│  W_av        — средний спин                                  │
│  Dp_av       — средняя скорость пласт. деформации            │
│  F_av        — средний градиент деформации                   │
│  acc_strain  — накопленная эквивалентная деформация          │
└─────────────────────────────────────────────────────────────┘
```

## Быстрый старт

### Компиляция

```
Открой EigenProject.sln в Visual Studio 2017
Build → Build Solution (Release x64)
```

Зависимости: **Eigen 3.4+** (путь в проекте: `D:\eigen-3.4.0`), **doctest** (`.h` в корне проекта)

### Запуск основного проекта

```
EigenProject.exe          # работа с моделью (EigenProject.cpp)
```

### Запуск валидации текстур

`rotation_validation.cpp` компилируется отдельно (не входит в .sln):

```bash
# Из папки проекта:
cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 rotation_validation.cpp Grain.cpp RVE.cpp MaterialParams.cpp

rotation_validation.exe           # генерация текстур (4 режима нагружения)
```

### Автоматизированные тесты

`tests.cpp` — автономные тесты с фреймворком [doctest](https://github.com/doctest/doctest) (`.h` в корне проекта):

```bash
# Из папки проекта:
cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 /D_ENABLE_EXTENDED_ALIGNED_STORAGE tests.cpp Grain.cpp RVE.cpp MaterialParams.cpp

tests.exe                         # запуск всех тестов
tests.exe -s                      # с выводом успешных проверок
tests.exe --test-case="no rotation"  # фильтр по имени теста
```

Тесты: проверка поворота зёрен, базовой функциональности RVE, верификация кривых напряжение-деформация, неявная схема с консистентным тангенсом, multi-element UMAT, thread-safety (18 тестов, 256 checks, допуск 5%, два [slow]).

### Построение полюсных фигур (MATLAB + MTEX)

```matlab
cd 'rotation_validation_test_results'
plot_polefigures   # построит PF для всех 4 режимов и сравнение
```

## Формат входного файла материала

```
hard_mode rot_mode
num_slip_sys
b1_x b1_y b1_z n1_x n1_y n1_z    (система 1: b=направление сдвига, n=нормаль)
...
c11 c12 c44                        (упругие константы, Па)
[если hard_mode==2: A q h0 ts alfa]  (параметры Kocks-Mecking)
[если hard_mode==1: a b]             (параметры изотропного упрочнения)
tc0                                (начальное критическое напряжение, Па)
strain_rate_infl                   (γ̇_0 для закона Хатчинсона)
n_power                            (степень скоростного упрочнения)
[integ_mode]```

### Пример (FCC, grain.txt)

```
2 2
12
0 1 1 1 -1 1
0 1 1 1 1 -1
...
107300000000 60080000000 28300000000
1 2 115e+06 34e+06 2.25
6e+06
1e-2
60
```

## Режимы нагружения

### Жёсткое (Simulate_Step)

Задан градиент скорости `L`:

```cpp
RVE rve("rve_validation.txt");

double eps_dot = 0.1;
Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
L(0,0) =  0.5 * eps_dot;   // поперечное растяжение
L(1,1) =  0.5 * eps_dot;
L(2,2) = -1.0 * eps_dot;   // сжатие вдоль Z

rve.Simulate_Step(dt, L);
```

### Мягкое (Simulate_Step_Soft)

Задана скорость изменения напряжений `dS_dt`:

```cpp
Eigen::Matrix3d dS_dt = Eigen::Matrix3d::Zero();
dS_dt(0,0) = 1e8;  // Па/с

rve.Simulate_Step_Soft(dt, dS_dt);
```
(Мягкое нагружение ещё в разработке)

## Структура проекта

```
EigenProject/
├── tensor_utils.h                — общие функции тензорной алгебры (Voigt, 9d)
├── MaterialParams.h/cpp          — материальные параметры (shared_ptr)
├── Grain.h/cpp                   — класс зерна
├── RVE.h/cpp                     — представительский объём (RVE)
├── EigenProject.cpp              — main() проекта (работа с моделью)
├── grain.txt                     — FCC материал (12 систем {111}<110>)
├── rve_validation.txt            — конфиг RVE (1000 зёрен)
├── grain_test_iso.txt            — материал с hard_mode=1 (изотропное упрочнение)
├── rve_validation_iso.txt        — конфиг RVE для изотропного упрочнения
├── dochtet.h                    — фреймворк для автотестов
├── grain_test_shear.txt          — материал для теста поворота (rot_mode=1)
├── grain_implicit.txt             — материал с integ_mode=1 (неявная схема, hard_mode=2)
├── rve_validation_impl.txt        — конфиг RVE для неявной схемы
├── grain_test_iso_impl.txt        — изотропное упрочнение + integ_mode=1
├── tests.cpp                     — автоматизированные тесты (7 тестов)
│
├── cp_model/                     — header-only библиотека (thread-safe)
│   ├── tensor_utils.hpp          — Voigt ↔ tensor conversion
│   ├── MaterialParams.hpp        — структура материала
│   └── Grain.hpp                 — класс зерна (интеграция)
│
├── cp_umat/                      — Abaqus/Standard UMAT
│   ├── cp_umat.cpp               — UMAT entry (thread-safe, __declspec(thread))
│   ├── cp_umat_interface.hpp     — CPUMAT class
│   ├── cp_umat_types.hpp         — A_to_V / V_to_A conversion
│   └── build_cp_umat.bat         — сборка (MSVC, /GR- /EHsc-)
│
├── cp_umat_validation/           — тесты и .inp файлы для Abaqus
│   ├── tension_8el_random.inp    — 8-элементный тест
│   ├── tension_343el_random.inp  — 343-элементный тест (7×7×7)
│   └── extract_s33.py            — извлечение S33 из ODB
│
├── rotation_validation_test_results/
│   ├── rotation_validation.cpp   — валидация: unit-тесты + текстуры
│   ├── plot_polefigures.m        — MTEX скрипт для полюсных фигур
│   ├── orientations_*.txt        — ориентации на 30/50/100%
│   └── PF_*.png                  — полюсные фигуры
│
├── CONTEXT.md                    — техническая документация
└── README.md                     — этот файл
```

## Математическая постановка

### Мезоуровень — одиночное зерно (класс Grain)

#### 1. Кинематика

Градиент скорости раскладывается на симметричную и антисимметричную части:

$$\mathbf{L} = \mathbf{D} + \mathbf{W}, \qquad \mathbf{D} = \frac{1}{2}(\mathbf{L} + \mathbf{L}^T), \qquad \mathbf{W} = \frac{1}{2}(\mathbf{L} - \mathbf{L}^T)$$

Градиент деформации обновляется: $\mathbf{F}(t + dt) = \exp(dt \cdot \mathbf{L}) \cdot \mathbf{F}(t)$

Для малых $||dt \cdot \mathbf{L}|| < 10^{-2}$ — быстрое 3-членное разложение Тейлора $\exp(\mathbf{A}) \approx \mathbf{I} + \mathbf{A} + \mathbf{A}^2/2!$. Для крупных шагов — полная `expm()` из Eigen.

#### 2. Системы скольжения

Каждая система скольжения α задаётся единичным направлением сдвига $\mathbf{\bar{b}}^\alpha$ и нормалью к плоскости скольжения $\mathbf{\bar{n}}^\alpha$ (в КСК). Тензор ориентации (Шмида):

$$\mathbf{M}^\alpha = \mathbf{\bar{b}}^\alpha \otimes \mathbf{\bar{n}}^\alpha$$

Симметричная и антисимметричная части:

$$\mathbf{M}_s^\alpha = \frac{1}{2}(\mathbf{M}^\alpha + \mathbf{M}^{\alpha T}), \qquad \mathbf{W}_{pm}^\alpha = \frac{1}{2}(\mathbf{M}^\alpha - \mathbf{M}^{\alpha T})$$

Системы удваиваются: для каждого $\mathbf{\bar{b}}^\alpha$ добавляется $-\mathbf{\bar{b}}^\alpha$.

#### 3. Упругость (закон Гука в КСК)

$$\dot{\boldsymbol{\sigma}} = \mathbf{P}{:}(\mathbf{D} - \mathbf{D}^p)$$

где $\mathbf{P}$ — тензор упругих модулей кубического кристалла: $P_{iiii} = C_{11}$, $P_{iijj} = C_{12}$, $P_{ijij} = C_{44}$ ($i \neq j$).

В нотации Фойгта используется 6-мерное представление: напряжения $\boldsymbol{\sigma}_v = [\sigma_{11}, \sigma_{22}, \sigma_{33}, \sigma_{23}, \sigma_{13}, \sigma_{12}]$, деформации — в инженерной нотации (сдвиговые компоненты удвоены).

#### 4. Вязкопластический закон (Хатчинсон)

Приведённое напряжение на системе α: $\tau^\alpha = \boldsymbol{\sigma}{:}\mathbf{M}_s^\alpha$

Скорость сдвига:

$$\dot{\gamma}^\alpha = \dot{\gamma}_0 \cdot D_{eq} \cdot \left|\frac{\tau^\alpha}{\tau_c^\alpha}\right|^n$$

где $D_{eq} = \sqrt{\frac{2}{3} \mathrm{dev}(\mathbf{D}){:}\mathrm{dev}(\mathbf{D})}$ — эквивалентная скорость деформации, $n$ — показатель скоростного упрочнения.

Скорость пластической деформации:

$$\mathbf{D}^p = \sum_{\alpha} \dot{\gamma}^\alpha \mathbf{M}_s^\alpha$$

Активные системы скольжения: $\tau^\alpha > \tau_c^\alpha$.

#### 5. Упрочнение

**Без упрочнения (hard_mode = 0):** $\tau_c^\alpha(t) = \tau_{c0}$

**Изотропное (hard_mode = 1):**

$$\tau_c(t) = \tau_{c0} + a \cdot \left(\sum_{\alpha} |\gamma^\alpha|\right)^b, \quad \forall \alpha$$

**Kocks-Mecking (hard_mode = 2):**

$$\dot{\tau}_c^\alpha = \sum_{\beta} H_{\alpha\beta} \cdot |\dot{\gamma}^\beta|$$

$$H_{ii} = A \cdot \theta_0 \cdot \left|1 - \frac{\tau_c^i}{\tau_s}\right|^\alpha, \qquad H_{ij} = q \cdot H_{ii} \quad (i \neq j)$$

где $\tau_s$ — напряжение насыщения, $\theta_0 = h_0 / \tau_{c0}$.

#### 6. Ротация кристаллической решётки

Ориентация решётки описывается кватернионом $\mathbf{q}$: вращение из ЛСК в КСК — $\mathbf{R} = \mathbf{R}(\mathbf{q})$.

Спин решётки $\boldsymbol{\Omega}$:

* **Материальный поворот (rot_mode = 1, Mandel):** $\boldsymbol{\Omega} = \mathbf{W}$
* **Стеснённый поворот (rot_mode = 2, Taylor):** $\boldsymbol{\Omega} = \mathbf{W} - \mathbf{W}^p$

Пластический спин:

$$\mathbf{W}^p = \sum_{\alpha} \dot{\gamma}^\alpha \mathbf{W}_{pm}^\alpha$$

Обновление кватерниона ($\boldsymbol{\omega}$ в КСК):

$$\boldsymbol{\omega} = \frac{1}{2}[\Omega_{32}-\Omega_{23},\; \Omega_{13}-\Omega_{31},\; \Omega_{21}-\Omega_{12}]^T$$

$$\mathbf{q}(t+dt) = \mathbf{q}(t) \cdot \delta\mathbf{q}(\boldsymbol{\omega} \cdot dt)$$

#### 7. Упругая проверка и адаптивное подшагивание

На каждом шаге проверяется критерий Треска: ${}^*\tau^\alpha \leq \tau_c^\alpha$. Если какая-либо система выходит за предел — определяется момент входа в пластичность и шаг дробится на упругую и пластическую части.

Подшагивание по деформации: если $D_{eq} \cdot dt > \Delta s_{crit}$, шаг разбивается на $n_{sub} = \lceil D_{eq} \cdot dt / \Delta s_{crit} \rceil$ подшагов.

---

### Макроуровень — RVE (класс RVE)

#### 1. Гипотеза Тейлора — Фойгта

Все зёрна испытывают одинаковый градиент скорости (жёсткое нагружение):

$$\mathbf{L}_g = \mathbf{L}_{RVE} \quad \forall g$$

#### 2. Гомогенизация (правило смесей)

Средние напряжения:

$$\bar{\boldsymbol{\sigma}} = \sum_{g=1}^{N} w_g \cdot \boldsymbol{\sigma}_g(\mathbf{L}), \qquad w_g = \frac{V_g}{V_{total}}$$

Средняя скорость пластической деформации:

$$\bar{\mathbf{D}}^p = \sum_{g=1}^{N} w_g \cdot \mathbf{D}^p_g$$

Накопленная эквивалентная деформация:

$$\varepsilon_{acc} = \int_0^t D_{eq}(\tau) \, d\tau, \qquad D_{eq} = \sqrt{\frac{2}{3} \mathrm{dev}(\mathbf{D}){:}\mathrm{dev}(\mathbf{D})}$$

#### 3. Мягкое нагружение (Simulate_Step_Soft)

Задана скорость изменения напряжений $\dot{\bar{\boldsymbol{\sigma}}}$. Вычисляется $\mathbf{D}$ из:

$$\mathbf{D} = \bar{\mathbf{P}}^{-1}{:}\dot{\bar{\boldsymbol{\sigma}}} + \bar{\mathbf{D}}^p_{\text{prev}}$$

$$\mathbf{L} = \mathbf{D} + \bar{\mathbf{W}}_{\text{prev}}$$

где $\bar{\mathbf{P}}$ — осреднённая упругость RVE в ЛСК:

$$\bar{\mathbf{P}} = \sum_{g=1}^{N} w_g \cdot \mathbf{P}_g^{lab}$$

#### 4. Интеграция для Abaqus VUMAT

*В разработке.* Версия для Abaqus/Explicit VUMAT будет опубликована в ветке `vumat_implementation`.

## Исправленные баги

| # | Баг | Влияние |
|---|---|---|
| 1 | **Wpm: перепутаны индексы** — `W+Wp` вместо `W-Wp` | Неверная ротация, неправильный тип текстуры |
| 2 | Порядок кватернионов `dq*q` вместо `q*dq` | Неверное направление вращения |
| 3 | `sp_shifts` не занулялись при упругом шаге | Мусор в пластическом спине |
| 4 | `check_elastic` не возвращал значение | Undefined behavior |
| 5 | `orient_q` не инициализирован в `init_common()` | Случайный поворот зерна при старте |

Подробности в [CONTEXT.md](CONTEXT.md).

## Зависимости

- [Eigen](https://eigen.tuxfamily.org/) 3.4+
- [doctest](https://github.com/doctest/doctest) (`.h` в корне, для автотестов)
- C++14
- OpenMP (опционально)
- MATLAB + MTEX (для валидации текстур)

## Лицензия

Apache License, Version 2.0

## Автор

Pavel Gladkikh
