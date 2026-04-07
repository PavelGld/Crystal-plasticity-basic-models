# Context: Crystal Plasticity Grain Model

## Общее описание проекта

C++ реализация модели кристаллической пластичности для одиночного зерна (кристаллита) и представительского объёма (RVE). Проект использует библиотеку Eigen для линейной алгебры. Все вычисления ведутся в нотации Фойгта (6-мерное пространство для симметричных тензоров).

**Текущий статус:** Базовая упруговязкопластическая модель + ротация кристаллической решётки + RVE + мягкое нагружение + валидация текстур — реализованы и протестированы. Полюсные фигуры для 4 режимов нагружения совпадают с литературными данными.

## Файловая структура

```
EigenProject/
├── MaterialParams.h/cpp          — структура материальных параметров (shared_ptr)
├── Grain.h/cpp                   — класс зерна (ориентация, системы скольжения, интеграция)
├── RVE.h/cpp                     — представительский объём (множество зёрен, осреднение)
├── EigenProject.cpp              — main() проекта (работа с моделью)
├── grain.txt                     — параметры материала (FCC, 12 систем {111}<110>)
├── rve_validation.txt            — конфиг RVE (1000 зёрен)
│
├── rotation_validation_test_results/
│   ├── rotation_validation.cpp   — файл валидации: unit-тесты + генерация текстур
│   ├── plot_polefigures.m        — скрипт MTEX для построения полюсных фигур
│   ├── orientations_*.txt        — выходные ориентации (30%, 50%, 100%)
│   └── PF_*.png                  — полюсные фигуры {111}, {110}, {100}
│
├── CONTEXT.md                    — этот файл (техническая документация)
└── README.md                     — обзор проекта
```

## Структура MaterialParams

`MaterialParams` — неизменяемые материальные параметры, общие для нескольких зёрен одного материала.

```cpp
struct MaterialParams {
    // Системы скольжения
    std::vector<Eigen::Matrix3d> M;      // тензоры ориентации M = b⊗n
    std::vector<Eigen::Matrix3d> Wpm;    // антисимметричная часть Wpm = 0.5*(b⊗n - n⊗b)
    std::vector<Vector6d> Ms_v;          // симметричная часть в Фойгте
    int num_slyde_sys;

    // Упругость
    Matrix6d P_v;                        // матрица упругости в Фойгте

    // Режимы модели
    unsigned hard_mode;                  // 0=изотроп, 1=степенное, 2=Kocks-Mecking
    unsigned rot_mode;                   // 0=нет, 1=материальный, 2=Тейлор

    // Параметры пластичности
    double tc0;                          // начальное критическое напряжение
    double strain_rate_infl;             // γ̇_0 для закона Хатчинсона
    double n_power;                      // степень в законе Хатчинсона

    // Упрочнение (Kocks-Mecking)
    Eigen::MatrixXd H;                   // матрица взаимодействия
    double A, q, h0, ts, alfa;
};
```

**Важно:** Используется `shared_ptr<const MaterialParams>` для разделения параметров между зёрнами.

## Класс RVE

```cpp
class RVE {
    std::vector<MaterialParamsPtr> materials;
    std::vector<Grain> grains;

    // Усреднённые величины
    Eigen::Matrix3d stress_av;           // средние напряжения
    Eigen::Matrix3d W_av;                // средний спин
    Eigen::Matrix3d Dp_av;               // средняя скорость пласт. деформации
    Eigen::Matrix3d F_av;                // средний градиент деформации
    double acc_strain;                   // накопленная эквивалентная деформация

    // Методы
    void Simulate_Step(double dt, Matrix3d L);           // жёсткое нагружение
    void Simulate_Step_Soft(double dt, Matrix3d dS_dt);  // мягкое нагружение
};
```

### Формат файла конфигурации RVE

```
num_materials
material_file.txt  num_grains
```

Зёрна инициализируются случайными ориентациями.

### Мягкое нагружение (Simulate_Step_Soft)

Итеративный алгоритм: задана скорость изменения напряжений `dS_dt`, вычисляется `D` через среднюю упругость RVE.

## Системы координат

### Конвенция
- **ЛСК** — лабораторная система (неподвижная)
- **КСК** — кристаллическая система (подвижная, связана с решёткой)

### Кватернион orient_q
- Описывает поворот **КСК → ЛСК**
- Обновление: `orient_q = orient_q * dq` (умножение справа, omega в КСК)

### Преобразования тензоров
- КСК → ЛСК: `T_lab = R * T_crystal * R^T`
- ЛСК → КСК: `T_crystal = R^T * T_lab * R`

## Нотация Фойгта

```
Voigt:    0    1    2    3      4      5
Tensor:  00   11   22   23     13     12
                        32     31     21
```

Сдвиговые: `v(3) = A_23+A_32`, `v(4) = A_13+A_31`, `v(5) = A_12+A_21`

## Алгоритм Simulate_Step

```
1. clear() — last_step = new_step
2. L → КСК (если rot_mode ≠ 0): L_crystal = R^T * L_lab * R
3. D = 0.5*(L + L^T), D_v = tensorToVoigt(D)
4. Цикл подшагов:
   a. check_elastic(dt) — упругая проверка (Решение по условию ТРЕШ)
   b. find_SS(dt) — пластический шаг (поиск γ̇)
   c. integrate(dt) — упрочнение + ротация
5. F = exp(dt*L) * F
```

## Ротация решётки

### Материальный поворот (rot_mode=1)
```
Ω = W = 0.5 * (L - L^T)
```

### Стеснённый поворот Тейлора (rot_mode=2)
```
W = 0.5 * (L - L^T)
Wp = Σ γ̇^α · Wpm^α       // пластический спин
Ω = W - Wp
```

### Обновление кватерниона
```
omega = 0.5 * [Ω_32-Ω_23, Ω_13-Ω_31, Ω_21-Ω_12]
dq = AngleAxisd(|omega|*dt, omega/|omega|)
orient_q = orient_q * dq   // СПРАВА (omega в КСК)
```

## Исправленные баги

### Критические (влияющие на результат)
1. **Wpm: перепутаны индексы** — `Wpm = 0.5*(b⊗n - n⊗b)` имел swapped indices, из-за чего `Ω = W - Wp` вычислялось как `W + Wp`. **Исправлено.**
2. **Порядок кватернионов** — было `dq * orient_q` (слева), стало `orient_q * dq` (справа).
3. **sp_shifts не занулялись** при упругом шаге.
4. **check_elastic не возвращал значение** — UB.

### Не критические
5. Неверный аргумент `check_elastic(dT + Sum_time)` → `check_elastic(dT)`.
6. `rot_mode` перенесён в `MaterialParams` (единый источник из файла материала).
7. Const-корректность: `mutable` для буферов в const-методах.

## Валидация текстур

4 режима нагружения проверяются через `rotation_validation.cpp`:

| Нагружение | L | Ожидаемая текстура FCC |
|---|---|---|
| Simple Shear | L₁₂ = γ̇ | Компоненты простого сдвига |
| Uniaxial Tension | diag(-ε̇/2, -ε̇/2, ε̇) | `<111>||ND` |
| Uniaxial Compression | diag(ε̇/2, ε̇/2, -ε̇) | `<110>||ND` |
| Plane Strain | diag(0, -ε̇, ε̇) | Brass-type `{110}<112>` |

Результаты: все 4 типа текстур совпали с литературными данными (проверка через MTEX pole figures).

## Формат входного файла материала

```
hard_mode rot_mode
num_slyde_sys
b1_x b1_y b1_z n1_x n1_y n1_z
...
c11 c12 c44
[A q h0 ts alfa]   (если hard_mode==2)
tc0
strain_rate_infl
n_power
```

Системы скольжения удваиваются автоматически (±b).

## Оптимизации

- Все промежуточные Eigen-объекты — буферы-члены класса (без аллокаций)
- `.noalias()` везде где возможно
- OpenMP для параллелизации зёрен в RVE
- Кэш матрицы поворота Фойгта (T_voigt_valid)

## Методы класса Grain

```cpp
Grain(std::string info)                    // конструктор из файла
Grain(MaterialParamsPtr, Quaterniond q)    // конструктор для RVE
void Simulate_Step(double dt, Matrix3d L)  // шаг интегрирования

Vector6d get_stress()                      // напряжения в КСК
Matrix3d get_stress_lab()                  // напряжения в ЛСК
Quaterniond get_orientation()              // кватернион ориентации
void set_orientation(const Quaterniond& q) // установка ориентации
Matrix6d get_elasticity_lab()              // упругость в ЛСК (Фойгт)
```

## Методы класса RVE

```cpp
RVE(const std::string& config_file)
void Simulate_Step(double dt, Matrix3d L)
void Simulate_Step_Soft(double dt, Matrix3d dS_dt)

Matrix3d get_stress()
Matrix3d get_spin()
Matrix3d get_plastic_strain_rate()
Matrix3d get_deformation_gradient()
double get_time()
double get_accumulated_strain()
size_t num_grains()
const Grain& get_grain(size_t i)
```
