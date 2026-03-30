# RVE и MaterialParams Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Добавить класс MaterialParams для хранения параметров материала и класс RVE для моделирования представительного объёма с множеством зёрен.

**Architecture:** 
- `MaterialParams` — immutable структура с тензорами M, Wpm, P_v и параметрами упрочнения
- `Grain` использует `shared_ptr<const MaterialParams>` вместо локальных копий
- `RVE` владеет MaterialParams и вектором зёрен, вычисляет осреднённые характеристики

**Tech Stack:** C++14, Eigen 3.4+, `<memory>` для shared_ptr

---

## File Structure

```
EigenProject/
├── MaterialParams.h         — структура параметров материала (NEW)
├── MaterialParams.cpp        — реализация парсинга из файла (NEW)
├── Grain.h                   — модифицированный класс зерна
├── Grain.cpp                 — модифицированная реализация
├── RVE.h                     — класс представительного объёма (NEW)
├── RVE.cpp                   — реализация RVE (NEW)
├── rve_test.txt              — тестовый конфиг для RVE (NEW)
└── EigenProject.cpp          — тесты
```

---

### Task 1: Create MaterialParams struct

**Files:**
- Create: `MaterialParams.h`
- Create: `MaterialParams.cpp`

- [ ] **Step 1: Write MaterialParams.h**

```cpp
#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <memory>

struct MaterialParams {
    using Vector6d = Eigen::Matrix<double, 6, 1>;
    using Matrix6d = Eigen::Matrix<double, 6, 6>;

    // Системы скольжения (удвоенные для ±b)
    std::vector<Eigen::Matrix3d> M;        // Тензоры ориентации M = b⊗n
    std::vector<Eigen::Matrix3d> Wpm;      // Антисимметричная часть (пластический спин)
    std::vector<Vector6d> Ms_v;            // Симметричная часть в Фойгте
    int num_slyde_sys;                     // Число систем (до удвоения)

    // Упругость
    Matrix6d P_v;                          // Матрица упругости в Фойгте

    // Пластичность
    double tc0;                            // Начальное критическое напряжение
    double strain_rate_infl;               // Параметр γ̇_0

    // Упрочнение
    unsigned hard_mode;
    Eigen::MatrixXd H;                     // Матрица упрочнения
    // Параметры Kocks-Mecking (hard_mode=2)
    double A = 0, q = 0, h0 = 0, ts = 1e30, alfa = 0;

    // Конструктор из файла
    explicit MaterialParams(const std::string& filename);
};

using MaterialParamsPtr = std::shared_ptr<const MaterialParams>;
```

- [ ] **Step 2: Write MaterialParams.cpp**

Извлечь логику парсинга из Grain.cpp в отдельный файл. Код для чтения систем скольжения, упругих констант и параметров упрочнения.

```cpp
#include "MaterialParams.h"
#include <iostream>
#include <fstream>
#include <cmath>

static int tensorToVoigtIndex(int i, int j) {
    static const int voigt_index[3][3] = {
        {0, 5, 4},
        {5, 1, 3},
        {4, 3, 2}
    };
    return voigt_index[i][j];
}

static std::pair<int, int> voigtToTensorIndex(int voigt_idx) {
    static const std::pair<int, int> tensor_index[6] = {
        {0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}
    };
    return tensor_index[voigt_idx];
}

static Eigen::Matrix<double, 6, 1> tensorToVoigt(const Eigen::Matrix3d& A) {
    Eigen::Matrix<double, 6, 1> v;
    v.setZero();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            v(tensorToVoigtIndex(i, j)) += A(i, j);
    return v;
}

MaterialParams::MaterialParams(const std::string& filename) {
    std::ifstream file(filename);
    
    unsigned rot_mode_temp;
    file >> hard_mode >> rot_mode_temp;  // rot_mode игнорируется для материала
    file >> num_slyde_sys;

    std::vector<double> b(3), n(3);
    M.resize(num_slyde_sys * 2);
    Wpm.resize(num_slyde_sys * 2);
    Ms_v.resize(num_slyde_sys * 2);

    for (int i = 0; i < num_slyde_sys; i++) {
        for (int j = 0; j < 3; j++) file >> b[j];
        for (int j = 0; j < 3; j++) file >> n[j];

        // Проверка ортогональности b и n
        if (std::abs(n[0]*b[0] + n[1]*b[1] + n[2]*b[2]) > 1e-10)
            std::cout << "Warning: b and n not orthogonal for system " << i << std::endl;

        // Нормировка
        double norm_n = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        double norm_b = std::sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
        for (int j = 0; j < 3; j++) { n[j] /= norm_n; b[j] /= norm_b; }

        // +b система
        M[i*2].setZero();
        Wpm[i*2].setZero();
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) {
                M[i*2](j, k) = b[k] * n[j];
                Wpm[i*2](j, k) = 0.5 * (b[k]*n[j] - b[j]*n[k]);
            }
        Eigen::Matrix3d Ms = 0.5 * (M[i*2] + M[i*2].transpose());
        Ms_v[i*2] = tensorToVoigt(Ms);

        // -b система
        for (int j = 0; j < 3; j++) b[j] = -b[j];
        M[i*2+1].setZero();
        Wpm[i*2+1].setZero();
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) {
                M[i*2+1](j, k) = b[k] * n[j];
                Wpm[i*2+1](j, k) = 0.5 * (b[k]*n[j] - b[j]*n[k]);
            }
        Ms = 0.5 * (M[i*2+1] + M[i*2+1].transpose());
        Ms_v[i*2+1] = tensorToVoigt(Ms);
    }

    // Упругие константы
    double c1, c2, c3;
    file >> c1 >> c2 >> c3;

    P_v.setZero();
    P_v(0,0) = P_v(1,1) = P_v(2,2) = c1;
    P_v(0,1) = P_v(1,0) = P_v(1,2) = P_v(2,1) = P_v(0,2) = P_v(2,0) = c2;
    P_v(3,3) = P_v(4,4) = P_v(5,5) = c3;

    // Параметры упрочнения
    H.resize(num_slyde_sys * 2, num_slyde_sys * 2);
    H.setZero();

    if (hard_mode == 2) {
        file >> A >> q >> h0 >> ts >> alfa;
        for (int i = 0; i < num_slyde_sys * 2; i++)
            for (int j = 0; j < num_slyde_sys * 2; j++)
                H(j, i) = (i == j) ? A : A * q;
    }

    file >> tc0;
    file >> strain_rate_infl;
}
```

- [ ] **Step 3: Verify compilation**

Run: build command in VS or verify syntax
Expected: MaterialParams.h and MaterialParams.cpp compile without errors

---

### Task 2: Modify Grain class to use MaterialParams

**Files:**
- Modify: `Grain.h`
- Modify: `Grain.cpp`

- [ ] **Step 1: Modify Grain.h**

Заменить локальные копии M, Wpm, P_v, tc0, strain_rate_infl, hard_mode, H, A, q, h0, ts, alfa на указатель на MaterialParams. Добавить новый конструктор.

```cpp
#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/CXX11/Tensor>
#include <string>
#include <memory>
#include "MaterialParams.h"

class Grain
{
    using Vector6d = Eigen::Matrix<double, 6, 1>;
    using Matrix6d = Eigen::Matrix<double, 6, 6>;
    using Vector9d = Eigen::Matrix<double, 9, 1>;

    // Параметры материала (общие для всех зёрен этого материала)
    MaterialParamsPtr params;
    
    // Индивидуальное состояние зерна
    unsigned rot_mode;
    unsigned hard_mode;  // Кэшируется из params для быстрого доступа

    Eigen::Matrix3d L;
    Vector6d D_v;
    Vector6d Dp_v;
    Vector6d Stress_v;
    Eigen::Matrix3d F;

    // Ориентация зерна: кватернион КСК → ЛСК
    Eigen::Quaterniond orient_q;
    mutable Eigen::Matrix3d R_buf;

    // Критические напряжения (индивидуальные для зерна)
    std::vector<double> tc;
    std::vector<double> sp_shifts;
    std::vector<bool> last_step;
    std::vector<bool> new_step;
    
    // Матрица упрочнения (копия, т.к. обновляется в integrate)
    Eigen::MatrixXd H;

    // Буферы
    mutable Eigen::Vector3d omega_buf;
    mutable Eigen::Quaterniond dq_buf;
    mutable Eigen::Matrix3d W_mat_buf;
    mutable Eigen::Matrix3d Omega_buf;
    mutable Eigen::Matrix3d Wp_mat_buf;
    mutable Eigen::Matrix3d S_buf;
    Eigen::Matrix3d Buf;
    Vector6d Buf_v;

    double dT;
    double Sum_time;
    double acc_def, time, stress_inten;

    void clear();
    bool check_elastic(double dt);
    void find_SS(double dT);
    void integrate(double dt);
    void computePlasticSpin();
    void updateOrientation(const Eigen::Matrix3d& Omega, double dt);

public:
    // Конструктор для автономного использования (создаёт MaterialParams)
    explicit Grain(const std::string& filename);
    
    // Конструктор для использования в RVE
    Grain(MaterialParamsPtr p, const Eigen::Quaterniond& orientation, unsigned rot_mode);

    void Simulate_Step(double dt, Eigen::Matrix3d Lnew);

    Vector6d get_stress() const;
    Eigen::Quaterniond get_orientation() const;
    void set_orientation(const Eigen::Quaterniond& q);
    Eigen::Matrix3d get_stress_lab() const;
    
    // Доступ к параметрам материала (read-only)
    MaterialParamsPtr get_params() const { return params; }
};
```

- [ ] **Step 2: Modify Grain.cpp constructors**

```cpp
// Конструктор из файла (автономное использование)
Grain::Grain(const std::string& filename) {
    params = std::make_shared<MaterialParams>(filename);
    hard_mode = params->hard_mode;
    rot_mode = 1;  // По умолчанию материальный поворот
    
    // Инициализация индивидуального состояния
    int n = params->num_slyde_sys * 2;
    tc.resize(n, params->tc0);
    sp_shifts.resize(n, 0.0);
    last_step.resize(n, false);
    new_step.resize(n, false);
    H = params->H;  // Копия матрицы упрочнения
    
    // Остальная инициализация
    orient_q = Eigen::Quaterniond::Identity();
    R_buf.setIdentity();
    F.setIdentity();
    L.setZero();
    D_v.setZero();
    Dp_v.setZero();
    Stress_v.setZero();
    Buf.setZero();
    Buf_v.setZero();
    omega_buf.setZero();
    dq_buf = Eigen::Quaterniond::Identity();
    W_mat_buf.setZero();
    Omega_buf.setZero();
    Wp_mat_buf.setZero();
    S_buf.setZero();
    acc_def = time = stress_inten = 0;
}

// Конструктор для RVE
Grain::Grain(MaterialParamsPtr p, const Eigen::Quaterniond& orientation, unsigned rot)
    : params(p)
    , rot_mode(rot)
    , hard_mode(p->hard_mode)
    , H(p->H)
{
    int n = params->num_slyde_sys * 2;
    tc.resize(n, params->tc0);
    sp_shifts.resize(n, 0.0);
    last_step.resize(n, false);
    new_step.resize(n, false);
    
    set_orientation(orientation);
    F.setIdentity();
    L.setZero();
    D_v.setZero();
    Dp_v.setZero();
    Stress_v.setZero();
    Buf.setZero();
    Buf_v.setZero();
    omega_buf.setZero();
    dq_buf = Eigen::Quaterniond::Identity();
    W_mat_buf.setZero();
    Omega_buf.setZero();
    Wp_mat_buf.setZero();
    S_buf.setZero();
    acc_def = time = stress_inten = 0;
}
```

- [ ] **Step 3: Update Grain methods to use params->**

Заменить все обращения к `M`, `Wpm`, `Ms_v`, `P_v`, `num_slyde_sys` на `params->M`, `params->Wpm`, `params->Ms_v`, `params->P_v`, `params->num_slyde_sys`.

Ключевые замены:
- `Ms_v[i]` → `params->Ms_v[i]`
- `P_v` → `params->P_v`
- `Wpm0[i]` → `params->Wpm[i]`
- `num_slyde_sys` → `params->num_slyde_sys`
- `strain_rate_infl` → `params->strain_rate_infl`
- `ts` → `params->ts`

- [ ] **Step 4: Verify existing tests still pass**

Run: `EigenProject.exe`
Expected: All 3 tests pass (Test 1, Test 2, Test 3)

---

### Task 3: Create RVE class

**Files:**
- Create: `RVE.h`
- Create: `RVE.cpp`

- [ ] **Step 1: Write RVE.h**

```cpp
#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <string>
#include <memory>
#include "Grain.h"

class RVE {
    using Vector6d = Eigen::Matrix<double, 6, 1>;
    using Matrix6d = Eigen::Matrix<double, 6, 6>;

    // Материалы (RVE владеет)
    std::vector<MaterialParamsPtr> materials;
    
    // Зёрна
    std::vector<Grain> grains;
    
    // Режим ротации для всех зёрен
    unsigned rot_mode;

    // Осреднённые характеристики в ЛСК
    Eigen::Matrix3d stress_av;       // Осреднённый тензор напряжений
    Eigen::Matrix3d W_av;            // Осреднённый спин
    Eigen::Matrix3d Dp_av;           // Осреднённая пластическая скорость деформации
    Eigen::Matrix3d F_av;            // Осреднённый градиент деформации
    
    // Скалярные параметры
    double time;
    double acc_strain;               // Накопленная интенсивность деформации
    Eigen::Matrix3d A_affine;        // Аффинор

public:
    // Конструктор из файла конфигурации RVE
    explicit RVE(const std::string& config_file);
    
    // Шаг интегрирования (L — градиент скорости в ЛСК, одинаковый для всех зёрен)
    void Simulate_Step(double dt, const Eigen::Matrix3d& L);
    
    // Геттеры для осреднённых характеристик
    Eigen::Matrix3d get_stress() const { return stress_av; }
    Eigen::Matrix3d get_spin() const { return W_av; }
    Eigen::Matrix3d get_plastic_strain_rate() const { return Dp_av; }
    Eigen::Matrix3d get_deformation_gradient() const { return F_av; }
    double get_time() const { return time; }
    double get_accumulated_strain() const { return acc_strain; }
    
    // Доступ к отдельным зёрнам
    size_t num_grains() const { return grains.size(); }
    const Grain& get_grain(size_t i) const { return grains[i]; }
};
```

- [ ] **Step 2: Write RVE.cpp**

```cpp
#include "RVE.h"
#include <iostream>
#include <fstream>
#include <random>

RVE::RVE(const std::string& config_file) {
    std::ifstream file(config_file);
    
    // Читаем rot_mode
    file >> rot_mode;
    
    // Читаем число материалов
    int num_materials;
    file >> num_materials;
    
    // Читаем материалы
    materials.resize(num_materials);
    for (int m = 0; m < num_materials; m++) {
        std::string mat_file;
        int num_grains_of_material;
        file >> mat_file >> num_grains_of_material;
        
        auto mat = std::make_shared<MaterialParams>(mat_file);
        materials[m] = mat;
        
        // Создаём зёрна с этим материалом
        for (int g = 0; g < num_grains_of_material; g++) {
            // Случайная ориентация (равномерное распределение на SO(3))
            Eigen::Quaterniond random_q = Eigen::Quaterniond::UnitRandom();
            grains.emplace_back(mat, random_q, rot_mode);
        }
    }
    
    // Инициализация осреднённых величин
    stress_av.setZero();
    W_av.setZero();
    Dp_av.setZero();
    F_av.setIdentity();
    A_affine.setIdentity();
    time = 0;
    acc_strain = 0;
}

void RVE::Simulate_Step(double dt, const Eigen::Matrix3d& L) {
    // Сброс осреднённых величин
    stress_av.setZero();
    W_av.setZero();
    Dp_av.setZero();
    
    // Шаг для каждого зерна
    for (auto& grain : grains) {
        grain.Simulate_Step(dt, L);
    }
    
    // Осреднение (простое среднее по объёму, вес = 1/N)
    double n = static_cast<double>(grains.size());
    for (const auto& grain : grains) {
        // Напряжения в ЛСК
        stress_av += grain.get_stress_lab() / n;
        
        // TODO: добавить осреднение W, Dp, F
    }
    
    // Обновление времени
    time += dt;
    
    // Обновление F_av (экспонента от L)
    F_av = (dt * L).exp() * F_av;
}
```

- [ ] **Step 3: Create test config file rve_test.txt**

```
1
1
grain_test.txt 10
```

Формат:
- Строка 1: rot_mode
- Строка 2: число материалов
- Для каждого материала: путь к файлу материала, число зёрен

- [ ] **Step 4: Add RVE test to EigenProject.cpp**

```cpp
void test_rve_basic()
{
    std::cout << "=== Test 4: RVE basic functionality ===" << std::endl;

    RVE rve("rve_test.txt");
    
    std::cout << "Number of grains: " << rve.num_grains() << std::endl;
    
    double eps = 0.001;
    Eigen::Matrix3d L;
    L.setZero();
    L(0, 0) = eps;
    L(1, 1) = -0.5 * eps;
    L(2, 2) = -0.5 * eps;

    double dt = 1e-3;
    for (int i = 0; i < 10; i++) {
        rve.Simulate_Step(dt, L);
    }

    Eigen::Matrix3d stress = rve.get_stress();
    std::cout << "Average stress:\n" << stress << std::endl;
    std::cout << "Time: " << rve.get_time() << std::endl;
    
    if (rve.num_grains() == 10) {
        std::cout << "TEST 4 PASSED: RVE created 10 grains" << std::endl;
    } else {
        std::cout << "TEST 4 FAILED: Wrong number of grains" << std::endl;
    }
}
```

- [ ] **Step 5: Run tests**

Run: `EigenProject.exe`
Expected: All 4 tests pass

---

### Task 4: Update project files and documentation

**Files:**
- Modify: `EigenProject.vcxproj` (add new files)
- Modify: `CONTEXT.md`
- Modify: `README.md`

- [ ] **Step 1: Add new files to vcxproj**

Добавить MaterialParams.h, MaterialParams.cpp, RVE.h, RVE.cpp в проект Visual Studio.

- [ ] **Step 2: Update CONTEXT.md**

Добавить секции о MaterialParams и RVE.

- [ ] **Step 3: Commit**

```bash
git add MaterialParams.h MaterialParams.cpp RVE.h RVE.cpp Grain.h Grain.cpp EigenProject.cpp
git commit -m "feat: add MaterialParams and RVE classes"
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | MaterialParams struct | MaterialParams.h/cpp |
| 2 | Modify Grain | Grain.h/cpp |
| 3 | RVE class | RVE.h/cpp, rve_test.txt |
| 4 | Documentation | vcxproj, CONTEXT.md, README.md |
