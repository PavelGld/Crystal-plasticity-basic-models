# Implicit Integration (Newton–Raphson) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Добавить в `Grain` неявную схему интегрирования (Newton–Raphson по напряжениям, 6×6) с консистентным тангенсом, сохранив явную схему как регрессионный эталон (переключение через `integ_mode` в файле материала).

**Architecture:** Невязка R(σ) = σ − σ_tr + Δt·P:Dᵖ(σ) = 0 решается Ньютоном с line search Армихо и fallback-делением шага пополам. Упрочнение — fully implicit (связка ∂Δγ/∂σ с предыдущей итерации). Тангенс ∂Δσ/∂Δε = J⁻¹·P на сошедшемся состоянии; RVE-тангенс — осреднение зернных в ЛСК. Спека: `docs/superpowers/specs/2026-07-26-implicit-integration-design.md`.

**Tech Stack:** C++14, Eigen 3.4 (`D:\eigen-3.4.0`), doctest (`doctest.h` в корне), OpenMP, MSVC `cl`.

## Global Constraints

- Ветка: `implicit-integration` (уже создана, спека закоммичена).
- Сборка тестов (из корня проекта, VS Developer prompt):
  `cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 /D_ENABLE_EXTENDED_ALIGNED_STORAGE tests.cpp Grain.cpp RVE.cpp MaterialParams.cpp`
- Запуск всех тестов: `tests.exe`
- Запуск без медленных: `tests.exe -tce="*[slow]*"`
- Явная схема (`check_elastic`/`find_SS`) не изменяется по поведению; `Simulate_Step_Soft` не трогаем (известно нерабочий, вне скоупа).
- Обновление F в Grain через `exp(dt·L_КСК)` — осознанная конвенция, не менять.
- Стиль кода: буферы — члены класса (без аллокаций в горячем цикле), `.noalias()` где возможно, `static constexpr` для констант.
- Коммиты на русском/английском в стиле репозитория (`feat:`, `test:`, `refactor:`, `docs:`).

### Отступления от спеки (согласовано на уровне плана)

1. **Юнит-тест якобиана (спека, п.5) слит с тестом тангенса (п.4):** FD-проверка тангенса верифицирует всю цепочку производных end-to-end (R, J, рекурсия G, порог); отдельный тест J потребовал бы раскрытия private-членов ради маргинальной диагностики.
2. **`DS_CRITICAL_IMPL = 1e-4`** — отдельный порог подшагивания для неявной схемы. Без него внешнее подшагивание (`DS_CRITICAL = 1.7e-6`) режет любой большой шаг и тест стабильности на dt×50 теряет смысл.
3. **`NEWTON_TOL = 1e-8`** (относительный, к ‖σ_tr‖/τc) вместо 1e-10 из спеки — практичнее по итерациям; константа, легко подкрутить.
4. Эквивалентность неявной схемы на полных кривых RVE (5%) остаётся, но помечена тегом `[slow]` (долгие прогоны); всегда активны быстрые тесты большого шага (×10, ×50, допуск 10%).
5. **Чистый степенной закон без порога в неявной схеме** (решение пользователя от 2026-07-26, после обнаружения безкорневости невязки на пороге при исполнении Task 3): γ̇ = γ̇₀·D_eq·(τ/τc)^n для всех τ > 0. Явная схема сохраняет порог. Детали — в спеке §2.1.

---

### Task 1: `integ_mode` и `MsE_v` в MaterialParams + конфиги

**Files:**
- Modify: `MaterialParams.h`, `MaterialParams.cpp`
- Create: `grain_implicit.txt`, `grain_test_iso_impl.txt`, `rve_validation_impl.txt`
- Test: `tests.cpp`

**Interfaces:**
- Produces: `MaterialParams::integ_mode` (unsigned, 0=явная, 1=неявная; default 0); `MaterialParams::MsE_v` (`std::vector<Vector6d>`, симметричная часть в **инженерной** нотации — сдвиги ×2; используется в Dᵖ и якобиане неявной схемы).

- [ ] **Step 1: Baseline — собрать и прогнать текущие тесты**

```
cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 /D_ENABLE_EXTENDED_ALIGNED_STORAGE tests.cpp Grain.cpp RVE.cpp MaterialParams.cpp
tests.exe
```
Expected: все 7 тестов PASS (это точка отсчёта).

- [ ] **Step 2: Failing test — парсинг integ_mode**

Добавить в конец `tests.cpp`:

```cpp
TEST_CASE("integ_mode parsing") {
    MaterialParams p_exp("grain.txt");
    CHECK(p_exp.integ_mode == 0);

    MaterialParams p_imp("grain_implicit.txt");
    CHECK(p_imp.integ_mode == 1);

    REQUIRE(p_imp.MsE_v.size() == p_imp.Ms_v.size());
    for (size_t i = 0; i < p_imp.Ms_v.size(); i++) {
        CHECK(p_imp.MsE_v[i](0) == doctest::Approx(p_imp.Ms_v[i](0)));
        CHECK(p_imp.MsE_v[i](3) == doctest::Approx(2.0 * p_imp.Ms_v[i](3)));
        CHECK(p_imp.MsE_v[i](4) == doctest::Approx(2.0 * p_imp.Ms_v[i](4)));
        CHECK(p_imp.MsE_v[i](5) == doctest::Approx(2.0 * p_imp.Ms_v[i](5)));
    }
}
```

Run: `tests.exe --test-case="integ_mode parsing"`
Expected: FAIL (компиляция: нет членов `integ_mode`, `MsE_v`; либо рантайм-ошибка открытия `grain_implicit.txt`).

- [ ] **Step 3: Создать конфиги**

`grain_implicit.txt` — копия `grain.txt` + строка `1` в конце:

```
2 2
12
0 1 1 1 -1 1
0 1 1 1 1 -1
0 -1 1 1 1 1
0 -1 1 -1 1 1
-1 1 0 1 1 1
-1 1 0 1 1 -1
1 1 0 -1 1 1
1 1 0 1 -1 1
1 0 1 1 1 -1
1 0 1 -1 1 1
-1 0 1 1 1  1
-1 0 1 1 -1 1
106750000000 60410000000 28340000000
1 2 115e+06 34e+06 2.25
6e+06
0.588235
50
1
```

`grain_test_iso_impl.txt`:

```
1 2
12
0 1 1 1 -1 1
0 1 1 1 1 -1
0 -1 1 1 1 1
0 -1 1 -1 1 1
-1 1 0 1 1 1
-1 1 0 1 1 -1
1 1 0 -1 1 1
1 1 0 1 -1 1
1 0 1 1 1 -1
1 0 1 -1 1 1
-1 0 1 1 1  1
-1 0 1 1 -1 1
106750000000 60410000000 28340000000
29e6 0.4
3e6
0.588235
50
1
```

`rve_validation_impl.txt`:

```
1
grain_implicit.txt 1000
```

- [ ] **Step 4: Реализация**

В `MaterialParams.h` внутри `struct MaterialParams` после строки `std::vector<Vector6d> Ms_v;` добавить:

```cpp
    std::vector<Vector6d> MsE_v;   // симметричная часть в инженерной нотации (сдвиги ×2)
```

и после `unsigned rot_mode;` добавить:

```cpp
    unsigned integ_mode = 0;       // 0=явная схема, 1=неявная (Newton-Raphson)
```

В `MaterialParams.cpp`:
- после `Ms_v.resize(num_slip_sys * 2);` добавить `MsE_v.resize(num_slip_sys * 2);`
- после блока заполнения `Ms_v[i*2](5) = Ms(0, 1);` (система +b) добавить:

```cpp
        MsE_v[i*2] = Ms_v[i*2];
        MsE_v[i*2](3) *= 2.0; MsE_v[i*2](4) *= 2.0; MsE_v[i*2](5) *= 2.0;
```

- аналогично после блока `Ms_v[i*2+1](5) = Ms(0, 1);` (система −b) добавить:

```cpp
        MsE_v[i*2+1] = Ms_v[i*2+1];
        MsE_v[i*2+1](3) *= 2.0; MsE_v[i*2+1](4) *= 2.0; MsE_v[i*2+1](5) *= 2.0;
```

- в конце конструктора после `file >> n_power;` добавить:

```cpp
    unsigned im;
    if (file >> im) integ_mode = im;   // опциональный токен; EOF → остаётся 0
```

- [ ] **Step 5: Тест проходит + регрессия**

```
tests.exe
```
Expected: 8 тестов PASS (7 старых + новый).

- [ ] **Step 6: Commit**

```bash
git add MaterialParams.h MaterialParams.cpp grain_implicit.txt grain_test_iso_impl.txt rve_validation_impl.txt tests.cpp
git commit -m "feat: add integ_mode and engineering-Voigt MsE_v to MaterialParams"
```

---

### Task 2: Рефакторинг `integrate()` — выделить `apply_rotation()` и `update_H_KM()`

Чистая подготовка под неявную схему: блоки ротации и пересчёта H нужны обеим схемам. Поведение явной схемы не меняется.

**Files:**
- Modify: `Grain.h`, `Grain.cpp`
- Test: `tests.cpp` (существующие тесты — регрессия)

**Interfaces:**
- Produces: `void Grain::apply_rotation(double dt)` (private); `void Grain::update_H_KM()` (private). Task 3 вызывает оба из `solve_implicit`.

- [ ] **Step 1: Рефакторинг**

В `Grain.h` в private-секцию (рядом с `void integrate(double dt);`) добавить:

```cpp
    void apply_rotation(double dt);
    void update_H_KM();
```

В `Grain.cpp` добавить новые методы (после `integrate()`):

```cpp
void Grain::update_H_KM()
{
    int n = 2 * params->num_slip_sys;
    if (H_buf.size() != n) H_buf.resize(n);
    for (int j = 0; j < n; j++)
    {
        double t_j = tc[j];
        double Hii = 0;
        if (t_j <= params->ts)
            Hii = params->A * params->h0 * std::pow(std::abs(1.0 - t_j / params->ts), params->alfa);
        H_buf(j) = Hii;
    }
    H.noalias() = params->q * Eigen::VectorXd::Ones(n) * H_buf.transpose();
    H.diagonal() = H_buf;
}

void Grain::apply_rotation(double dt)
{
    if (params->rot_mode == 0) return;

    W_mat_buf = 0.5 * (L - L.transpose());

    if (params->rot_mode == 1)
    {
        Omega_buf = W_mat_buf;
    }
    else if (params->rot_mode == 2)
    {
        computePlasticSpin();
        Omega_buf.noalias() = W_mat_buf - Wp_mat_buf;
    }

    updateOrientation(Omega_buf, dt);
}
```

В `Grain::integrate()` заменить блок Kocks-Mecking (текущие строки `if (params->hard_mode == 2) { ... H.diagonal() = H_buf; }`) на:

```cpp
    if (params->hard_mode == 2)
        update_H_KM();
```

и заменить весь блок ротации (`if (params->rot_mode != 0) { ... updateOrientation(Omega_buf, dt); }`) на:

```cpp
    apply_rotation(dt);
```

- [ ] **Step 2: Регрессия**

```
cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 /D_ENABLE_EXTENDED_ALIGNED_STORAGE tests.cpp Grain.cpp RVE.cpp MaterialParams.cpp
tests.exe
```
Expected: все 8 тестов PASS (поведение не изменилось).

- [ ] **Step 3: Commit**

```bash
git add Grain.h Grain.cpp
git commit -m "refactor: extract apply_rotation and update_H_KM from integrate"
```

---

### Task 3: Ядро Ньютона — `solve_implicit` + роутинг в `Simulate_Step`

**Files:**
- Modify: `Grain.h`, `Grain.cpp`
- Test: `tests.cpp`

**Interfaces:**
- Consumes: `MaterialParams::integ_mode`, `MaterialParams::MsE_v` (Task 1); `apply_rotation`, `update_H_KM` (Task 2).
- Produces: `void Grain::solve_implicit(double dt, int depth)` (private); `void Grain::evaluate_Dp(const Vector6d& sigma, double dt)` (private); `void Grain::assemble_J(double dt)` (private). Константы: `NEWTON_TOL`, `NEWTON_MAX_ITER`, `NEWTON_MAX_BACKTRACKS`, `RATIO_MAX`, `MAX_FALLBACK_DEPTH`, `DS_CRITICAL_IMPL`. После Task 5 сюда добавится заполнение `tangent_buf`.

**Математика (напоминание из спеки):**
- R(σ) = σ − σ_tr + Δt·P:Dᵖ(σ), σ_tr = σₙ + P:Δε
- J = I + Δt·P·Σ_α MsE^α ⊗ G^α, где G^α = ∂γ̇^α/∂σ (строка 1×6)
- γ̇^α = γ̇₀·D_eq·(τ/τc_eff)^n при τ > τc_eff, иначе 0 (порог как в явной схеме)
- τc_eff: fixed-point, 3 sweep'а (1 для hard_mode=0); G — рекурсия с предыдущего sweep'а

- [ ] **Step 1: Failing test — эквивалентность на одиночном зерне**

Добавить в конец `tests.cpp`:

```cpp
TEST_CASE("implicit single grain vs explicit") {
    Grain gr_e("grain.txt");            // явная схема
    Grain gr_i("grain_implicit.txt");   // неявная

    double eps_dot_local = 0.017;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 0) = -0.5 * eps_dot_local;
    L(1, 1) = -0.5 * eps_dot_local;
    L(2, 2) = eps_dot_local;

    double dt_local = 1e-3;
    for (int i = 0; i < 1000; i++) {
        gr_e.Simulate_Step(dt_local, L);
        gr_i.Simulate_Step(dt_local, L);
    }

    double s_e = gr_e.get_stress_intensity();
    double s_i = gr_i.get_stress_intensity();
    CHECK(s_i == doctest::Approx(s_e).epsilon(0.02));
}
```

Run: `tests.exe --test-case="implicit single grain vs explicit"`
Expected на этом шаге: PASS тривиально — до реализации роутинга (Step 4) `grain_implicit.txt` работает по явной схеме, обе кривые совпадают точно. Тест написан заранее и становится настоящей проверкой после Step 4: если неявная схема реализована с ошибкой, расхождение выйдет за 2% и тест упадёт. Гейт корректности — Step 5, не этот прогон.

- [ ] **Step 2: Failing test — изотропное упрочнение под неявной схемой**

Добавить в конец `tests.cpp`:

```cpp
TEST_CASE("implicit isotropic hardening") {
    Grain gr("grain_test_iso_impl.txt");

    double eps_dot_local = 0.017;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 0) = -0.5 * eps_dot_local;
    L(1, 1) = -0.5 * eps_dot_local;
    L(2, 2) = eps_dot_local;

    double tc0 = gr.get_tc(0);
    double dt_local = 1e-3;
    for (int i = 0; i < 1000; i++)
        gr.Simulate_Step(dt_local, L);

    double tc_after = gr.get_tc(0);
    CHECK(tc_after > tc0);

    size_t n = gr.num_slip_systems();
    for (size_t i = 1; i < n; i++)
        CHECK(gr.get_tc(i) == doctest::Approx(tc_after));
}
```

- [ ] **Step 3: Реализация — члены класса**

В `Grain.h`, private-секцию, после строки `static constexpr double DS_CRITICAL = 1.7e-6;` добавить:

```cpp
    // Неявная схема (integ_mode==1)
    static constexpr double NEWTON_TOL = 1e-8;
    static constexpr int NEWTON_MAX_ITER = 20;
    static constexpr int NEWTON_MAX_BACKTRACKS = 5;
    static constexpr double RATIO_MAX = 2.0;          // кламп τ/τc (2^50 ≈ 1e15, переполнения нет)
    static constexpr int MAX_FALLBACK_DEPTH = 8;
    static constexpr double DS_CRITICAL_IMPL = 1e-4;  // порог подшагивания для неявной схемы

    mutable Matrix6d J_buf;
    mutable Matrix6d dyad_buf;
    Vector6d res_v;                 // невязка R
    Vector6d res_new_v;
    Vector6d sigma_trial_v;
    Vector6d sigma_new_v;
    Vector6d delta_sigma_v;
    Vector6d trial_inc_v;           // P_v * D_v (trial-инкремент)
    Eigen::VectorXd gamma_dot;      // скорости сдвига (итерат/сошедшиеся), n
    Eigen::VectorXd tau_buf;        // τ^α, n
    Eigen::VectorXd tc_eff;         // эффективные τc с учётом упрочнения внутри шага, n
    Eigen::MatrixXd G_buf;          // ∂γ̇/∂σ, n×6
    Eigen::MatrixXd G_prev;         // G с предыдущего sweep'а, n×6
    double rate_factor_impl = 0.0;  // γ̇₀ · D_eq на текущий шаг

    void solve_implicit(double dt, int depth);
    void evaluate_Dp(const Vector6d& sigma, double dt);
    void assemble_J(double dt);
```

В `Grain.cpp` в `init_common()` после `sp_shifts.assign(n, 0.0);` добавить:

```cpp
    gamma_dot.setZero(n);
    tau_buf.setZero(n);
    tc_eff.setZero(n);
    G_buf.setZero(n, 6);
    G_prev.setZero(n, 6);
    J_buf.setZero();
    dyad_buf.setZero();
    res_v.setZero();
    res_new_v.setZero();
    sigma_trial_v.setZero();
    sigma_new_v.setZero();
    delta_sigma_v.setZero();
    trial_inc_v.setZero();
```

В `Grain.cpp` добавить методы (после `apply_rotation()`):

```cpp
// Оценка Dp и производных для заданного σ (детерминированная функция σ:
// фиксированное число sweep'ов fixed-point по γ̇↔τc).
// Заполняет: tau_buf, gamma_dot, Dp_v (инженерная), tc_eff, G_prev (∂γ̇/∂σ).
void Grain::evaluate_Dp(const Vector6d& sigma, double dt)
{
    int n = 2 * params->num_slip_sys;

    for (int i = 0; i < n; i++)
        tau_buf(i) = sigma.dot(params->Ms_v[i]);

    gamma_dot.setZero();

    int sweeps = (params->hard_mode == 0) ? 1 : 3;
    for (int s = 0; s < sweeps; s++)
    {
        // эффективные τc от текущего итерата γ̇
        if (params->hard_mode == 2) {
            for (int i = 0; i < n; i++) {
                double dtc = 0.0;
                for (int j = 0; j < n; j++)
                    dtc += H(i, j) * gamma_dot(j);
                tc_eff(i) = tc[i] + dtc * dt;
            }
        } else if (params->hard_mode == 1) {
            double gamma_tot = 0.0;
            for (int i = 0; i < n; i++)
                gamma_tot += shifts[i] + gamma_dot(i) * dt;
            double tc_iso = params->tc0 + params->iso_a * std::pow(gamma_tot, params->iso_b);
            tc_eff.setConstant(tc_iso);
        } else {
            for (int i = 0; i < n; i++)
                tc_eff(i) = tc[i];
        }

        for (int i = 0; i < n; i++) {
            if (tau_buf(i) > tc_eff(i)) {
                double ratio = tau_buf(i) / tc_eff(i);
                if (ratio > RATIO_MAX) ratio = RATIO_MAX;   // защита от переполнения pow
                gamma_dot(i) = rate_factor_impl * std::pow(ratio, params->n_power);
            } else {
                gamma_dot(i) = 0.0;
            }
        }
    }

    // Dp в инженерной нотации Фойгта (сдвиги ×2 — как в find_SS)
    Dp_v.setZero();
    for (int i = 0; i < n; i++)
        if (gamma_dot(i) > 0.0)
            Dp_v.noalias() += gamma_dot(i) * params->MsE_v[i];

    // G = ∂γ̇/∂σ (n×6), та же схема sweep'ов; связка упрочнения — с предыдущего sweep'а
    G_prev.setZero();
    for (int s = 0; s < sweeps; s++)
    {
        for (int i = 0; i < n; i++) {
            if (gamma_dot(i) > 0.0) {
                // производные в точке клампированного ratio (осознанное приближение за клампом)
                double a_i = params->n_power * gamma_dot(i) / tau_buf(i);
                double b_i = -params->n_power * gamma_dot(i) / tc_eff(i);
                G_buf.row(i) = a_i * params->Ms_v[i].transpose();
                if (params->hard_mode == 2) {
                    for (int j = 0; j < n; j++)
                        if (H(i, j) != 0.0)
                            G_buf.row(i) += (b_i * H(i, j) * dt) * G_prev.row(j);
                } else if (params->hard_mode == 1) {
                    double gamma_tot = 0.0;
                    for (int j = 0; j < n; j++)
                        gamma_tot += shifts[j] + gamma_dot(j) * dt;
                    double dh = (gamma_tot > 1e-30)
                        ? params->iso_a * params->iso_b * std::pow(gamma_tot, params->iso_b - 1.0)
                        : 0.0;
                    G_buf.row(i) += (b_i * dh * dt) * G_prev.colwise().sum();
                }
            } else {
                G_buf.row(i).setZero();
            }
        }
        G_prev.swap(G_buf);
    }
}

// J = I + Δt · P · Σ_α MsE^α ⊗ G^α
void Grain::assemble_J(double dt)
{
    int n = 2 * params->num_slip_sys;
    dyad_buf.setZero();
    for (int i = 0; i < n; i++)
        if (gamma_dot(i) > 0.0)
            dyad_buf.noalias() += params->MsE_v[i] * G_prev.row(i);
    J_buf.noalias() = Matrix6d::Identity() + params->P_v * dyad_buf * dt;
}

// Ньютон по σ с line search Армихо; fallback — деление шага пополам.
// Требует: D_v (КСК, инженерная) и L (КСК) уже установлены, clear() выполнен.
void Grain::solve_implicit(double dt, int depth)
{
    int n = 2 * params->num_slip_sys;

    rate_factor_impl = params->strain_rate_infl * computeDeviatoricDu();

    trial_inc_v.noalias() = params->P_v * D_v;
    sigma_trial_v.noalias() = Stress_v + trial_inc_v * dt;

    Buf_v = sigma_trial_v;   // текущий итерат σ

    evaluate_Dp(Buf_v, dt);
    res_v = Buf_v - sigma_trial_v + params->P_v * Dp_v * dt;
    double scale = std::max(sigma_trial_v.norm(), params->tc0);
    double rnorm = res_v.norm();

    bool converged = (rnorm < NEWTON_TOL * scale);
    for (int iter = 0; iter < NEWTON_MAX_ITER && !converged; iter++)
    {
        assemble_J(dt);
        Eigen::PartialPivLU<Matrix6d> lu(J_buf);
        delta_sigma_v = lu.solve(-res_v);

        double lambda = 1.0;
        double rnorm_new = rnorm;
        bool accepted = false;
        for (int bt = 0; bt <= NEWTON_MAX_BACKTRACKS; bt++)
        {
            sigma_new_v = Buf_v + lambda * delta_sigma_v;
            evaluate_Dp(sigma_new_v, dt);
            res_new_v = sigma_new_v - sigma_trial_v + params->P_v * Dp_v * dt;
            rnorm_new = res_new_v.norm();
            if (rnorm_new <= (1.0 - 1e-4 * lambda) * rnorm) { accepted = true; break; }
            lambda *= 0.5;
        }
        if (!accepted) break;   // → fallback

        Buf_v = sigma_new_v;
        res_v = res_new_v;
        rnorm = rnorm_new;
        converged = (rnorm < NEWTON_TOL * scale);
    }

    if (!converged)
    {
        if (depth >= MAX_FALLBACK_DEPTH)
            throw std::runtime_error("implicit Newton failed to converge at max fallback depth");
        solve_implicit(0.5 * dt, depth + 1);
        solve_implicit(0.5 * dt, depth + 1);
        return;
    }

    // Обновление состояния по сошедшимся γ̇
    // (последний evaluate_Dp был именно в Buf_v → gamma_dot/tc_eff/G_prev соответствуют точке сходимости)
    Stress_v = Buf_v;
    for (int i = 0; i < n; i++) {
        sp_shifts[i] = gamma_dot(i);
        new_step[i] = (gamma_dot(i) > 0.0) ? 1 : 0;
        shifts[i] += gamma_dot(i) * dt;
        tc[i] = tc_eff(i);
    }
    if (params->hard_mode == 2)
        update_H_KM();

    apply_rotation(dt);
}
```

- [ ] **Step 4: Роутинг в `Simulate_Step`**

В `Grain::Simulate_Step` заменить блок подсчёта подшагов:

```cpp
    int n_substeps = 1;
    if (ds > DS_CRITICAL) {
        n_substeps = static_cast<int>(std::ceil(ds / DS_CRITICAL));
    }
```

на:

```cpp
    double ds_crit = (params->integ_mode == 1) ? DS_CRITICAL_IMPL : DS_CRITICAL;
    int n_substeps = 1;
    if (ds > ds_crit) {
        n_substeps = static_cast<int>(std::ceil(ds / ds_crit));
    }
```

и заменить блок упруго-пластического сплита:

```cpp
        // Elastic-plastic split for this substep
        dT = actual_dt;
        Sum_time = 0;
        elastic_buf = 0;

        while (actual_dt > Sum_time + TIME_EPSILON) {
            if (check_elastic(dT))
                Sum_time += dT;
            else {
                double dt_plastic = dT - elastic_buf;
                Sum_time += dT;
                if (dt_plastic > TIME_EPSILON)
                    find_SS(dt_plastic);
            }
            dT = actual_dt - Sum_time;
        }
```

на:

```cpp
        // Elastic-plastic split for this substep
        if (params->integ_mode == 1) {
            solve_implicit(actual_dt, 0);
        } else {
            dT = actual_dt;
            Sum_time = 0;
            elastic_buf = 0;

            while (actual_dt > Sum_time + TIME_EPSILON) {
                if (check_elastic(dT))
                    Sum_time += dT;
                else {
                    double dt_plastic = dT - elastic_buf;
                    Sum_time += dT;
                    if (dt_plastic > TIME_EPSILON)
                        find_SS(dt_plastic);
                }
                dT = actual_dt - Sum_time;
            }
        }
```

- [ ] **Step 5: Прогон новых тестов**

```
cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 /D_ENABLE_EXTENDED_ALIGNED_STORAGE tests.cpp Grain.cpp RVE.cpp MaterialParams.cpp
tests.exe --test-case="implicit single grain vs explicit"
tests.exe --test-case="implicit isotropic hardening"
```
Expected: оба PASS. Отладка при расхождении: проверить знак невязки (R = σ − σ_tr **+** Δt·P·Dᵖ), порог τ > τc_eff, ×2 в MsE_v.

- [ ] **Step 6: Полная регрессия**

```
tests.exe
```
Expected: все 10 тестов PASS (7 старых + parsing + 2 новых).

- [ ] **Step 7: Commit**

```bash
git add Grain.h Grain.cpp tests.cpp
git commit -m "feat: implicit Newton-Raphson stress integration (integ_mode=1)"
```

---

### Task 4: Эквивалентность на уровне RVE + стабильность большого шага

**Files:**
- Test: `tests.cpp`

**Interfaces:**
- Consumes: `rve_validation_impl.txt` (Task 1), роутинг `Simulate_Step` (Task 3).

**Примечания:** Ориентации зёрен детерминированы (default seed `rand()`), поэтому RVE из `rve_validation_impl.txt` имеет тот же набор ориентаций, что и RVE из `rve_validation.txt`, на котором получены референсные значения. `[slow]`-тесты — долгие (тот же порядок длительности, что существующие кривые, × стоимость Ньютона).

- [ ] **Step 1: Тесты**

Добавить в конец `tests.cpp`:

```cpp
TEST_CASE("implicit uniaxial tension equivalence [slow]") {
    RVE rve("rve_validation_impl.txt");
    Eigen::Matrix3d L = make_uniaxial_tension_L();

    struct { double strain; double stress_mpa; } points[] = {
        {0.040001,  47.5336},
        {0.0800003, 61.1466}
    };
    size_t cp_idx = 0;

    while (cp_idx < 2 && rve.get_accumulated_strain() < 0.08) {
        rve.Simulate_Step(dt, L);

        double acc = rve.get_accumulated_strain();
        if (acc >= points[cp_idx].strain) {
            double actual = rve.get_stress_intensity() / 1e6;
            CHECK(actual == doctest::Approx(points[cp_idx].stress_mpa).epsilon(0.05));
            cp_idx++;
        }
    }
}

TEST_CASE("implicit hard loading equivalence [slow]") {
    RVE rve("rve_validation_impl.txt");
    Eigen::Matrix3d L = make_uniaxial_tension_L();
    Eigen::Matrix3d L1 = make_plane_strain_L();

    struct { double strain; double stress_mpa; } points[] = {
        {0.0900014, 64.1854},
        {0.0910006, 55.8186},
        {0.3,        80.2589}
    };
    size_t cp_idx = 0;

    while (cp_idx < 3 && rve.get_accumulated_strain() < 0.3) {
        if (rve.get_accumulated_strain() < 0.09)
            rve.Simulate_Step(dt, L);
        else
            rve.Simulate_Step(dt, L1);

        double acc = rve.get_accumulated_strain();
        if (acc >= points[cp_idx].strain) {
            double actual = rve.get_stress_intensity() / 1e6;
            CHECK(actual == doctest::Approx(points[cp_idx].stress_mpa).epsilon(0.05));
            cp_idx++;
        }
    }
}

TEST_CASE("implicit big step x10") {
    RVE rve("rve_validation_impl.txt");
    Eigen::Matrix3d L = make_uniaxial_tension_L();
    double dt_big = 10e-3;

    struct { double strain; double stress_mpa; } points[] = {
        {0.040001,  47.5336},
        {0.0800003, 61.1466}
    };
    size_t cp_idx = 0;

    while (cp_idx < 2 && rve.get_accumulated_strain() < 0.08) {
        rve.Simulate_Step(dt_big, L);

        double acc = rve.get_accumulated_strain();
        if (acc >= points[cp_idx].strain) {
            double actual = rve.get_stress_intensity() / 1e6;
            CHECK(actual == actual);  // not NaN
            CHECK(actual == doctest::Approx(points[cp_idx].stress_mpa).epsilon(0.10));
            cp_idx++;
        }
    }
}

TEST_CASE("implicit big step x50") {
    RVE rve("rve_validation_impl.txt");
    Eigen::Matrix3d L = make_uniaxial_tension_L();
    double dt_big = 50e-3;

    while (rve.get_accumulated_strain() < 0.04) {
        rve.Simulate_Step(dt_big, L);
    }

    double actual = rve.get_stress_intensity() / 1e6;
    CHECK(actual == actual);  // not NaN
    CHECK(actual == doctest::Approx(47.5336).epsilon(0.10));
}
```

- [ ] **Step 2: Прогон быстрых тестов**

```
cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 /D_ENABLE_EXTENDED_ALIGNED_STORAGE tests.cpp Grain.cpp RVE.cpp MaterialParams.cpp
tests.exe --test-case="implicit big step*"
```
Expected: оба PASS. Если Ньютон не сходится на ×50 — не ослаблять допуски, а проверять кламп RATIO_MAX и line search.

- [ ] **Step 3: Прогон медленных тестов**

```
tests.exe --test-case="implicit * equivalence*"
```
Expected: оба PASS (допуск 5%, как у явных аналогов). Прогон долгий — порядка длительности исходных тестов кривых или дольше.

- [ ] **Step 4: Commit**

```bash
git add tests.cpp
git commit -m "test: RVE equivalence and big-step stability for implicit scheme"
```

---

### Task 5: Консистентный тангенс (Grain + RVE) и FD-верификация

**Files:**
- Modify: `Grain.h`, `Grain.cpp`, `RVE.h`, `RVE.cpp`
- Test: `tests.cpp`

**Interfaces:**
- Consumes: `solve_implicit`, `assemble_J`, `G_prev`, `gamma_dot` (Task 3).
- Produces: `void Grain::get_tangent_crystal_into(Matrix6d& result) const`; `void Grain::get_tangent_lab_into(Matrix6d& result) const`; `void RVE::get_tangent_lab_into(Matrix6d& result) const`. Тангенс = ∂σ_{n+1}/∂Δε = J⁻¹·P на сошедшемся состоянии (КСК, Фойгт).

- [ ] **Step 1: Failing test — FD-верификация тангенса**

Добавить в конец `tests.cpp`:

```cpp
TEST_CASE("implicit tangent FD") {
    Grain gr("grain_implicit.txt");

    double eps_dot_local = 0.017;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 0) = -0.5 * eps_dot_local;
    L(1, 1) = -0.5 * eps_dot_local;
    L(2, 2) = eps_dot_local;

    double dt_local = 1e-3;
    for (int i = 0; i < 49; i++)
        gr.Simulate_Step(dt_local, L);   // установившееся пластическое состояние

    Grain gr_base = gr;                  // состояние S49 (копия)
    gr.Simulate_Step(dt_local, L);       // шаг S49→S50, тангенс соответствует этому шагу

    Matrix6d T = Matrix6d::Zero();
    gr.get_tangent_crystal_into(T);

    double h = 1e-6;                     // возмущение скорости деформации
    for (int k = 0; k < 6; k++) {
        auto ij = voigtToTensorIndex(k);
        Eigen::Matrix3d dL = Eigen::Matrix3d::Zero();
        if (ij.first == ij.second)
            dL(ij.first, ij.second) = h;
        else {
            dL(ij.first, ij.second) = 0.5 * h;   // симметрично → D_v(k) = h
            dL(ij.second, ij.first) = 0.5 * h;
        }

        Grain gp = gr_base;
        gp.Simulate_Step(dt_local, L + dL);
        Grain gm = gr_base;
        gm.Simulate_Step(dt_local, L - dL);

        Vector6d col_fd = (gp.get_stress() - gm.get_stress()) / (2.0 * h * dt_local);
        for (int r = 0; r < 6; r++) {
            double diff = std::abs(col_fd(r) - T(r, k));
            CHECK(diff < 1e-3 * std::abs(T(r, k)) + 1e7);   // отн. 1e-3, абс. 1e7 Па
        }
    }
}
```

Run: `tests.exe --test-case="implicit tangent FD"`
Expected: FAIL (нет метода `get_tangent_crystal_into` — ошибка компиляции).

- [ ] **Step 2: Реализация — тангенс в Grain**

В `Grain.h` private-секцию (рядом с буферами неявной схемы из Task 3) добавить:

```cpp
    mutable Matrix6d tangent_buf;   // консистентный тангенс последнего неявного шага (КСК)
    bool tangent_valid = false;
```

В public-секцию `Grain.h` (рядом с `get_elasticity_lab_into`) добавить:

```cpp
    void get_tangent_crystal_into(Matrix6d& result) const;  // J^-1·P сошедшегося шага (КСК)
    void get_tangent_lab_into(Matrix6d& result) const;      // то же в ЛСК
```

В `Grain.cpp` в `init_common()` добавить (рядом с обнулением буферов неявной схемы):

```cpp
    tangent_buf.setZero();
    tangent_valid = false;
```

В `Grain.cpp` в `solve_implicit()` после блока `if (params->hard_mode == 2) update_H_KM();` и ПЕРЕД `apply_rotation(dt);` вставить:

```cpp
    // Консистентный тангенс: J^{-1}·P на сошедшемся состоянии
    // (gamma_dot/G_prev соответствуют точке сходимости — последний evaluate_Dp был в Buf_v)
    assemble_J(dt);
    {
        Eigen::PartialPivLU<Matrix6d> lu(J_buf);
        tangent_buf = lu.solve(params->P_v);
    }
    tangent_valid = true;
```

В `Grain.cpp` добавить методы (после `get_elasticity_lab_into`):

```cpp
void Grain::get_tangent_crystal_into(Matrix6d& result) const
{
    if (!tangent_valid)
        throw std::runtime_error("get_tangent_crystal_into: no implicit step performed yet");
    result = tangent_buf;
}

void Grain::get_tangent_lab_into(Matrix6d& result) const
{
    if (!tangent_valid)
        throw std::runtime_error("get_tangent_lab_into: no implicit step performed yet");
    if (!T_voigt_valid) {
        computeVoigtRotationInPlace(T_voigt_buf);
        T_voigt_valid = true;
    }
    Buf_P.noalias() = T_voigt_buf * tangent_buf;
    result.noalias() = Buf_P * T_voigt_buf.transpose();
}
```

- [ ] **Step 3: Прогон FD-теста**

```
cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 /D_ENABLE_EXTENDED_ALIGNED_STORAGE tests.cpp Grain.cpp RVE.cpp MaterialParams.cpp
tests.exe --test-case="implicit tangent FD"
```
Expected: PASS. При расхождении > порога: типичные причины — пропущены члены упрочнения в G (b_i-часть), не тот dt в assemble_J, G_prev не из точки сходимости.

- [ ] **Step 4: Тангенс RVE + тест**

В `RVE.h` public-секцию (после `get_stress_intensity`) добавить:

```cpp
    void get_tangent_lab_into(Matrix6d& result) const;  // осреднённый тангенс RVE (Taylor, ЛСК)
```

В `RVE.cpp` добавить (после `get_stress_intensity`):

```cpp
void RVE::get_tangent_lab_into(Matrix6d& result) const
{
    result.setZero();
    int n = static_cast<int>(grains.size());
    for (int i = 0; i < n; i++) {
        double w = grains[i].get_volume() / total_volume;
        grains[i].get_tangent_lab_into(temp_P_buf);
        result.noalias() += w * temp_P_buf;
    }
}
```

Добавить в конец `tests.cpp`:

```cpp
TEST_CASE("implicit RVE tangent") {
    RVE rve("rve_validation_impl.txt");
    Eigen::Matrix3d L = make_uniaxial_tension_L();

    for (int i = 0; i < 50; i++)
        rve.Simulate_Step(dt, L);

    Matrix6d T;
    rve.get_tangent_lab_into(T);

    for (int r = 0; r < 6; r++)
        for (int c = 0; c < 6; c++)
            CHECK(T(r, c) == T(r, c));   // not NaN

    CHECK(T(0, 0) > 1e10);                              // упругий порядок величины
    CHECK(std::abs(T(0, 1) - T(1, 0)) < 1e-2 * std::abs(T(0, 1)) + 1e6);  // грубая симметрия
}
```

Run: `tests.exe --test-case="implicit RVE tangent"`
Expected: PASS.

- [ ] **Step 5: Полная регрессия**

```
tests.exe -tce="*[slow]*"
```
Expected: все не-slow тесты PASS.

- [ ] **Step 6: Commit**

```bash
git add Grain.h Grain.cpp RVE.h RVE.cpp tests.cpp
git commit -m "feat: consistent tangent for implicit scheme (Grain + RVE)"
```

---

### Task 6: Документация (CONTEXT.md + README.md)

**Files:**
- Modify: `CONTEXT.md`, `README.md`

**Interfaces:**
- Consumes: вся реализация Tasks 1–5.

- [ ] **Step 1: CONTEXT.md**

Обновить:
- «Текущий статус»: добавить «+ неявная схема интегрирования (integ_mode=1) с консистентным тангенсом».
- Файловая структура: добавить `grain_implicit.txt`, `grain_test_iso_impl.txt`, `rve_validation_impl.txt` с пояснениями.
- Структура MaterialParams: добавить поля `MsE_v` (инженерная нотация, сдвиги ×2), `integ_mode` (0=явная, 1=неявная).
- Формат входного файла материала: добавить опциональную последнюю строку `integ_mode`.
- Новый раздел «Неявная схема (integ_mode=1)»:

```markdown
## Неявная схема (integ_mode=1)

Newton–Raphson по напряжениям (6×6): R(σ) = σ − σ_tr + Δt·P:Dᵖ(σ) = 0,
J = I + Δt·P·Σ MsE^α⊗G^α. Пороговый закон Хатчинсона как в явной схеме
(γ̇=0 при τ≤τc_eff). Упрочнение fully implicit: fixed-point 3 sweep'а
γ̇↔τc_eff, связка ∂Δγ/∂σ с предыдущего sweep'а (приближение Калидинди).

- Line search Армихо (до 5 откатов), кламп τ/τc ≤ RATIO_MAX=2.0
- Fallback: деление шага пополам, глубина ≤ 8
- NEWTON_TOL=1e-8 (относительный), NEWTON_MAX_ITER=20
- Подшагивание: DS_CRITICAL_IMPL=1e-4 (отдельно от явной DS_CRITICAL=1.7e-6)
- Консистентный тангенс: ∂Δσ/∂Δε = J⁻¹·P; RVE — осреднение в ЛСК (Taylor)
- Методы: Grain::get_tangent_crystal_into, Grain::get_tangent_lab_into, RVE::get_tangent_lab_into
- Верификация: FD-тест тангенса (отн. 1e-3), эквивалентность кривых 5% ([slow]), big-step ×10/×50 10%
```

- Таблица автотестов: добавить 6 новых тестов.
- Методы Grain/RVE: добавить `get_tangent_crystal_into`, `get_tangent_lab_into`.

- [ ] **Step 2: README.md**

Обновить:
- «Возможности»: добавить неявную схему и консистентный тангенс.
- Математическая постановка: подраздел про неявную интеграцию (формулы из CONTEXT.md выше, сжато).
- Формат входного файла: строка `integ_mode` (опционально).
- Структура проекта: новые конфиг-файлы.
- Число тестов: «13 тестов» (7 старых + parsing + single grain + iso + 2 big-step + tangent FD + RVE tangent; 2 equivalence — [slow], итого 15 кейсов).

- [ ] **Step 3: Финальная регрессия**

```
cl /EHsc /O2 /openmp /I D:\eigen-3.4.0 /D_ENABLE_EXTENDED_ALIGNED_STORAGE tests.cpp Grain.cpp RVE.cpp MaterialParams.cpp
tests.exe
```
Expected: все тесты PASS (включая [slow]).

- [ ] **Step 4: Commit**

```bash
git add CONTEXT.md README.md
git commit -m "docs: implicit integration scheme in CONTEXT and README"
```

---

## Definition of Done

- [ ] Все тесты зелёные: `tests.exe` (включая `[slow]`)
- [ ] Тангенс верифицирован конечными разностями (отн. 1e-3)
- [ ] Неявная схема воспроизводит референсные кривые (5%, same dt) и стабильна при dt ×10/×50 (10%)
- [ ] Явная схема не изменена (7 исходных тестов без правок проходят)
- [ ] CONTEXT.md и README.md отражают новую схему
