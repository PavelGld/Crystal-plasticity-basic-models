// Замена тестов тангенса: сравнение одноосного пластического модуля E_t = 1/(C⁻¹)₂₂
// аналитического тангенса и полного 6×6 тангенса по центральным разностям.
// Компонента T(2,2) ≈ K + (малая пластическая добавка) слепа к упрочнению: тангенс
// идеальной пластичности отличается от правильного по T(2,2) на ~0.6%.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "cp_model/Grain.hpp"
#include <cmath>

// Полный FD-тангенс в КСК: столбец k = dσ_v/dD_v(k), D_v(k) возмущается на ±h
// (сдвиги инженерные, как в D_v). Ориентация зерна должна быть единичной (ЛСК = КСК).
static Matrix6d fd_tangent_crystal(const Grain& base, const Eigen::Matrix3d& L, double dtl, double h)
{
    Matrix6d C;
    for (int k = 0; k < 6; k++) {
        auto ij = voigtToTensorIndex(k);
        Eigen::Matrix3d dL = Eigen::Matrix3d::Zero();
        if (ij.first == ij.second) dL(ij.first, ij.second) = h;
        else { dL(ij.first, ij.second) = 0.5 * h; dL(ij.second, ij.first) = 0.5 * h; }
        Grain gp = base; gp.Simulate_Step(dtl, L + dL);
        Grain gm = base; gm.Simulate_Step(dtl, L - dL);
        C.col(k) = (gp.get_stress() - gm.get_stress()) / (2.0 * h * dtl);
    }
    return C;
}
// Одноосный касательный модуль: dσ33/dε33 при σ_остальные = 0
static double uniaxial_Et(const Matrix6d& C) { return 1.0 / C.inverse()(2, 2); }

static Eigen::Matrix3d uniax_L(double e) {
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 0) = -0.5 * e; L(1, 1) = -0.5 * e; L(2, 2) = e;
    return L;
}

TEST_CASE("implicit tangent FD: uniaxial plastic modulus (KM)") {
    Grain gr("grain_implicit.txt");
    const double ed = 0.017, dtl = 1e-3;
    Eigen::Matrix3d L = uniax_L(ed);
    for (int i = 0; i < 49; i++) gr.Simulate_Step(dtl, L);
    Grain base = gr;
    gr.Simulate_Step(dtl, L);
    Matrix6d T; gr.get_tangent_crystal_into(T);
    Matrix6d C = fd_tangent_crystal(base, L, dtl, 1e-3 * ed);

    double Et_a = uniaxial_Et(T), Et_fd = uniaxial_Et(C);
    MESSAGE("KM: E_t analytic = " << Et_a << ", FD = " << Et_fd << ", T22 = " << T(2, 2));
    CHECK(Et_fd > 0.0);
    CHECK(Et_fd < 0.05 * T(2, 2));                           // действительно пластическое состояние
    CHECK(Et_a == doctest::Approx(Et_fd).epsilon(1e-3));
}

TEST_CASE("implicit tangent FD: uniaxial plastic modulus, isotropic, multi-substep") {
    Grain gr("grain_test_iso_impl.txt");
    const double ed = 0.017;
    Eigen::Matrix3d L = uniax_L(ed);
    for (int i = 0; i < 200; i++) gr.Simulate_Step(1e-3, L);
    const double dtl = 5e-3 / ed;                            // Δε33 = 5e-3 → 3 подшага при DS_CRIT = 2e-3
    Grain base = gr;
    gr.Simulate_Step(dtl, L);
    Matrix6d T; gr.get_tangent_crystal_into(T);
    Matrix6d C = fd_tangent_crystal(base, L, dtl, 1e-4 * ed);

    double Et_a = uniaxial_Et(T), Et_fd = uniaxial_Et(C);
    MESSAGE("iso: E_t analytic = " << Et_a << ", FD = " << Et_fd << ", T22 = " << T(2, 2));
    CHECK(Et_fd > 0.0);
    CHECK(Et_fd < 0.05 * T(2, 2));
    CHECK(Et_a == doctest::Approx(Et_fd).epsilon(1e-3));
}
