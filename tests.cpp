#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "Grain.h"
#include "RVE.h"
#include <cmath>

static constexpr double eps_dot = 0.0017;
static constexpr double dt = 1e-3;

static Eigen::Matrix3d make_uniaxial_tension_L() {
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 0) = -0.5 * eps_dot;
    L(1, 1) = -0.5 * eps_dot;
    L(2, 2) = eps_dot;
    return L;
}

static Eigen::Matrix3d make_plane_strain_L() {
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 0) = eps_dot;
    L(1, 1) = -eps_dot;
    return L;
}

TEST_CASE("no rotation") {
    Grain gr("grain.txt");
    Eigen::AngleAxisd aa_init(3.14159265358979323846 / 6.0, Eigen::Vector3d::UnitZ());
    gr.set_orientation(Eigen::Quaterniond(aa_init));
    Eigen::Matrix3d R_before = gr.get_orientation().toRotationMatrix();

    Eigen::Matrix3d L = make_uniaxial_tension_L();
    for (int i = 0; i < 100; i++)
        gr.Simulate_Step(dt, L);

    double diff_norm = (gr.get_orientation().toRotationMatrix() - R_before).norm();
    CHECK(diff_norm < 1e-10);
}

TEST_CASE("simple shear rotation") {
    // Uses rot_mode=1 (material rotation, no plastic spin)
    Grain gr("grain_test_shear.txt");

    double gamma = 0.001;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 1) = gamma;

    int N = 100;
    double expected_angle = gamma * N * dt / 2.0;

    for (int i = 0; i < N; i++)
        gr.Simulate_Step(dt, L);

    double actual_angle = 2.0 * std::acos(std::abs(gr.get_orientation().w()));
    CHECK(std::abs(actual_angle - expected_angle) < 1e-5);
}

TEST_CASE("RVE basic") {
    RVE rve("rve_validation.txt");
    CHECK(rve.num_grains() == 1000);

    Eigen::Matrix3d L = make_uniaxial_tension_L();
    for (int i = 0; i < 1000; i++)
        rve.Simulate_Step(dt, L);

    Eigen::Matrix3d stress = rve.get_stress();
    CHECK(stress(0, 0) == stress(0, 0)); // not NaN
}

TEST_CASE("hard loading") {
    RVE rve("rve_validation.txt");
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

TEST_CASE("uniaxial tension") {
    RVE rve("rve_validation.txt");
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

TEST_CASE("isotropic hardening") {
    Grain gr("grain_test_iso.txt");

    double eps_dot_local = 0.017;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 0) = -0.5 * eps_dot_local;
    L(1, 1) = -0.5 * eps_dot_local;
    L(2, 2) = eps_dot_local;

    double tc0 = gr.get_tc(0);
    double dt_local = 1e-3;
    int steps = 1000;
    for (int i = 0; i < steps; i++)
        gr.Simulate_Step(dt_local, L);

    double tc_after = gr.get_tc(0);
    CHECK(tc_after > tc0);

    size_t n = gr.num_slip_systems();
    for (size_t i = 1; i < n; i++)
        CHECK(gr.get_tc(i) == doctest::Approx(tc_after));
}
