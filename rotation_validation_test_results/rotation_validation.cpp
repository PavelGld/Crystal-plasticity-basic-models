#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "RVE.h"

#define M_PI 3.14159265358979323846

// ============================================================
// Unit tests (single grain / RVE sanity checks)
// ============================================================

static void test_uniaxial_stretch_no_rotation()
{
    std::cout << "=== Test 1: Uniaxial stretch (no rotation expected) ===\n";

    Grain gr("grain.txt");

    Eigen::AngleAxisd aa_init(M_PI / 6.0, Eigen::Vector3d::UnitZ());
    Eigen::Quaterniond q_init(aa_init);
    gr.set_orientation(q_init);

    Eigen::Quaterniond q_before = gr.get_orientation();
    Eigen::Matrix3d R_before = q_before.toRotationMatrix();

    double eps = 0.001;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 0) = eps;
    L(1, 1) = -0.5 * eps;
    L(2, 2) = -0.5 * eps;

    double dt = 1e-3;
    int N = 100;
    for (int i = 0; i < N; i++)
        gr.Simulate_Step(dt, L);

    Eigen::Quaterniond q_after = gr.get_orientation();
    Eigen::Matrix3d R_after = q_after.toRotationMatrix();

    double diff_norm = (R_after - R_before).norm();

    std::cout << "Difference norm: " << diff_norm << "\n";
    if (diff_norm < 1e-10)
        std::cout << "TEST 1 PASSED: Orientation unchanged\n\n";
    else
        std::cout << "TEST 1 FAILED: Orientation changed unexpectedly!\n\n";
}

static void test_simple_shear_rotation()
{
    std::cout << "=== Test 2: Simple shear (rotation expected) ===\n";

    Grain gr("grain.txt");

    Eigen::Quaterniond q_before = gr.get_orientation();
    Eigen::Matrix3d R_before = q_before.toRotationMatrix();

    double eps = 0.001;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 1) = eps;

    double dt = 1e-3;
    int N = 100;
    double total_time = N * dt;
    double expected_angle = eps * total_time / 2.0;

    for (int i = 0; i < N; i++)
        gr.Simulate_Step(dt, L);

    Eigen::Quaterniond q_after = gr.get_orientation();
    double actual_angle = 2.0 * std::acos(std::abs(q_after.w()));

    std::cout << "Expected rotation angle: " << expected_angle << " rad ("
              << expected_angle * 180.0 / M_PI << " deg)\n";
    std::cout << "Actual rotation angle: " << actual_angle << " rad ("
              << actual_angle * 180.0 / M_PI << " deg)\n";

    double angle_error = std::abs(actual_angle - expected_angle);
    if (angle_error < 1e-5)
        std::cout << "TEST 2 PASSED: Rotation angle matches expected\n\n";
    else
        std::cout << "TEST 2 FAILED: Rotation angle mismatch!\n\n";
}

static void test_rve_basic()
{
    std::cout << "=== Test 3: RVE basic functionality ===\n";

    RVE rve("rve_test.txt");
    std::cout << "Number of grains: " << rve.num_grains() << "\n";

    double eps = 0.001;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 0) = eps;
    L(1, 1) = -0.5 * eps;
    L(2, 2) = -0.5 * eps;

    double dt = 1e-3;
    for (int i = 0; i < 1000; i++)
        rve.Simulate_Step(dt, L);

    Eigen::Matrix3d stress = rve.get_stress();
    std::cout << "Average stress:\n" << stress / 1e6 << " MPa\n";
    std::cout << "Time: " << rve.get_time() << "\n";

    if (rve.num_grains() == 1000)
        std::cout << "TEST 3 PASSED: RVE created 1000 grains\n\n";
    else
        std::cout << "TEST 3 FAILED: Wrong number of grains\n\n";
}

static void test_soft_loading()
{
    std::cout << "=== Test 4: Soft loading (uniaxial tension) ===\n";

    RVE rve("rve_test.txt");

    double stress_rate = 1e8;
    Eigen::Matrix3d dS_dt = Eigen::Matrix3d::Zero();
    dS_dt(0, 0) = stress_rate;

    double dt = 1e-5;
    int n_steps = 100;
    for (int i = 0; i < n_steps; i++)
        rve.Simulate_Step_Soft(dt, dS_dt);

    Eigen::Matrix3d stress = rve.get_stress();
    std::cout << "Final stress:\n" << stress / 1e6 << " MPa\n";

    double s11 = std::abs(stress(0, 0));
    double off_diag_max = std::max({
        std::abs(stress(0, 1)), std::abs(stress(0, 2)),
        std::abs(stress(1, 2)), std::abs(stress(1, 1)),
        std::abs(stress(2, 2))
    });

    if (off_diag_max < 0.1 * s11)
        std::cout << "TEST 4 PASSED: Stress is approximately uniaxial\n\n";
    else
        std::cout << "TEST 4 FAILED: Off-diagonal components too large\n\n";
}

void run_all_tests()
{
    //test_uniaxial_stretch_no_rotation();
    //test_simple_shear_rotation();
    //test_rve_basic();
    test_soft_loading();
    std::cout << "=== All tests completed ===\n";
}

// ============================================================
// Orientation export helpers
// ============================================================

static void export_orientations(const RVE& rve, const std::string& filename, double strain, const std::string& loading_type) {
    std::ofstream out(filename);
    out << "# Loading: " << loading_type << ", Strain: " << strain << "\n";
    out << "# qw qx qy qz\n";
    for (size_t i = 0; i < rve.num_grains(); i++) {
        Eigen::Quaterniond q = rve.get_grain(i).get_orientation();
        out << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << "\n";
    }
    out.close();
}

// ============================================================
// Texture simulation (4 loading modes, checkpoints at 30/50/100%)
// ============================================================

static void run_loading(const std::string& config_file, const Eigen::Matrix3d& L, double target_strain,
                        double dt, const std::string& loading_type,
                        const std::vector<double>& checkpoints) {

    RVE rve(config_file);
    size_t current_checkpoint = 0;

    std::cout << "\n--- " << loading_type << " ---\n";
    std::cout << "RVE grains: " << rve.num_grains() << "\n";
    std::cout << "Target strain (equiv): " << target_strain * 100 << "%\n";
    std::cout << "Step size: dt = " << dt << "\n";

    int step_count = 0;
    while (rve.get_accumulated_strain() < target_strain) {
        rve.Simulate_Step(dt, L);
        step_count++;

        double acc = rve.get_accumulated_strain();

        if (current_checkpoint < checkpoints.size()) {
            double cp = checkpoints[current_checkpoint];
            if (acc >= cp) {
                std::cout << "  Reached " << cp * 100 << "% equiv strain at step " << step_count
                          << " (actual: " << acc * 100 << "%)\n";

                std::string result_file = "orientations_" + loading_type + "_" + std::to_string(static_cast<int>(cp * 100)) + "pct.txt";
                export_orientations(rve, result_file, acc, loading_type);
                current_checkpoint++;
            }
        }

        if (step_count % 500 == 0)
            std::cout << "  Step " << step_count << ": " << acc * 100 << "%\n";
    }

    std::cout << "Completed: " << loading_type
              << " final_equiv_strain = " << rve.get_accumulated_strain() * 100 << "%\n";

    Eigen::Matrix3d stress = rve.get_stress();
    std::cout << "Final stress:\n" << stress / 1e6 << " MPa\n";

    Eigen::Matrix3d Dp = rve.get_plastic_strain_rate();
    std::cout << "Plastic strain rate:\n" << Dp << "\n";
}

void run_texture_simulation()
{
    std::cout << "=== Texture Validation (4 loading modes) ===\n";

    const std::string config_file = "rve_validation.txt";
    const double target_strain = 1.0;
    const double dt = 1e-4;
    std::vector<double> checkpoints = {0.3, 0.5, 1.0};

    // 1. Simple Shear
    {
        double gamma_dot = 0.1;
        Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
        L(0, 1) = gamma_dot;
        run_loading(config_file, L, target_strain, dt, "simple_shear", checkpoints);
    }

    // 2. Uniaxial Tension
    {
        double eps_dot = 0.1;
        Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
        L(0, 0) = -0.5 * eps_dot;
        L(1, 1) = -0.5 * eps_dot;
        L(2, 2) = eps_dot;
        run_loading(config_file, L, target_strain, dt, "uniaxial_tension", checkpoints);
    }

    // 3. Uniaxial Compression
    {
        double eps_dot = 0.1;
        Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
        L(0, 0) = 0.5 * eps_dot;
        L(1, 1) = 0.5 * eps_dot;
        L(2, 2) = -eps_dot;
        run_loading(config_file, L, target_strain, dt, "uniaxial_compression", checkpoints);
    }

    // 4. Plane Strain: L = diag(0, -eps_dot, eps_dot)
    {
        double eps_dot = 0.1;
        Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
        L(1, 1) = -eps_dot;
        L(2, 2) = eps_dot;
        run_loading(config_file, L, target_strain, dt, "plane_strain", checkpoints);
    }

    std::cout << "\n=== All simulations completed ===\n";
    std::cout << "Output files: orientations_<loading>_<strain>.txt\n";
    std::cout << "Run plot_polefigures.m in MATLAB to build pole figures\n";
}

// ============================================================
// Main
// ============================================================

int main(int argc, char* argv[])
{
    
	run_all_tests();
    //run_texture_simulation();
    
    return 0;
}
