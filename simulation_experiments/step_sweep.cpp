#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "Grain.h"
#include "RVE.h"

using namespace std::chrono;

// Traced from solve_implicit
static int   newton_iters_trace    = 0;
static int   newton_fallbacks_trace = 0;
static int   newton_depth_sum       = 0;
static double newton_rnorm_trace    = 0;

// Hooks: these are defined in Grain.cpp — compile-time switch
// We use a simple approach: run with dt and capture via internal prints.
// Alternative: extract stats from Grain after each step.
// For this experiment we measure externally: stress after fixed strain.

static Eigen::Matrix3d make_L(double e) {
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0,0) = -0.5*e; L(1,1) = -0.5*e; L(2,2) = e;
    return L;
}

struct Result {
    double dt;
    double stress_mpa;
    double ms_per_step;
    int    steps_to_target;
};

void sweep_grain(const char* file, double eps_dot, double target_strain,
                 const double* dts, int ndt, std::vector<Result>& out)
{
    Eigen::Matrix3d L = make_L(eps_dot);
    for (int i = 0; i < ndt; i++) {
        double dt = dts[i];
        if (dt * eps_dot * 0.5 > 1.0) continue;  // sanity: not too huge
        Grain gr(file);
        int steps = 0;
        auto t0 = high_resolution_clock::now();
        while (gr.get_acc_def() < target_strain) {
            gr.Simulate_Step(dt, L);
            steps++;
        }
        auto t1 = high_resolution_clock::now();
        double stress = gr.get_stress_intensity() / 1e6;
        double ms = duration_cast<milliseconds>(t1 - t0).count();
        Result r = {dt, stress, ms / steps, steps};
        out.push_back(r);
        std::cout << "  grain " << file << " dt=" << std::setw(8) << dt
                  << " stress=" << std::fixed << std::setprecision(4) << stress
                  << " MPa  steps=" << steps << "  ms/step=" << ms/steps << "\n";
    }
}

void sweep_rve(const char* config, double eps_dot, double target_strain,
               const double* dts, int ndt, std::vector<Result>& out)
{
    Eigen::Matrix3d L = make_L(eps_dot);
    for (int i = 0; i < ndt; i++) {
        double dt = dts[i];
        if (dt * eps_dot * 0.5 > 1.0) continue;
        RVE rve(config);
        int steps = 0;
        auto t0 = high_resolution_clock::now();
        while (rve.get_accumulated_strain() < target_strain) {
            rve.Simulate_Step(dt, L);
            steps++;
        }
        auto t1 = high_resolution_clock::now();
        double stress = rve.get_stress_intensity() / 1e6;
        double ms = duration_cast<milliseconds>(t1 - t0).count();
        Result r = {dt, stress, ms / steps, steps};
        out.push_back(r);
        std::cout << "  RVE   " << config << " dt=" << std::setw(8) << dt
                  << " stress=" << std::fixed << std::setprecision(4) << stress
                  << " MPa  steps=" << steps << "  ms/step=" << ms/steps << "\n";
    }
}

int main() {
    double dts[] = { 1e-3, 2.5e-3, 5e-3, 1e-2, 2.5e-2, 5e-2, 0.1, 0.25, 0.5 };
    int ndt = sizeof(dts)/sizeof(dts[0]);

    std::ofstream f("simulation_experiments/step_sweep_results.txt");

    // --- Single grain, reference rate ---
    std::cout << "=== Single grain, ref rate (eps=0.0017), target eps=0.04 ===\n";
    std::vector<Result> res_g;
    sweep_grain("grain_implicit.txt", 0.0017, 0.04, dts, ndt, res_g);
    f << "# Single grain implicit, ref rate, target eps=0.04\n";
    f << "# dt  stress(MPa)  ms_per_step  steps\n";
    for (auto& r : res_g) f << r.dt << " " << r.stress_mpa << " " << r.ms_per_step << " " << r.steps_to_target << "\n";

    // --- Single grain, 10x rate ---
    std::cout << "\n=== Single grain, 10x rate (eps=0.017), target eps=0.04 ===\n";
    std::vector<Result> res_gf;
    sweep_grain("grain_implicit.txt", 0.017, 0.04, dts, ndt, res_gf);
    f << "\n# Single grain implicit, 10x rate, target eps=0.04\n";
    for (auto& r : res_gf) f << r.dt << " " << r.stress_mpa << " " << r.ms_per_step << " " << r.steps_to_target << "\n";

    // --- RVE, reference rate ---
    std::cout << "\n=== RVE (1000 grains), ref rate, target eps=0.04 ===\n";
    std::vector<Result> res_rve;
    sweep_rve("rve_validation_impl.txt", 0.0017, 0.04, dts, ndt, res_rve);
    f << "\n# RVE implicit 1000 grains, ref rate, target eps=0.04\n";
    f << "# dt  stress(MPa)  ms_per_step  steps\n";
    for (auto& r : res_rve) f << r.dt << " " << r.stress_mpa << " " << r.ms_per_step << " " << r.steps_to_target << "\n";

    // --- RVE, 10x rate ---
    std::cout << "\n=== RVE (1000 grains), 10x rate, target eps=0.02 ===\n";
    std::vector<Result> res_rvef;
    sweep_rve("rve_validation_impl.txt", 0.017, 0.02, dts, ndt, res_rvef);
    f << "\n# RVE implicit 1000 grains, 10x rate, target eps=0.02\n";
    for (auto& r : res_rvef) f << r.dt << " " << r.stress_mpa << " " << r.ms_per_step << " " << r.steps_to_target << "\n";

    // --- Explicit reference at dt=1e-3 for comparison ---
    std::cout << "\n=== Reference: explicit RVE dt=1e-3, ref rate, target eps=0.04 ===\n";
    {
        RVE rve("rve_validation.txt");
        auto L = make_L(0.0017);
        int steps = 0;
        auto t0 = high_resolution_clock::now();
        while (rve.get_accumulated_strain() < 0.04) { rve.Simulate_Step(1e-3, L); steps++; }
        auto t1 = high_resolution_clock::now();
        double stress = rve.get_stress_intensity() / 1e6;
        double ms = duration_cast<milliseconds>(t1 - t0).count();
        f << "\n# Reference explicit RVE dt=1e-3: stress=" << stress << " MPa, ms/step=" << ms/steps << ", steps=" << steps << "\n";
        std::cout << "  stress=" << stress << " MPa, ms/step=" << ms/steps << ", steps=" << steps << "\n";
    }

    // --- Stress-strain curve at dt=0.1 (for visualization) ---
    std::cout << "\n=== Stress-strain curve: implicit RVE dt=0.1, ref rate ===\n";
    {
        RVE rve("rve_validation_impl.txt");
        auto L = make_L(0.0017);
        std::ofstream fc("simulation_experiments/curve_dt0.1.txt");
        double dt = 0.1;
        int step = 0;
        while (rve.get_accumulated_strain() < 0.1) {
            rve.Simulate_Step(dt, L);
            step++;
            fc << rve.get_accumulated_strain() << " " << rve.get_stress_intensity()/1e6 << "\n";
        }
        std::cout << "  Written " << step << " points to curve_dt0.1.txt\n";
    }

    f.close();
    std::cout << "\nResults written to simulation_experiments/step_sweep_results.txt\n";
    return 0;
}
