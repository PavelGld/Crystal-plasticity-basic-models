#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <string>
#include <cmath>
#include <fstream>
#include "Grain.h"
#include "RVE.h"

#define M_PI 3.14159265358979323846

int main()
{
    std::ofstream out;          // ����� ��� ������
    out.open("Stress_strain_hard_loading.txt");      // ��������� ���� ��� ������

    const std::string config_file = "rve_validation.txt";
    const std::string loading_type = "========= HARD LOADING =========";
    const double target_strain = 0.3;
    const double dt = 1e-3;
    std::vector<double> checkpoints = {};
    checkpoints.resize(1000);
    for (unsigned i = 0; i < 1000; i++)
        checkpoints[i] = 0.001*i;
    double eps_dot = 0.0017;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0, 0) = -0.5 * eps_dot;
    L(1, 1) = -0.5 * eps_dot;
    L(2, 2) = eps_dot;

    Eigen::Matrix3d L1 = Eigen::Matrix3d::Zero();
    L1(0, 0) = eps_dot;
    L1(1, 1) = -eps_dot;


    RVE rve(config_file);
    size_t current_checkpoint = 0;

    std::cout << "\n--- " << loading_type << " ---\n";
    std::cout << "RVE grains: " << rve.num_grains() << "\n";
    std::cout << "Target strain (equiv): " << target_strain * 100 << "%\n";
    std::cout << "Step size: dt = " << dt << "\n";

    int step_count = 0;
    while (rve.get_accumulated_strain() < target_strain) {
        if (rve.get_accumulated_strain() < 0.09)
            rve.Simulate_Step(dt, L);
        else
            rve.Simulate_Step(dt, L1);
        step_count++;

        double acc = rve.get_accumulated_strain();

        if (current_checkpoint < checkpoints.size()) {
            double cp = checkpoints[current_checkpoint];
            if (acc >= cp) {
                Eigen::Matrix3d str = rve.get_stress();
                str -= str.trace() * Eigen::Matrix3d::Identity();

                if (str(0, 0) != str(0, 0))
                    std::cout << "NAN\n";

                out << acc << "    " << rve.get_stress_intensity() / 1e6 << '\n';

                std::cout << "  Reached " << cp * 100 << "% equiv strain at step " << step_count
                    << " (actual: " << acc * 100 << "%)\n";

                //std::string result_file = "orientations_" + loading_type + "_" + std::to_string(static_cast<int>(cp * 100)) + "pct.txt";
                //export_orientations(rve, result_file, acc, loading_type);
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
    out.close();

    return 0;
}
