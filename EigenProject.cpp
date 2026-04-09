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
	const std::string config_file = "rve_validation.txt";
	const std::string loading_type = "========= UNAXIAL TENSION =========";
	const double target_strain = 1.0;
	const double dt = 1e-4;
	std::vector<double> checkpoints = { 0.3, 0.5, 1.0 };
	double eps_dot = 0.1;
	Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
	L(0, 0) = -0.5 * eps_dot;
	L(1, 1) = -0.5 * eps_dot;
	L(2, 2) = eps_dot;


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

    return 0;
}
