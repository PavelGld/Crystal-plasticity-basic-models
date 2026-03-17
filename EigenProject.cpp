#include <iostream>
#include <Eigen/Dense>
#include <string>
#include "Grain.h"

using namespace Eigen;

int main()
{
	std::string info  = "grain.txt";
	Grain gr(info);
	
	double de = 0.001;
	Eigen::Matrix3d L;
	L(0, 0) = de;
	L(0, 1) = 0;
	L(0, 2) = 0;
	
	L(1, 0) = 0;
	L(1, 1) = -0.5*de;
	L(1, 2) = 0;
	
	L(2, 0) = 0;
	L(2, 1) = 0;
	L(2, 2) = -0.5*de;
	
	Eigen::VectorXd Stress;
	
	double dt = 1e-2;
	for (unsigned i = 0; i < 25000; i++)
	{
		if (i == 20)
			std::cout << "here" << '\n';
		gr.Simulate_Step(dt, L);
		Stress = gr.get_stress();
		std::cout << Stress(0)/1e6 << '\n';
		
	}

	return 0;
}
