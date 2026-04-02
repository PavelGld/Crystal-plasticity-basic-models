#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <string>
#include <cmath>
#include <fstream>
#include "Grain.h"
#include "RVE.h"

#define M_PI 3.14159265358979323846

void test_uniaxial_stretch_no_rotation()
{
	std::cout << "=== Test 1: Uniaxial stretch (no rotation expected) ===" << std::endl;

	Grain gr("grain_test.txt");

	// Начальная ориентация - поворот на 30 градусов вокруг Z
	Eigen::AngleAxisd aa_init(M_PI / 6.0, Eigen::Vector3d::UnitZ());
	Eigen::Quaterniond q_init(aa_init);
	gr.set_orientation(q_init);

	Eigen::Quaterniond q_before = gr.get_orientation();
	Eigen::Matrix3d R_before = q_before.toRotationMatrix();

	// L = diag(eps, -eps/2, -eps/2) - симметричный, W = 0
	double eps = 0.001;
	Eigen::Matrix3d L;
	L.setZero();
	L(0, 0) = eps;
	L(1, 1) = -0.5 * eps;
	L(2, 2) = -0.5 * eps;

	double dt = 1e-3;
	int N = 100;
	for (int i = 0; i < N; i++)
	{
		gr.Simulate_Step(dt, L);
	}

	Eigen::Quaterniond q_after = gr.get_orientation();
	Eigen::Matrix3d R_after = q_after.toRotationMatrix();

	// Проверяем что ориентация не изменилась
	Eigen::Matrix3d R_diff = R_after - R_before;
	double diff_norm = R_diff.norm();

	std::cout << "R_before:\n" << R_before << std::endl;
	std::cout << "R_after:\n" << R_after << std::endl;
	std::cout << "Difference norm: " << diff_norm << std::endl;

	if (diff_norm < 1e-10)
	{
		std::cout << "TEST 1 PASSED: Orientation unchanged (as expected)" << std::endl;
	}
	else
	{
		std::cout << "TEST 1 FAILED: Orientation changed unexpectedly!" << std::endl;
	}
	std::cout << std::endl;
}

void test_simple_shear_rotation()
{
	std::cout << "=== Test 2: Simple shear (rotation expected) ===" << std::endl;

	Grain gr("grain_test.txt");

	Eigen::Quaterniond q_before = gr.get_orientation();
	Eigen::Matrix3d R_before = q_before.toRotationMatrix();

	// L = [[0, eps, 0], [0, 0, 0], [0, 0, 0]] - простой сдвиг
	// W = 0.5*(L - L^T) = [[0, eps/2, 0], [-eps/2, 0, 0], [0, 0, 0]]
	// omega = (0, 0, -eps/2) - вращение вокруг Z (отрицательное = по часовой)
	// theta = |omega| * t = eps * t / 2
	double eps = 0.001;
	Eigen::Matrix3d L;
	L.setZero();
	L(0, 1) = eps;

	double dt = 1e-3;
	int N = 100;
	double total_time = N * dt;
	double expected_angle = eps * total_time / 2.0;

	for (int i = 0; i < N; i++)
	{
		gr.Simulate_Step(dt, L);
	}

	Eigen::Quaterniond q_after = gr.get_orientation();
	Eigen::Matrix3d R_after = q_after.toRotationMatrix();

	// Извлекаем угол поворота из кватерниона
	// q = [w, x, y, z] = [cos(theta/2), n*sin(theta/2)]
	// Для вращения вокруг Z: q = [cos(theta/2), 0, 0, ±sin(theta/2)]
	double actual_angle = 2.0 * std::acos(std::abs(q_after.w()));

	// Проверяем ось вращения - должна быть Z (или -Z)
	Eigen::Vector3d vec = q_after.vec();
	double vec_norm = vec.norm();
	Eigen::Vector3d axis = vec_norm > 1e-10 ? vec / vec_norm : Eigen::Vector3d(0, 0, 1);

	std::cout << "R_before:\n" << R_before << std::endl;
	std::cout << "R_after:\n" << R_after << std::endl;
	std::cout << "Expected rotation angle: " << expected_angle << " rad (" 
	          << expected_angle * 180.0 / M_PI << " deg)" << std::endl;
	std::cout << "Actual rotation angle: " << actual_angle << " rad ("
	          << actual_angle * 180.0 / M_PI << " deg)" << std::endl;
	std::cout << "Rotation axis: " << axis.transpose() << std::endl;
	std::cout << "Angle difference: " << std::abs(actual_angle - expected_angle) << std::endl;

	double angle_error = std::abs(actual_angle - expected_angle);
	if (angle_error < 1e-5)
	{
		std::cout << "TEST 2 PASSED: Rotation angle matches expected" << std::endl;
	}
	else
	{
		std::cout << "TEST 2 FAILED: Rotation angle mismatch!" << std::endl;
	}
	std::cout << std::endl;
}

void test_no_rotation_mode()
{
	std::cout << "=== Test 3: rot_mode=0 (no rotation) ===" << std::endl;

	// Зерно без ротации
	Grain gr_no_rot("grain_test_no_rot.txt");

	// Зерно с материальной ротацией
	Grain gr_rot("grain_test.txt");

	// Начальная ориентация - поворот на 45 градусов вокруг Z
	Eigen::AngleAxisd aa_init(M_PI / 4.0, Eigen::Vector3d::UnitZ());
	Eigen::Quaterniond q_init(aa_init);
	gr_no_rot.set_orientation(q_init);
	gr_rot.set_orientation(q_init);

	// Простой сдвиг
	double eps = 0.001;
	Eigen::Matrix3d L;
	L.setZero();
	L(0, 1) = eps;

	double dt = 1e-3;
	int N = 50;
	for (int i = 0; i < N; i++)
	{
		gr_no_rot.Simulate_Step(dt, L);
		gr_rot.Simulate_Step(dt, L);
	}

	Eigen::Matrix3d R_no_rot = gr_no_rot.get_orientation().toRotationMatrix();
	Eigen::Matrix3d R_rot = gr_rot.get_orientation().toRotationMatrix();
	Eigen::Matrix3d R_init = q_init.toRotationMatrix();

	std::cout << "Initial R:\n" << R_init << std::endl;
	std::cout << "R (rot_mode=0) after " << N << " steps:\n" << R_no_rot << std::endl;
	std::cout << "R (rot_mode=1) after " << N << " steps:\n" << R_rot << std::endl;

	double diff_no_rot = (R_no_rot - R_init).norm();
	double diff_rot = (R_rot - R_init).norm();

	std::cout << "Difference (no rotation): " << diff_no_rot << std::endl;
	std::cout << "Difference (with rotation): " << diff_rot << std::endl;

	if (diff_no_rot < 1e-10)
	{
		std::cout << "TEST 3 PASSED: rot_mode=0 keeps orientation fixed" << std::endl;
	}
	else
	{
		std::cout << "TEST 3 FAILED: rot_mode=0 changed orientation!" << std::endl;
	}
	std::cout << std::endl;
}

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

	std::ofstream out;          // поток для записи
	out.open("hello.txt");      // открываем файл для записи
	for (int i = 0; i < 1000; i++) {
		rve.Simulate_Step(dt, L);

		Eigen::Matrix3d stress = rve.get_stress();
		//std::cout << "Average stress:\n" << stress/1e6 << std::endl;
		if (i%100==0)
			out << dt * eps*i*100 << "	" << stress(0,0)/1e6<<"	" << stress(1, 1) / 1e6 << "	" << stress(2, 2) / 1e6 << "\n";
	}

	Eigen::Matrix3d stress = rve.get_stress();
	std::cout << "Average stress:\n" << stress << std::endl;
	std::cout << "Time: " << rve.get_time() << std::endl;
	
	if (rve.num_grains() == 1000) {
		std::cout << "TEST 4 PASSED: RVE created 10 grains" << std::endl;
	} else {
		std::cout << "TEST 4 FAILED: Wrong number of grains" << std::endl;
	}
	std::cout << std::endl;
}

void test_soft_loading()
{
	std::cout << "=== Test 5: Soft loading (uniaxial tension) ===" << std::endl;

	RVE rve("rve_test.txt");
	
	double stress_rate = 1e8;  // Па/с
	Eigen::Matrix3d dS_dt = Eigen::Matrix3d::Zero();
	dS_dt(0, 0) = stress_rate;
	
	double dt = 1e-5;
	int n_steps = 100;
	
	for (int i = 0; i < n_steps; i++) {
		rve.Simulate_Step_Soft(dt, dS_dt);
		Eigen::Matrix3d stress = rve.get_stress();
		std::cout << "\n" << stress << std::endl;
	}
	
	Eigen::Matrix3d stress = rve.get_stress();
	std::cout << "Final stress:\n" << stress << std::endl;
	std::cout << "Expected s11: " << stress_rate * n_steps * dt << " Pa" << std::endl;
	
	double s11 = std::abs(stress(0, 0));
	double s22 = std::abs(stress(1, 1));
	double s33 = std::abs(stress(2, 2));
	double s12 = std::abs(stress(0, 1));
	double s13 = std::abs(stress(0, 2));
	double s23 = std::abs(stress(1, 2));
	
	double off_diag_max = std::max({s12, s13, s23, s22, s33});
	
	double tol = 0.1;
	if (off_diag_max < tol * s11) {
		std::cout << "TEST 5 PASSED: Stress is approximately uniaxial" << std::endl;
	} else {
		std::cout << "TEST 5 FAILED: Off-diagonal components too large" << std::endl;
		std::cout << "  s11 = " << s11 << ", max_off_diag = " << off_diag_max << std::endl;
	}
	std::cout << std::endl;
}

int main()
{
	test_uniaxial_stretch_no_rotation();
	test_simple_shear_rotation();
	test_no_rotation_mode();
	test_rve_basic();
	test_soft_loading();

	std::cout << "=== All tests completed ===" << std::endl;
	return 0;
}
