#include <iostream>
#include <Eigen/Dense>
#include <string>
#include "Grain.h"

#include "trajectory_generator.hpp"

using namespace Eigen;

class Model
{
private:
	unsigned counter;
	Eigen::Matrix3d checking;


public:
	Model() {
		counter = 0;
	}

	void set_L(Eigen::Matrix3d L)
	{
		L = 0.5*(L + L.transpose());
		L = L - 1.0 / 3.0 * L.trace()*MatrixXd::Identity(3, 3);
		checking = L;
		//std::cout << sqrt((L*L.transpose()).trace()) << std::endl;
		//std::cout << L(0,0)<<"	" << L(1, 0) << "	" << L(2, 0) << "	" << std::endl;
		//std::cout << L.norm() << "	" << std::endl;
	}
	void evolve()
	{
		counter++;
	}
	void output_semax()
	{
	}
	void write_counter()
	{
		std::cout << counter << std::endl;
	}
};

int fmain()
{
	double mean_ds = 0;
	double s = 0;

	Matrix3d F = MatrixXd::Identity(3, 3);
	Model model;

	unsigned N = 20000;
	unsigned i = N-1;
	double dt = 5*1e-5;
	// ----------------------------------------------------------
	// 1. Определяем тип траектории и seed
	// ----------------------------------------------------------
	unsigned int rw_count = (N * 2) / 3;  // первые 2/3 — RW, остальные — SM
	bool         is_rw = (i < rw_count);
	unsigned int seed = i;             // уникальный seed для каждого эксперимента

	// ----------------------------------------------------------
	// 2. Генерируем последовательность тензоров деформации ε
	//    RW: 1001 точка (включая нулевую)
	//    SM: 1000 точек
	// ----------------------------------------------------------
	std::vector<Matrix3d> eps_seq = is_rw
		? generateRW(seed, /*n_final=*/1001)
		: generateSM(seed, /*n_final=*/1000);

	// ----------------------------------------------------------
	// 3. Переписанный while-цикл
	//
	//    Оригинал:
	//      while (s_sum < s_max) {
	//          update_L();
	//          for (unsigned k = 0; k < n; k++) model.evolve();
	//          output_semax();
	//      }
	//
	//    Теперь вместо s_sum < s_max — итерация по интервалам траектории.
	//    Вместо update_L() — вычисление L_norm из соседних ε и вызов set_L.
	//    n шагов определяется через длину интервала и dt.
	// ----------------------------------------------------------
	for (size_t seg = 0; seg + 1 < eps_seq.size(); ++seg)
	{
		// --- Аналог: update_L() ---
		//
		// Вычисляем нормированный градиент скорости L и число шагов n
		// для перехода от eps_seq[seg] к eps_seq[seg+1].
		//
		// L_norm * dt  →  за n шагов модель пройдёт ~s_int деформации.
		// Остаток (s_int mod dt) < dt — недохождение меньше одного шага.
		IntervalData iv = computeInterval(eps_seq[seg], eps_seq[seg + 1], dt);


		//std::cout << dt * iv.n_steps << std::endl;
		mean_ds += (dt * iv.n_steps)/double(eps_seq.size());
		s += dt * iv.n_steps;
		std::cout << iv.n_steps << std::endl;
		if (iv.n_steps == 0)
			continue;  // интервал короче одного шага dt — пропускаем

		model.set_L(iv.L_norm);  // <-- вместо update_L()

		// --- Аналог: for (k < n) model.evolve() ---
		for (int k = 0; k < iv.n_steps; ++k)
			model.evolve();

		//F = ((iv.L_norm)*iv.n_steps*dt).exp()*F;
		//Matrix3d U = F.transpose()*F;
		//Matrix3d Deformation = 0.5*U.log();
		//Matrix3d check1 = Deformation.exp();
		//Matrix3d check2 = check1.log();

		//std::cout << eps_seq[seg + 1](0, 0) - Deformation(0, 0) << std::endl;
		//std::cout << check2(0, 0) - Deformation(0, 0) << std::endl;
		// --- Аналог: output_semax() ---
		model.output_semax();
	}





	//std::string info  = "grain.txt";
	//Grain gr(info);
	//
	//double de = 0.001;
	//Eigen::Matrix3d L;
	//L(0, 0) = de;
	//L(0, 1) = 0;
	//L(0, 2) = 0;
	//
	//L(1, 0) = 0;
	//L(1, 1) = -0.5*de;
	//L(1, 2) = 0;
	//
	//L(2, 0) = 0;
	//L(2, 1) = 0;
	//L(2, 2) = -0.5*de;
	//
	//Eigen::VectorXd Stress;
	//
	//double dt = 1e-2;
	//for (unsigned i = 0; i < 1000; i++)
	//{
	//	if (i == 20)
	//		std::cout << "here" << '\n';
	//	gr.Simulate_Step(dt, L);
	//	Stress = gr.get_stress();
	//	std::cout << Stress(0)/1e6 << '\n';
	//	
	//}

	model.write_counter();
	std::cout << s << std::endl;
	std::cout << mean_ds << std::endl;

	return 0;
}
