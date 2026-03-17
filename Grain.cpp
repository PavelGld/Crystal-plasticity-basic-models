#include "Grain.h"
#include <iostream>
#include <fstream>
#include <unsupported/Eigen/MatrixFunctions>

// Низкоуровневое сопоставление индексов (как в исходной реализации)
static int tensorToVoigtIndex(int i, int j) {
	static const int voigt_index[3][3] = {
		{0, 5, 4},  // (0,0), (0,1)/(1,0), (0,2)/(2,0)
		{5, 1, 3},
		{4, 3, 2}
	};
	return voigt_index[i][j];
}

static std::pair<int, int> voigtToTensorIndex(int voigt_idx) {
	static const std::pair<int, int> tensor_index[6] = {
		{0, 0},
		{1, 1},
		{2, 2},
		{1, 2},
		{0, 2},
		{0, 1}
	};
	return tensor_index[voigt_idx];
}

// Преобразование 3x3 → Vector6d: для симметричных тензоров
// вне диагонали получается инженерная сдвиговая компонента (сумма пар).
static Eigen::Matrix<double, 6, 1> tensorToVoigt(const Eigen::Matrix3d& A) {
	Eigen::Matrix<double, 6, 1> v;
	v.setZero();
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			v(tensorToVoigtIndex(i, j)) += A(i, j);
	return v;
}

// Обратное преобразование Vector6d → симметричный 3x3
static Eigen::Matrix3d voigtToTensor(const Eigen::Matrix<double, 6, 1>& v) {
	Eigen::Matrix3d A;
	A.setZero();
	for (int k = 0; k < 6; ++k) {
		auto ij = voigtToTensorIndex(k);
		int i = ij.first;
		int j = ij.second;
		if (i == j) {
			A(i, j) = v(k);
		} else {
			A(i, j) += 0.5 * v(k);
			A(j, i) += 0.5 * v(k);
		}
	}
	return A;
}

// Обратное преобразование Vector6d → симметричный 3x3 (in-place, без аллокаций)
static void voigtToTensorInPlace(const Eigen::Matrix<double, 6, 1>& v, Eigen::Matrix3d& A) {
	A.setZero();
	for (int k = 0; k < 6; ++k) {
		auto ij = voigtToTensorIndex(k);
		int i = ij.first;
		int j = ij.second;
		if (i == j) {
			A(i, j) = v(k);
		} else {
			A(i, j) += 0.5 * v(k);
			A(j, i) += 0.5 * v(k);
		}
	}
}

// Преобразование 3x3 → Vector9d (строчно-мажорная развёртка)
static Eigen::Matrix<double, 9, 1> tensorTo9d(const Eigen::Matrix3d& A) {
	Eigen::Matrix<double, 9, 1> v;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			v(3 * i + j) = A(i, j);
	return v;
}

// Обратное преобразование Vector9d → 3x3
static Eigen::Matrix3d nineDToTensor(const Eigen::Matrix<double, 9, 1>& v) {
	Eigen::Matrix3d A;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			A(i, j) = v(3 * i + j);
	return A;
}

Grain::Grain(std::string info)
{


	std::ifstream file(info);
	file >> hard_mode >> rot_mode;
	file >> num_slyde_sys;

	std::vector<double> b = { 0,0,0 };
	std::vector<double> n = { 0,0,0 };

	//M Wpm;
	M.resize(num_slyde_sys * 2);
	Wpm.resize(num_slyde_sys * 2);
	M0.resize(num_slyde_sys * 2);
	Wpm0.resize(num_slyde_sys * 2);
	Ms_v.resize(num_slyde_sys * 2);
	M_9.resize(num_slyde_sys * 2);
	Wpm_9.resize(num_slyde_sys * 2);

	for (int i = 0; i < num_slyde_sys; i++)
	{
		double norm = 0;
		for (int j = 0; j < 3; j++)
			file >> b[j];
	
	
	
		for (int j = 0; j < 3; j++)
			file >> n[j];
	
		if (n[0] * b[0] + n[1] * b[1] + n[2] * b[2] != 0)
			std::cout << "NORMIRATION FALSE" << '\n';
	
		norm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		for (int j = 0; j < 3; j++)
			n[j] = n[j] / norm;
	
		norm = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
		for (int j = 0; j < 3; j++)
			b[j] = b[j] / norm;

		M[i * 2].setZero();
		Wpm[i * 2].setZero();

		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
			{
				M[i * 2](j, k) = b[k]*n[j];
				Wpm[i * 2](j, k) = 0.5*(b[k] * n[j]- b[j] * n[k]);
			}
		Eigen::Matrix3d Ms = 0.5 * (M[i * 2] + M[i * 2].transpose());
		Ms_v[i * 2] = tensorToVoigt(Ms);
		M_9[i * 2]   = tensorTo9d(M[i * 2]);
		Wpm_9[i * 2] = tensorTo9d(Wpm[i * 2]);
		// Сохраняем начальные тензоры в КСК (постоянны на протяжении всего расчёта)
		M0[i * 2]   = M[i * 2];
		Wpm0[i * 2] = Wpm[i * 2];
	
		for (int j = 0; j < 3; j++)
			b[j] = -1 * b[j];

		M[i * 2 + 1].setZero();
		Wpm[i * 2 + 1].setZero();

		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
			{
				M[i * 2 + 1](j, k) += b[k] * n[j];
				Wpm[i * 2 + 1](j, k) += 0.5*(b[k] * n[j] - b[j] * n[k]);
			}
		Eigen::Matrix3d Ms2 = 0.5 * (M[i * 2 + 1] + M[i * 2 + 1].transpose());
		Ms_v[i * 2 + 1] = tensorToVoigt(Ms2);
		M_9[i * 2 + 1]   = tensorTo9d(M[i * 2 + 1]);
		Wpm_9[i * 2 + 1] = tensorTo9d(Wpm[i * 2 + 1]);
		// Сохраняем начальные тензоры в КСК (постоянны на протяжении всего расчёта)
		M0[i * 2 + 1]   = M[i * 2 + 1];
		Wpm0[i * 2 + 1] = Wpm[i * 2 + 1];
	
	
	
	}
	
	double c1, c2, c3;
	file >> c1 >> c2 >> c3;

	P = Eigen::Tensor<double, 4>(3, 3, 3, 3);

	for (int j = 0; j < 3; j++)
		for (int k = 0; k < 3; k++)
			for (int g = 0; g < 3; g++)
				for (int f = 0; f < 3; f++)
					P(j,k,g,f) = 0;

	P(0, 0, 0, 0) = P(1, 1, 1, 1) = P(2, 2, 2, 2) = c1;
	P(0, 0, 1, 1) = P(1, 1, 0, 0) = P(1, 1, 2, 2) = P(2, 2, 1, 1) = P(2, 2, 0, 0) = P(0, 0, 2, 2) = c2;
	P(0, 1, 0, 1) = P(1, 0, 0, 1) = P(0, 1, 1, 0) = P(1, 0, 1, 0) = P(1, 2, 1, 2) = P(1, 2, 2, 1) = P(2, 1, 1, 2) = P(2, 1, 2, 1) = P(2, 0, 2, 0) = P(0, 2, 2, 0) = P(2, 0, 0, 2) = P(0, 2, 0, 2) = c3;
	
	P_v.setZero();

	for (unsigned i = 0; i < 6; i++)
	{
		for (unsigned j = 0; j < 6; j++)
		{
			std::pair<int, int> ind1 = voigtToTensorIndex(i);
			std::pair<int, int> ind2 = voigtToTensorIndex(j);
			P_v(i, j) = P(ind1.first, ind1.second, ind2.first, ind2.second);
		}
	}

	H.resize(num_slyde_sys * 2, num_slyde_sys * 2);
	D_v.setZero();
	Dp_v.setZero();
	Stress_v.setZero();

	if (hard_mode == 2)
	{
		file >> A >> q >> h0 >> ts >> alfa;
	
		for (int i = 0; i < num_slyde_sys * 2; i++)
			for (int j = 0; j < num_slyde_sys * 2; j++)
				if (i == j)
					H(j,i) = A;
				else
					H(j, i) = A * q;
		h0 = h0;	}
	
	file >> tc0;

	file >> strain_rate_infl;
	
	time = 0;
	acc_def = 0;
	stress_inten = 0;

	Buf.setZero();// тензор на случай необходимости доп памяти

	W.setZero();
	W_9 = tensorTo9d(W);
	Spin.setZero();
	Wp.setZero();
	Z.setZero();
	Z_9  = tensorTo9d(Z);
	Ze.setZero();
	Ze_9 = tensorTo9d(Ze);
	Zp.setZero();
	Zp_9 = tensorTo9d(Zp);
	L.setZero();
	L_9  = tensorTo9d(L);
	Deformation.setZero();
	F.setIdentity();

	Stress.setZero();
	Orient.setZero();
	orient_q = Eigen::Quaterniond::Identity();
	R_buf.setIdentity();
	omega_buf.setZero();
	dq_buf = Eigen::Quaterniond::Identity();
	W_mat_buf.setZero();
	Omega_buf.setZero();
	Wp_mat_buf.setZero();
	S_buf.setZero();

	last_step.resize(num_slyde_sys * 2);
	new_step.resize(num_slyde_sys * 2);

	tc.resize(num_slyde_sys * 2);
	tc_copy.resize(num_slyde_sys * 2);

	shifts.resize(num_slyde_sys * 2);
	shifts_copy.resize(num_slyde_sys * 2);
	sp_shifts.resize(num_slyde_sys * 2);

	for (int i = 0; i < num_slyde_sys*2; i++)
	{
		last_step[i] = 0;
		new_step[i] = 0;

		tc[i] = tc0;
		tc_copy[i] = tc0;

		shifts[i] = 0;
		shifts_copy[i] = 0;
		sp_shifts[i] = 0;
	}

}

void Grain::clear()
{
	Buf.setZero();
	Buf_v.setZero();
	W.setZero();
	W_9 = tensorTo9d(W);
	Wp.setZero();
	Z.setZero();
	Z_9  = tensorTo9d(Z);
	Ze.setZero();
	Ze_9 = tensorTo9d(Ze);
	Zp.setZero();
	Zp_9 = tensorTo9d(Zp);
	Dp_v.setZero();
	D_v.setZero();
	L_9 = tensorTo9d(L);
	for (unsigned i = 0; i < 2 * num_slyde_sys; i++)
	{
		last_step[i] = new_step[i];
		new_step[i] = 0;
		//sp_shifts[i] = 0;
	}
}

bool Grain::check_elastic(double dt)
{

	bool TF = true;

	Buf_v = Stress_v + P_v * D_v * dt;

	for (unsigned i = 0; i < 2 * num_slyde_sys; i++)
	{
		double t = Buf_v.dot(Ms_v[i]);
		if (t > tc[i])
			TF = false;
	}

	if (!TF)
	{
		bool was_in_plastic(false);

		for (unsigned i = 0; i < 2 * num_slyde_sys; i++)
		{
			double t = Stress_v.dot(Ms_v[i]);
			if (t > tc[i])
				was_in_plastic = true;
		}


		if (was_in_plastic)
		{
			Buf_v.setZero();
			return TF;
		}
		else
		{
			double buft(1e300);
			for (unsigned i = 0; i < 2 * num_slyde_sys; i++)
			{
				double ai = (P_v*D_v).dot(Ms_v[i]);
				double t = Stress_v.dot(Ms_v[i]);
				if (ai > 0)
				{
					if ((tc[i] - t) / ai < buft)
						buft = (tc[i] - t) / ai;
				}
			}
			Sum_time += buft;
			dT -= buft;
			Stress_v += P_v * D_v * dT;
			Buf_v.setZero();
			integrate(dT);
			return false;
		}
	}

	// Упругий шаг: sp_shifts и new_step должны быть обнулены,
	// чтобы integrate() не использовал мусорные значения предыдущего шага
	for (unsigned i = 0; i < 2 * num_slyde_sys; i++)
	{
		sp_shifts[i] = 0.0;
		new_step[i] = 0;
	}
	Stress_v = Buf_v;
	Buf_v.setZero();
	integrate(dt);
	return true;
}

void Grain::find_SS(double dT)
{ 

	// Интенсивность скорости деформаций по J2 в 6D-нотации Фойгта
	double trD = D_v(0) + D_v(1) + D_v(2);
	double mean = trD / 3.0;
	double d1 = D_v(0) - mean;
	double d2 = D_v(1) - mean;
	double d3 = D_v(2) - mean;
	double j2 = d1 * d1 + d2 * d2 + d3 * d3
		+ 0.5 * (D_v(3) * D_v(3) + D_v(4) * D_v(4) + D_v(5) * D_v(5));
	double du = std::sqrt(2.0 / 3.0 * j2);

	for (unsigned i = 0; i < 2 * num_slyde_sys; i++)
		if (last_step[i])
			Dp_v += Ms_v[i] * sp_shifts[i];

	Buf_v = Stress_v + P_v * (D_v - Dp_v) * dT;

	for (unsigned i = 0; i < 2 * num_slyde_sys; i++)
	{
		double t = Buf_v.dot(Ms_v[i]);
		if (t > tc[i])
		{
			new_step[i] = 1;
			sp_shifts[i] = strain_rate_infl *du*pow(fabs(t / tc[i]), 100.0);
		}
		else
		{
			new_step[i] = 0;
			sp_shifts[i] = 0;
		}

	}

	Stress_v = Buf_v;

	integrate(dT);
}

void Grain::integrate(double dt)
{
	for (unsigned i = 0; i < 2 * num_slyde_sys; i++)
		for (unsigned j = 0; j < 2 * num_slyde_sys; j++)
			if ((new_step[j]) &&(tc[i]<ts))
				tc[i] += H(i, j) * sp_shifts[j];

	if (hard_mode == 2)
		for (int j = 0; j < num_slyde_sys * 2; j++)
		{
			double t_j = tc[j];
			double Hii = A * h0*pow(fabs(1 - t_j / ts), alfa);

			for (int i = 0; i < num_slyde_sys * 2; i++)
			{
				if (i == j)
					H(i, j) = Hii;
				else
					H(i, j) = Hii * q;
			}
			if (Hii != Hii)
				std::cout << "here";
		}

	// Обновляем ориентацию зерна в конце каждого подшага
	if (rot_mode != 0)
	{
		W_mat_buf = 0.5 * (L - L.transpose());

		if (rot_mode == 1)
		{
			Omega_buf = W_mat_buf;
		}
		else if (rot_mode == 2)
		{
			computePlasticSpin();
			Omega_buf.noalias() = W_mat_buf - Wp_mat_buf;
		}

		updateOrientation(Omega_buf, dt);
	}
}

void Grain::computePlasticSpin()
{
	Wp_mat_buf.setZero();
	for (unsigned i = 0; i < 2 * num_slyde_sys; i++)
		if (new_step[i])
			Wp_mat_buf.noalias() += sp_shifts[i] * Wpm0[i];
}

void Grain::updateOrientation(const Eigen::Matrix3d& Omega, double dt)
{
	omega_buf(0) = Omega(2, 1) - Omega(1, 2);
	omega_buf(1) = Omega(0, 2) - Omega(2, 0);
	omega_buf(2) = Omega(1, 0) - Omega(0, 1);
	omega_buf *= 0.5;

	double angle = omega_buf.norm() * dt;
	if (angle > 1e-14)
	{
		dq_buf = Eigen::AngleAxisd(angle, omega_buf / omega_buf.norm());
		orient_q = dq_buf * orient_q;
		orient_q.normalize();
		R_buf = orient_q.toRotationMatrix();
	}
}

void Grain::Simulate_Step(double dt, Eigen::Matrix3d Lnew)
{
	clear();

	if (rot_mode != 0)
	{
		L.noalias() = R_buf.transpose() * Lnew * R_buf;
	}
	else
	{
		L = Lnew;
	}

	L_9 = tensorTo9d(L);
	Buf = 0.5 * (L + L.transpose());
	D_v = tensorToVoigt(Buf);

	dT = 0;
	Sum_time = 0;

	dT = dt - Sum_time;

	while (dt > Sum_time)
	{
		if (check_elastic(dT))
			Sum_time += dT;
		else
		{
			find_SS(dT);
			Sum_time += dT;
		}
		dT = dt - Sum_time;
	}	

	F = (dt*L).exp() * F;
	Deformation = 0.5 * (F.transpose() * F - Eigen::Matrix3d::Identity());

	// При желании можно вызывать здесь пересчёт всех 9D-представлений:
	// Update9D();

}

Eigen::Matrix<double, 6, 1> Grain::get_stress()
{
	return Stress_v;
}

Eigen::Quaterniond Grain::get_orientation() const
{
	return orient_q;
}

Eigen::Matrix3d Grain::get_stress_lab() const
{
	voigtToTensorInPlace(Stress_v, S_buf);
	return R_buf * S_buf * R_buf.transpose();
}