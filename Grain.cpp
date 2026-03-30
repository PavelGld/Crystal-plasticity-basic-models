#include "Grain.h"
#include <iostream>
#include <fstream>
#include <unsupported/Eigen/MatrixFunctions>

static int tensorToVoigtIndex(int i, int j) {
	static const int voigt_index[3][3] = {
		{0, 5, 4},
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

static Eigen::Matrix<double, 6, 1> tensorToVoigt(const Eigen::Matrix3d& A) {
	Eigen::Matrix<double, 6, 1> v;
	v.setZero();
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			v(tensorToVoigtIndex(i, j)) += A(i, j);
	return v;
}

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

static Eigen::Matrix<double, 9, 1> tensorTo9d(const Eigen::Matrix3d& A) {
	Eigen::Matrix<double, 9, 1> v;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			v(3 * i + j) = A(i, j);
	return v;
}

Grain::Grain(std::string info)
{
	params = std::make_shared<MaterialParams>(info);

	std::ifstream file(info);
	unsigned hard_mode_temp;
	file >> hard_mode_temp >> rot_mode;
	file.close();

	H = params->H;

	time = 0;
	acc_def = 0;
	stress_inten = 0;

	Buf.setZero();

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

	D_v.setZero();
	Dp_v.setZero();
	Stress_v.setZero();

	last_step.resize(params->num_slyde_sys * 2);
	new_step.resize(params->num_slyde_sys * 2);

	tc.resize(params->num_slyde_sys * 2);
	tc_copy.resize(params->num_slyde_sys * 2);

	shifts.resize(params->num_slyde_sys * 2);
	shifts_copy.resize(params->num_slyde_sys * 2);
	sp_shifts.resize(params->num_slyde_sys * 2);

	for (int i = 0; i < params->num_slyde_sys * 2; i++)
	{
		last_step[i] = 0;
		new_step[i] = 0;

		tc[i] = params->tc0;
		tc_copy[i] = params->tc0;

		shifts[i] = 0;
		shifts_copy[i] = 0;
		sp_shifts[i] = 0;
	}
}

Grain::Grain(MaterialParamsPtr p, const Eigen::Quaterniond& orientation, unsigned rot_mode)
	: params(p), rot_mode(rot_mode), H(p->H)
{
	time = 0;
	acc_def = 0;
	stress_inten = 0;

	Buf.setZero();

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
	orient_q = orientation.normalized();
	R_buf = orient_q.toRotationMatrix();
	omega_buf.setZero();
	dq_buf = Eigen::Quaterniond::Identity();
	W_mat_buf.setZero();
	Omega_buf.setZero();
	Wp_mat_buf.setZero();
	S_buf.setZero();

	D_v.setZero();
	Dp_v.setZero();
	Stress_v.setZero();

	last_step.resize(params->num_slyde_sys * 2);
	new_step.resize(params->num_slyde_sys * 2);

	tc.resize(params->num_slyde_sys * 2);
	tc_copy.resize(params->num_slyde_sys * 2);

	shifts.resize(params->num_slyde_sys * 2);
	shifts_copy.resize(params->num_slyde_sys * 2);
	sp_shifts.resize(params->num_slyde_sys * 2);

	for (int i = 0; i < params->num_slyde_sys * 2; i++)
	{
		last_step[i] = 0;
		new_step[i] = 0;

		tc[i] = params->tc0;
		tc_copy[i] = params->tc0;

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
	for (int i = 0; i < 2 * params->num_slyde_sys; i++)
	{
		last_step[i] = new_step[i];
		new_step[i] = 0;
	}
}

bool Grain::check_elastic(double dt)
{

	bool TF = true;

	Buf_v = Stress_v + params->P_v * D_v * dt;

	for (int i = 0; i < 2 * params->num_slyde_sys; i++)
	{
		double t = Buf_v.dot(params->Ms_v[i]);
		if (t > tc[i])
			TF = false;
	}

	if (!TF)
	{
		bool was_in_plastic(false);

		for (int i = 0; i < 2 * params->num_slyde_sys; i++)
		{
			double t = Stress_v.dot(params->Ms_v[i]);
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
			for (int i = 0; i < 2 * params->num_slyde_sys; i++)
			{
				double ai = (params->P_v*D_v).dot(params->Ms_v[i]);
				double t = Stress_v.dot(params->Ms_v[i]);
				if (ai > 0)
				{
					if ((tc[i] - t) / ai < buft)
						buft = (tc[i] - t) / ai;
				}
			}
			Sum_time += buft;
			dT -= buft;
			Stress_v += params->P_v * D_v * dT;
			Buf_v.setZero();
			integrate(dT);
			return false;
		}
	}

	for (int i = 0; i < 2 * params->num_slyde_sys; i++)
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

	double trD = D_v(0) + D_v(1) + D_v(2);
	double mean = trD / 3.0;
	double d1 = D_v(0) - mean;
	double d2 = D_v(1) - mean;
	double d3 = D_v(2) - mean;
	double j2 = d1 * d1 + d2 * d2 + d3 * d3
		+ 0.5 * (D_v(3) * D_v(3) + D_v(4) * D_v(4) + D_v(5) * D_v(5));
	double du = std::sqrt(2.0 / 3.0 * j2);

	for (int i = 0; i < 2 * params->num_slyde_sys; i++)
		if (last_step[i])
			Dp_v += params->Ms_v[i] * sp_shifts[i];

	Buf_v = Stress_v + params->P_v * (D_v - Dp_v) * dT;

	for (int i = 0; i < 2 * params->num_slyde_sys; i++)
	{
		double t = Buf_v.dot(params->Ms_v[i]);
		if (t > tc[i])
		{
			new_step[i] = 1;
			sp_shifts[i] = params->strain_rate_infl *du*pow(fabs(t / tc[i]), 100.0);
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
	for (int i = 0; i < 2 * params->num_slyde_sys; i++)
		for (int j = 0; j < 2 * params->num_slyde_sys; j++)
			if ((new_step[j]) &&(tc[i]<params->ts))
				tc[i] += H(i, j) * sp_shifts[j];

	if (params->hard_mode == 2)
		for (int j = 0; j < params->num_slyde_sys * 2; j++)
		{
			double t_j = tc[j];
			double Hii = params->A * params->h0*pow(fabs(1 - t_j / params->ts), params->alfa);

			for (int i = 0; i < params->num_slyde_sys * 2; i++)
			{
				if (i == j)
					H(i, j) = Hii;
				else
					H(i, j) = Hii * params->q;
			}
			if (Hii != Hii)
				std::cout << "here";
		}

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
	for (int i = 0; i < 2 * params->num_slyde_sys; i++)
		if (new_step[i])
			Wp_mat_buf.noalias() += sp_shifts[i] * params->Wpm[i];
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
		orient_q = orient_q * dq_buf;
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

}

Eigen::Matrix<double, 6, 1> Grain::get_stress()
{
	return Stress_v;
}

Eigen::Quaterniond Grain::get_orientation() const
{
	return orient_q;
}

void Grain::set_orientation(const Eigen::Quaterniond& q)
{
	orient_q = q.normalized();
	R_buf = orient_q.toRotationMatrix();
}

Eigen::Matrix3d Grain::get_stress_lab() const
{
	voigtToTensorInPlace(Stress_v, S_buf);
	return R_buf * S_buf * R_buf.transpose();
}
