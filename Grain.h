#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/CXX11/Tensor>
#include <string>
#include <memory>
#include "MaterialParams.h"

class Grain
{
	using Vector6d = Eigen::Matrix<double, 6, 1>;
	using Matrix6d = Eigen::Matrix<double, 6, 6>;
	using Vector9d = Eigen::Matrix<double, 9, 1>;

	MaterialParamsPtr params;

	unsigned rot_mode;

	Eigen::MatrixXd H;

	Eigen::Matrix3d Buf;
	Vector6d Buf_v;

	Eigen::Matrix3d W;
	Eigen::Matrix3d Spin;
	Eigen::Matrix3d Wp;
	Eigen::Matrix3d Z;
	Eigen::Matrix3d Ze;
	Eigen::Matrix3d Zp;
	Vector6d Dp_v;
	Eigen::Matrix3d L;
	Vector6d D_v;
	Vector9d L_9;
	Vector9d Z_9;
	Vector9d Ze_9;
	Vector9d Zp_9;
	Vector9d W_9;
	Eigen::Matrix3d Deformation;
	Eigen::Matrix3d F;

	Eigen::Matrix3d Stress;
	Vector6d Stress_v;
	Eigen::Matrix3d Orient;

	Eigen::Quaterniond orient_q;
	mutable Eigen::Matrix3d R_buf;
	mutable Eigen::Vector3d omega_buf;
	mutable Eigen::Quaterniond dq_buf;
	mutable Eigen::Matrix3d W_mat_buf;
	mutable Eigen::Matrix3d Omega_buf;
	mutable Eigen::Matrix3d Wp_mat_buf;
	mutable Eigen::Matrix3d S_buf;

	std::vector<bool> last_step;
	std::vector<bool> new_step;

	long int random_state;
	long int num_crystall;

	std::vector<double> tc;
	std::vector<double> tc_copy;

	std::vector<double> shifts;
	std::vector<double> shifts_copy;
	std::vector<double> sp_shifts;

	int linear_mode;
	double sum_shifts;

	double acc_def, time, stress_inten;
	double volume = 1.0;

	void clear();
	bool check_elastic(double dt);

	double dT;
	double Sum_time;
	void find_SS(double dT);
	void integrate(double dt);

	void computePlasticSpin();
	void updateOrientation(const Eigen::Matrix3d& Omega, double dt);

public:
	Grain(std::string info);
	Grain(MaterialParamsPtr p, const Eigen::Quaterniond& orientation, unsigned rot_mode);
	void Simulate_Step(double dt, Eigen::Matrix3d Lnew);

	Grain::Vector6d get_stress();

	Eigen::Quaterniond get_orientation() const;

	void set_orientation(const Eigen::Quaterniond& q);

	Eigen::Matrix3d get_stress_lab() const;

	double get_volume() const { return volume; }
	void set_volume(double v) { volume = v; }
};
