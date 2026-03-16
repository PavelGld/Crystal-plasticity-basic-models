#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <string>

class Grain
{
	using Vector6d = Eigen::Matrix<double, 6, 1>;
	using Matrix6d = Eigen::Matrix<double, 6, 6>;

	unsigned hard_mode, rot_mode;

	Eigen::Tensor<double, 4> P; // (3, 3, 3, 3);
	Matrix6d P_v;

	std::vector<Eigen::Matrix3d> M;
	std::vector<Eigen::Matrix3d> Wpm;

	std::vector<Vector6d> Ms_v;

	Eigen::MatrixXd H;
	double degree, factor;// ��������� � ��������� ���������� ����������

	double A, q, h0, ts, alfa;// ��������� � ������������� ����������

	Eigen::Matrix3d Buf;// ������ �� ������ ������������� ��� ������
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
	Eigen::Matrix3d Deformation;
	Eigen::Matrix3d F;

	Eigen::Matrix3d Stress;
	Vector6d Stress_v;
	Eigen::Matrix3d Orient;

	std::vector<bool> last_step;
	std::vector<bool> new_step;

	long int random_state;// ��� ��������� ���������
	long int num_crystall;// ����� ���������

	std::vector<double> tc;
	std::vector<double> tc_copy;
	double tc0;

	std::vector<double> shifts;
	std::vector<double> shifts_copy;
	std::vector<double> sp_shifts;

	int num_slyde_sys;
	int linear_mode;// ����� �������������
	double sum_shifts;

	double acc_def, time, stress_inten;

	void clear();
	bool check_elastic(double dt);

	double dT;
	double Sum_time;
	void find_SS(double dT);
	void integrate();

	//
	////Tensor find_We(vector<double> Sp_shifts);
	////Tensor find_Dp(vector<double> Sp_shifts);
	//
	//bool Flexibility_test(double dT);// �������� �� �������� ���� �������. ���� ��� ����� ������� �������, �� ������ 1 � ����� ��������� ������, ����� ������� �� ������ ��������� ����� � �������� �� ������ ���� � new_step
	//int checking_old_SS(double dT);// �������� ����� �� ������� ��� � ��������� � �������� ���� ��, ���� ��, �� ������ "-1" � ����������� � ����, ����� ������ ����� �� ������� ����� ������������ (�� � ����. �����������) 
	//double find_dt(double sum_t, double dT, int num);// ����� ������� ��� ��������� �������� ������� ���������� � ���������� ����. ���������� ���������� ��������� ��. ���������� ��������� ��� �� �������
	//
	//
	//void new_tc(double dT);
	//void new_shifts(double dT);
	//
	//void old_tc(double dT);
	//void old_shifts(double dT);
	//
	//void update_orient(double dt);
	//
	//double get_MPM(unsigned i, unsigned j);
	//double get_MPD(unsigned ind1);
	//
	//void Reset(double tc0);// �������� ��� �������� �� � ���������� �� � ���������� ���������
	//
	//void Set_H_Params(double E, double k, double newtc);
	//
	//void Set_const_model(double E, double newtc);
	//
	//void Set_st_Params(double factor, double degree, double newtc);

public:
	//Grain();
	//~Grain();
	Grain(std::string info);
	void Simulate_Step(double dt, Eigen::Matrix3d Lnew);

	Grain::Vector6d get_stress();

};

