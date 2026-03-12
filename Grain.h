#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <string>

class Grain
{
	unsigned hard_mode, rot_mode;

	Eigen::Tensor<double, 4> P; // (3, 3, 3, 3);
	Eigen::MatrixXd P_v;

	std::vector<Eigen::Matrix3d> M;
	std::vector<Eigen::Matrix3d> Wpm;

	std::vector<Eigen::VectorXd> M_v;

	Eigen::MatrixXd H;
	double degree, factor;// параметры в степенном изотропном упрочнении

	double A, q, h0, ts, alfa;// параметры в квазилинейном упрочнении

	Eigen::Matrix3d Buf;// тензор на случай необходимости доп памяти
	Eigen::VectorXd Buf_v;

	Eigen::Matrix3d W;
	Eigen::Matrix3d Spin;
	Eigen::Matrix3d Wp;
	Eigen::Matrix3d Z;
	Eigen::Matrix3d Ze;
	Eigen::Matrix3d Zp;
	Eigen::VectorXd Zp_v;
	Eigen::Matrix3d L;
	Eigen::VectorXd L_v;
	Eigen::Matrix3d Deformation;
	Eigen::Matrix3d F;

	Eigen::Matrix3d Stress;
	Eigen::VectorXd Stress_v;
	Eigen::Matrix3d Orient;

	std::vector<bool> last_step;
	std::vector<bool> new_step;

	long int random_state;// сид начальной генерации
	long int num_crystall;// номер кристалла

	std::vector<double> tc;
	std::vector<double> tc_copy;
	double tc0;

	std::vector<double> shifts;
	std::vector<double> shifts_copy;
	std::vector<double> sp_shifts;

	int num_slyde_sys;
	int linear_mode;// режим моделирования
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
	//bool Flexibility_test(double dT);// Проверка на описание шага упругим. Если его можно описать упругим, то вернет 1 и сразу посчитает упруго, иначе отметит СС бывшие активными ранее и активные на данном шаге в new_step
	//int checking_old_SS(double dT);// Проверка можно ли описать шаг с активными с прошлого шага СС, если да, то вернет "-1" и смоделирует с ними, иначе вернет номер СС которую нужно активировать (СС с наиб. напряжением) 
	//double find_dt(double sum_t, double dT, int num);// поиск подшага для активации заданной системы скольжения и совершение шага. Использует отмеченные активными СС. Возвращает пройденый шаг по времени
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
	//void Reset(double tc0);// обнуляет все свойства ПО и откатывает их к начальному состоянию
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

	Eigen::VectorXd get_stress();

};

