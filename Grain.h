#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/CXX11/Tensor>
#include <string>

class Grain
{
	using Vector6d = Eigen::Matrix<double, 6, 1>;
	using Matrix6d = Eigen::Matrix<double, 6, 6>;
	using Vector9d = Eigen::Matrix<double, 9, 1>;

	unsigned hard_mode, rot_mode;

	Eigen::Tensor<double, 4> P; // (3, 3, 3, 3);
	Matrix6d P_v;

	std::vector<Eigen::Matrix3d> M;
	std::vector<Eigen::Matrix3d> Wpm;

	// Начальные тензоры систем скольжения в КСК (постоянны, не меняются)
	std::vector<Eigen::Matrix3d> M0;
	std::vector<Eigen::Matrix3d> Wpm0;

	std::vector<Vector6d> Ms_v;

	// Полные 9D-представления (строчно-мажорная развёртка 3x3)
	std::vector<Vector9d> M_9;
	std::vector<Vector9d> Wpm_9;

	Eigen::MatrixXd H;
	double degree, factor;// параметры степенного закона упрочнения

	double A, q, h0, ts, alfa;// параметры нелинейного упрочнения
	double strain_rate_infl;
	Eigen::Matrix3d Buf;// буфер на случай необходимости доп памяти
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
	// 9D-представления тензоров шага
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

	// Ориентация зерна: кватернион, задающий поворот КСК → ЛСК
	Eigen::Quaterniond orient_q;
	// Кэшированная матрица поворота (обновляется при изменении orient_q)
	mutable Eigen::Matrix3d R_buf;
	// Буфер для вектора угловой скорости
	mutable Eigen::Vector3d omega_buf;
	// Буфер для временного кватерниона при обновлении ориентации
	mutable Eigen::Quaterniond dq_buf;
	// Буферы для промежуточных вычислений в integrate()
	mutable Eigen::Matrix3d W_mat_buf;
	mutable Eigen::Matrix3d Omega_buf;
	mutable Eigen::Matrix3d Wp_mat_buf;
	// Буфер для тензора напряжений при переводе в ЛСК
	mutable Eigen::Matrix3d S_buf;

	std::vector<bool> last_step;
	std::vector<bool> new_step;

	long int random_state;// для генерации случайных
	long int num_crystall;// номер кристаллита

	std::vector<double> tc;
	std::vector<double> tc_copy;
	double tc0;

	std::vector<double> shifts;
	std::vector<double> shifts_copy;
	std::vector<double> sp_shifts;

	int num_slyde_sys;
	int linear_mode;// режим линеаризации
	double sum_shifts;

	double acc_def, time, stress_inten;

	void clear();
	bool check_elastic(double dt);

	double dT;
	double Sum_time;
	void find_SS(double dT);
	void integrate(double dt);

	// Вычисляет пластический спин Wp = sum_alpha (gamma_dot_alpha * Wpm[alpha])
	// Результат сохраняется в Wp_mat_buf
	void computePlasticSpin();

	// Обновляет кватернион ориентации по спину Omega за время dt
	void updateOrientation(const Eigen::Matrix3d& Omega, double dt);

public:
	Grain(std::string info);
	void Simulate_Step(double dt, Eigen::Matrix3d Lnew);

	Grain::Vector6d get_stress();

	// Возвращает текущий кватернион ориентации (КСК → ЛСК)
	Eigen::Quaterniond get_orientation() const;

	// Возвращает тензор напряжений Коши в лабораторной системе координат
	Eigen::Matrix3d get_stress_lab() const;
};

