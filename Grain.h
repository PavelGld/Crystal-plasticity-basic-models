#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <string>
#include <memory>
#include <cstdint>
#include "tensor_utils.h"
#include "MaterialParams.h"

class Grain
{

    MaterialParamsPtr params;

    Eigen::MatrixXd H;
    mutable Eigen::VectorXd H_buf;

    mutable Eigen::Matrix3d Buf;
    mutable Vector6d Buf_v;
    mutable Vector6d Buf_PD;
    mutable Matrix6d Buf_P;

    Eigen::Matrix3d Wp;
    Eigen::Matrix3d L;
    Vector6d Dp_v;
    Vector6d D_v;
    Eigen::Matrix3d strain_gl;  // Green-Lagrange strain tensor
    Eigen::Matrix3d F;

    Eigen::Matrix3d Stress;
    Vector6d Stress_v;

    Eigen::Quaterniond orient_q;
    mutable Eigen::Matrix3d R_buf;
    mutable Eigen::Vector3d omega_buf;
    mutable Eigen::Quaterniond dq_buf;
    mutable Eigen::Matrix3d W_mat_buf;
    mutable Eigen::Matrix3d Omega_buf;
    mutable Eigen::Matrix3d Wp_mat_buf;
    mutable Eigen::Matrix3d S_buf;

    std::vector<uint8_t> last_step;
    std::vector<uint8_t> new_step;

    long int random_state;  // seed for grain orientation generation

    std::vector<double> tc;
    std::vector<double> shifts;
    std::vector<double> sp_shifts;

    double acc_def, sim_time;  // grain-level state (used when modeling a single grain)
    double volume = 1.0;

    static constexpr double ROTATION_ANGLE_THRESHOLD = 1e-14;
    static constexpr double TIME_EPSILON = 1e-15;
    static constexpr double DS_CRITICAL = 1.7e-6;

    // Неявная схема (integ_mode==1)
    static constexpr double NEWTON_TOL = 1e-8;
    static constexpr int NEWTON_MAX_ITER = 20;
    static constexpr int NEWTON_MAX_BACKTRACKS = 5;
    static constexpr double RATIO_MAX = 2.0;          // кламп τ/τc (2^50 ≈ 1e15, переполнения нет)
    static constexpr int MAX_FALLBACK_DEPTH = 8;
    static constexpr double DS_CRITICAL_IMPL = 1e-4;  // порог подшагивания для неявной схемы

    mutable Matrix6d tangent_buf;   // консистентный тангенс последнего неявного шага (КСК)
    bool tangent_valid = false;

    mutable Matrix6d J_buf;
    mutable Matrix6d dyad_buf;
    Vector6d res_v;                 // невязка R
    Vector6d res_new_v;
    Vector6d sigma_trial_v;
    Vector6d sigma_new_v;
    Vector6d delta_sigma_v;
    Vector6d trial_inc_v;           // P_v * D_v (trial-инкремент)
    Eigen::VectorXd gamma_dot;      // скорости сдвига (итерат/сошедшиеся), n
    Eigen::VectorXd tau_buf;        // τ^α, n
    Eigen::VectorXd tc_eff;         // эффективные τc с учётом упрочнения внутри шага, n
    Eigen::MatrixXd G_buf;          // ∂γ̇/∂σ, n×6
    Eigen::MatrixXd G_prev;         // G с предыдущего sweep'а, n×6
    double rate_factor_impl = 0.0;  // γ̇₀ · D_eq на текущий шаг

    void solve_implicit(double dt, int depth);
    void evaluate_Dp(const Vector6d& sigma, double dt);
    void assemble_J(double dt);

    void init_common();
    void clear();
    bool check_elastic(double dt);

    double dT;
    double Sum_time;
    double elastic_buf;
    void find_SS(double dT);
    void integrate(double dt);
    void apply_rotation(double dt);
    void update_H_KM();

    void computePlasticSpin();
    void updateOrientation(const Eigen::Matrix3d& Omega, double dt);
    double computeDeviatoricDu() const;

public:
    Grain(std::string info);
    Grain(MaterialParamsPtr p, const Eigen::Quaterniond& orientation);
    void Simulate_Step(double dt, const Eigen::Matrix3d& Lnew);

    const Vector6d& get_stress() const;

    Eigen::Quaterniond get_orientation() const;

    void set_orientation(const Eigen::Quaterniond& q);

    Eigen::Matrix3d get_stress_lab() const;

    double get_stress_intensity() const;  // von Mises equivalent stress

    // State variables
    double get_acc_def() const { return acc_def; }
    double get_time() const { return sim_time; }

    // Deformation
    const Eigen::Matrix3d& get_strain() const { return strain_gl; }
    void get_strain_v_into(Vector6d& out) const {
        tensorToVoigtInPlace(strain_gl, out);
    }

    // Slip system state (indexed)
    size_t num_slip_systems() const;  // total (doubled: В±b)
    double get_tc(size_t i) const;
    const std::vector<double>& get_tc_all() const { return tc; }
    double get_shifts(size_t i) const;
    const std::vector<double>& get_shifts_all() const { return shifts; }
    double get_sp_shifts(size_t i) const;
    const std::vector<double>& get_sp_shifts_all() const { return sp_shifts; }
    uint8_t get_last_step(size_t i) const;
    uint8_t get_new_step(size_t i) const;

    long int get_random_state() const { return random_state; }
    void set_random_state(long int seed) { random_state = seed; }
    void reset_state();

    // Hardening matrix
    const Eigen::MatrixXd& get_H() const { return H; }

    double get_volume() const { return volume; }
    void set_volume(double v) { volume = v; }

    const Vector6d& get_Dp_v() const { return Dp_v; }
    void get_elasticity_lab_into(Matrix6d& result) const;
    void get_tangent_crystal_into(Matrix6d& result) const;  // J^-1·P сошедшегося шага (КСК)
    void get_tangent_lab_into(Matrix6d& result) const;      // то же в ЛСК

private:
    mutable Matrix6d T_voigt_buf;
    mutable bool T_voigt_valid = false;
    mutable Eigen::Matrix3d exp_buf;
    void computeVoigtRotationInPlace(Matrix6d& T) const;
};
