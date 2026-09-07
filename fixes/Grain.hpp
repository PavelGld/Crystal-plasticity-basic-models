#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <string>
#include <memory>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <unsupported/Eigen/MatrixFunctions>
#ifndef ABAQUS_UMAT
#include <iostream>
#include <fstream>
#include <stdexcept>
#endif
#include "cp_model/tensor_utils.hpp"
#include "cp_model/MaterialParams.hpp"

class Grain
{

#ifdef ABAQUS_UMAT
    const MaterialParams* params;   // raw pointer (UMAT: no CRT-heavy shared_ptr)
#else
    MaterialParamsPtr params;       // shared_ptr (standalone)
#endif

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
#ifndef CP_DS_CRIT_IMPL
#define CP_DS_CRIT_IMPL 2.5e-4   // порог подшага (экв. деформация за подшаг). 2.5e-4: ошибка BE ~0.1% при Δε=5e-4/инкремент;
#endif                           // 1.25e-4: ~0.02% (+50% времени). Переопределяется /DCP_DS_CRIT_IMPL=...
    static constexpr double DS_CRITICAL_IMPL = CP_DS_CRIT_IMPL;  // порог подшагивания для неявной схемы

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

    void solve_implicit(double dt, int depth, double weight);
    bool step_failed = false;   // local integrator failed at max fallback depth
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

    void computePlasticSpin();
    void updateOrientation(const Eigen::Matrix3d& Omega, double dt);
    double computeDeviatoricDu() const;

public:
    Grain(std::string info);
#ifndef ABAQUS_UMAT
    Grain(MaterialParamsPtr p, const Eigen::Quaterniond& orientation);
#else
    Grain(const MaterialParams* p, const Eigen::Quaterniond& orientation);
#endif
    void Simulate_Step(double dt, const Eigen::Matrix3d& Lnew);

    const Vector6d& get_stress() const;

    Eigen::Quaterniond get_orientation() const;

    void set_orientation(const Eigen::Quaterniond& q);
    // NLGEOM=YES: приращение жёсткого вращения базиса Abaqus (DROT). Состояние хранится в КСК,
    // поэтому поворот базиса эквивалентен повороту ориентации решётки: R_lat <- dR * R_lat.
    void rotate_lattice(const Eigen::Matrix3d& dR) {
        Eigen::Quaterniond dq(dR);
        orient_q = (dq * orient_q).normalized();
        R_buf = orient_q.toRotationMatrix();
        T_voigt_valid = false;
    }

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
    void update_H_KM();
    bool get_step_failed() const { return step_failed; }
    void update_H_KM_if_needed() { if (params->hard_mode == 2) update_H_KM(); }

    // Hardening matrix
    const Eigen::MatrixXd& get_H() const { return H; }

    double get_volume() const { return volume; }
    void set_volume(double v) { volume = v; }

    const Vector6d& get_Dp_v() const { return Dp_v; }
    void get_elasticity_lab_into(Matrix6d& result) const;
    void get_tangent_crystal_into(Matrix6d& result) const;  // J^-1·P сошедшегося шага (КСК)
    void get_tangent_lab_into(Matrix6d& result) const;      // то же в ЛСК

    // State persistence for UMAT
    void set_stress_v(const Vector6d& sv) { Stress_v = sv; }
    void set_tc_all(const std::vector<double>& t) { tc = t; }
    void set_shifts_all(const std::vector<double>& s) { shifts = s; }
    void set_acc_def(double a) { acc_def = a; }
    void set_sim_time(double t) { sim_time = t; }

private:
    mutable Matrix6d T_voigt_buf;
    mutable Matrix6d T_voigt_inv_buf;
    mutable bool T_voigt_valid = false;
    mutable Eigen::Matrix3d exp_buf;
    void computeVoigtRotationInPlace(Matrix6d& T, const Eigen::Matrix3d& R) const;
};

#ifndef ABAQUS_UMAT
inline Grain::Grain(std::string info)
    : params(std::make_shared<MaterialParams>(info))
{
    H = params->H;

    init_common();
}
#endif

#ifdef ABAQUS_UMAT
inline Grain::Grain(const MaterialParams* p, const Eigen::Quaterniond& orientation)
    : params(p)
{
    H = params->H;
    init_common();
    orient_q = orientation.normalized();
    R_buf = orient_q.toRotationMatrix();
    T_voigt_valid = false;
}
#else
inline Grain::Grain(MaterialParamsPtr p, const Eigen::Quaterniond& orientation)
    : params(p)
{
    H = params->H;

    init_common();

    orient_q = orientation.normalized();
    R_buf = orient_q.toRotationMatrix();
    T_voigt_valid = false;
}
#endif

inline void Grain::init_common()
{
    sim_time = 0;
    acc_def = 0;

    Buf.setZero();
    Buf_v.setZero();
    Buf_PD.setZero();
    Buf_P.setZero();
    Wp.setZero();
    L.setZero();
    strain_gl.setZero();
    F.setIdentity();

    Stress.setZero();
    D_v.setZero();
    Dp_v.setZero();
    Stress_v.setZero();

    orient_q = Eigen::Quaterniond::Identity();
    R_buf.setIdentity();
    omega_buf.setZero();
    dq_buf = Eigen::Quaterniond::Identity();
    W_mat_buf.setZero();
    Omega_buf.setZero();
    Wp_mat_buf.setZero();
    S_buf.setZero();

    int n = 2 * params->num_slip_sys;
    last_step.assign(n, 0);
    new_step.assign(n, 0);
    tc.assign(n, params->tc0);
    shifts.assign(n, 0.0);
    sp_shifts.assign(n, 0.0);

    gamma_dot.setZero(n);
    tau_buf.setZero(n);
    tc_eff.setZero(n);
    G_buf.setZero(n, 6);
    G_prev.setZero(n, 6);
    J_buf.setZero();
    dyad_buf.setZero();
    res_v.setZero();
    res_new_v.setZero();
    sigma_trial_v.setZero();
    sigma_new_v.setZero();
    delta_sigma_v.setZero();
    trial_inc_v.setZero();

    tangent_buf.setZero();
    tangent_valid = false;

    T_voigt_valid = false;
    random_state = 0;
}

inline void Grain::clear()
{
    Buf.setZero();
    Buf_v.setZero();
    Wp.setZero();
    Dp_v.setZero();
    D_v.setZero();
    for (int i = 0; i < 2 * params->num_slip_sys; i++)
    {
        last_step[i] = new_step[i];
        new_step[i] = 0;
    }
}

inline bool Grain::check_elastic(double dt)
{
    int n = 2 * params->num_slip_sys;
    
    Buf_PD.noalias() = params->P_v * D_v;
    
    Buf_v.noalias() = Stress_v + Buf_PD * dt;

    bool TF = true;
    for (int i = 0; i < n; i++)
    {
        double t = Buf_v.dot(params->Ms_v[i]);
        if (t > tc[i])
            TF = false;
    }

    if (!TF)
    {
        bool was_in_plastic = false;
        for (int i = 0; i < n; i++)
        {
            double t = Stress_v.dot(params->Ms_v[i]);
            if (t > tc[i])
                was_in_plastic = true;
        }

        if (was_in_plastic)
        {
            elastic_buf = 0;  // зерно уже в пластике — весь шаг пластический
            Buf_v.setZero();
            return TF;
        }
        else
        {
            double buft = dt;
            for (int i = 0; i < n; i++)
            {
                double ai = Buf_PD.dot(params->Ms_v[i]);
                double t = Stress_v.dot(params->Ms_v[i]);
                if (ai > 0)
                    buft = std::min(buft, (tc[i] - t) / ai);
            }
            // Упругая фаза: интегрируем с полным шагом buft
            Stress_v.noalias() += Buf_PD * buft;
            Buf_v.setZero();
            elastic_buf = buft;
            integrate(buft);
            return false;
        }
    }

    for (int i = 0; i < n; i++)
    {
        sp_shifts[i] = 0.0;
        new_step[i] = 0;
    }
    Stress_v = Buf_v;
    Buf_v.setZero();
    integrate(dt);
    return true;
}

inline void Grain::find_SS(double dT)
{ 
    int n = 2 * params->num_slip_sys;
    
    double du = computeDeviatoricDu();

    for (int i = 0; i < n; i++)
        if (last_step[i]) {
            // Ms_v хранит тензорные компоненты сдвига (без ×2),
            // но D_v — в инженерной нотации Фойгта (с ×2 на сдвиге).
            // Конвертируем Dp_v в инженерную нотацию для корректного
            // вычитания (D_v - Dp_v) перед умножением на P_v.
            // Без этого сдвиговые компоненты Dp_v занижались в 2 раза,
            // что приводило к завышению напряжений на ~20-25%.
            Dp_v(0) += params->Ms_v[i](0) * sp_shifts[i];
            Dp_v(1) += params->Ms_v[i](1) * sp_shifts[i];
            Dp_v(2) += params->Ms_v[i](2) * sp_shifts[i];
            Dp_v(3) += 2.0 * params->Ms_v[i](3) * sp_shifts[i];
            Dp_v(4) += 2.0 * params->Ms_v[i](4) * sp_shifts[i];
            Dp_v(5) += 2.0 * params->Ms_v[i](5) * sp_shifts[i];
        }

    Buf_v.noalias() = params->P_v * (D_v - Dp_v) * dT;
    Buf_v += Stress_v;
    
    double rate_factor = params->strain_rate_infl * du;

    for (int i = 0; i < n; i++)
    {
        double t = Buf_v.dot(params->Ms_v[i]);
        if (t > tc[i])
        {
            new_step[i] = 1;
            sp_shifts[i] = rate_factor * std::pow(std::abs(t / tc[i]), params->n_power);
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

inline void Grain::integrate(double dt)
{
    int n = 2 * params->num_slip_sys;
    
    for (int j = 0; j < n; j++)
    {
        if (!new_step[j]) continue;
        double shift_j = sp_shifts[j];
        shifts[j] += shift_j * dt;

        if (params->hard_mode != 1) {
            for (int i = 0; i < n; i++)
            {
                if (tc[i] < params->ts)
                    tc[i] += H(i, j) * shift_j * dt;
            }
        }
    }

    if (params->hard_mode == 2)
        update_H_KM();
    else if (params->hard_mode == 1)
    {
        double total_gamma = 0.0;
        for (int i = 0; i < n; i++)
            total_gamma += shifts[i];
        double new_tc = params->tc0 + params->iso_a * std::pow(total_gamma, params->iso_b);
        for (int i = 0; i < n; i++)
            tc[i] = new_tc;
    }

    apply_rotation(dt);
}

inline void Grain::update_H_KM()
{
    int n = 2 * params->num_slip_sys;
    if (H_buf.size() != n) H_buf.resize(n);
    for (int j = 0; j < n; j++)
    {
        double t_j = tc[j];
        double Hii = 0;
        if (t_j <= params->ts)
            Hii = params->A * params->h0 * std::pow(std::abs(1.0 - t_j / params->ts), params->alfa);
        H_buf(j) = Hii;
    }
    H.noalias() = params->q * Eigen::VectorXd::Ones(n) * H_buf.transpose();
    H.diagonal() = H_buf;
}

inline void Grain::apply_rotation(double dt)
{
    if (params->rot_mode == 0) return;

    W_mat_buf = 0.5 * (L - L.transpose());

    if (params->rot_mode == 1)
    {
        Omega_buf = W_mat_buf;
    }
    else if (params->rot_mode == 2)
    {
        computePlasticSpin();
        Omega_buf.noalias() = W_mat_buf - Wp_mat_buf;
    }

    updateOrientation(Omega_buf, dt);
}

// Оценка Dp и производных для заданного σ (детерминированная функция σ:
// фиксированное число sweep'ов fixed-point по γ̇↔τc).
// Заполняет: tau_buf, gamma_dot, Dp_v (инженерная), tc_eff, G_prev (∂γ̇/∂σ).
inline void Grain::evaluate_Dp(const Vector6d& sigma, double dt)
{
    int n = 2 * params->num_slip_sys;

    for (int i = 0; i < n; i++)
        tau_buf(i) = sigma.dot(params->Ms_v[i]);

    gamma_dot.setZero();

    int sweeps = (params->hard_mode == 0) ? 1 : 3;
    for (int s = 0; s < sweeps; s++)
    {
        // эффективные τc от текущего итерата γ̇
        if (params->hard_mode == 2) {
            // H = q·1·hᵀ + (1−q)·diag(h)  →  Σ_j H(i,j)γ̇_j = q·Σ_j h_j γ̇_j + (1−q)·h_i γ̇_i  (O(n))
            double s_all = 0.0;
            for (int j = 0; j < n; j++) s_all += H(j, j) * gamma_dot(j);
            for (int i = 0; i < n; i++) {
                double dtc = params->q * s_all + (1.0 - params->q) * H(i, i) * gamma_dot(i);
                tc_eff(i) = tc[i] + dtc * dt;
            }
        } else if (params->hard_mode == 1) {
            double gamma_tot = 0.0;
            for (int i = 0; i < n; i++)
                gamma_tot += shifts[i] + gamma_dot(i) * dt;
            double tc_iso = params->tc0 + params->iso_a * std::pow(gamma_tot, params->iso_b);
            tc_eff.setConstant(tc_iso);
        } else {
            for (int i = 0; i < n; i++)
                tc_eff(i) = tc[i];
        }

        // Чистый степенной закон (без порога): γ̇ = rf·(τ/τc)^n для τ > 0.
        // Решение от 2026-07-26: пороговый закон разрывен на τ=τc → у невязки
        // Ньютона нет корня на упруго-пластическом переходе. Явная схема
        // сохраняет порог; эквивалентность проверяется тестами (5%).
        for (int i = 0; i < n; i++) {
            if (tau_buf(i) > 0.0) {
                double ratio = tau_buf(i) / tc_eff(i);
                if (ratio > RATIO_MAX) ratio = RATIO_MAX;   // защита от переполнения pow
                gamma_dot(i) = rate_factor_impl * std::pow(ratio, params->n_power);
            } else {
                gamma_dot(i) = 0.0;
            }
        }
    }

    // Dp в инженерной нотации Фойгта (сдвиги ×2 — как в find_SS)
    Dp_v.setZero();
    for (int i = 0; i < n; i++)
        if (gamma_dot(i) > 0.0)
            Dp_v.noalias() += gamma_dot(i) * params->MsE_v[i];

    // G = ∂γ̇/∂σ (n×6), та же схема sweep'ов; связка упрочнения — с предыдущего sweep'а
    G_prev.setZero();
    for (int s = 0; s < sweeps; s++)
    {
        // rank-1 часть связки KM: R = Σ_j h_j G_prev.row(j)  (один проход, O(n·6))
        Eigen::Matrix<double,1,6> R_all; R_all.setZero();
        if (params->hard_mode == 2)
            for (int j = 0; j < n; j++)
                if (gamma_dot(j) > 0.0) R_all += H(j, j) * G_prev.row(j);
        for (int i = 0; i < n; i++) {
            if (gamma_dot(i) > 0.0) {
                // производные в точке клампированного ratio (осознанное приближение за клампом)
                double a_i = params->n_power * gamma_dot(i) / tau_buf(i);
                double b_i = -params->n_power * gamma_dot(i) / tc_eff(i);
                G_buf.row(i) = a_i * params->Ms_v[i].transpose();
                if (params->hard_mode == 2) {
                    G_buf.row(i) += (b_i * dt) * (params->q * R_all + (1.0 - params->q) * H(i, i) * G_prev.row(i));
                } else if (params->hard_mode == 1) {
                    double gamma_tot = 0.0;
                    for (int j = 0; j < n; j++)
                        gamma_tot += shifts[j] + gamma_dot(j) * dt;
                    double dh = (gamma_tot > 1e-30)
                        ? params->iso_a * params->iso_b * std::pow(gamma_tot, params->iso_b - 1.0)
                        : 0.0;
                    G_buf.row(i) += (b_i * dh * dt) * G_prev.colwise().sum();
                }
            } else {
                G_buf.row(i).setZero();
            }
        }
        G_prev.swap(G_buf);
    }
}

// J = I + Δt · P · Σ_α MsE^α ⊗ G^α
inline void Grain::assemble_J(double dt)
{
    int n = 2 * params->num_slip_sys;
    dyad_buf.setZero();
    for (int i = 0; i < n; i++)
        if (gamma_dot(i) > 0.0)
            dyad_buf.noalias() += params->MsE_v[i] * G_prev.row(i);
    J_buf.noalias() = Matrix6d::Identity() + params->P_v * dyad_buf * dt;
}

// Ньютон по σ с line search Армихо; fallback — деление шага пополам.
// Требует: D_v (КСК, инженерная) и L (КСК) уже установлены, clear() выполнен.
inline void Grain::solve_implicit(double dt, int depth, double weight)
{
    int n = 2 * params->num_slip_sys;

    rate_factor_impl = params->strain_rate_infl * computeDeviatoricDu();

    trial_inc_v.noalias() = params->P_v * D_v;
    sigma_trial_v.noalias() = Stress_v + trial_inc_v * dt;

    Buf_v = sigma_trial_v;   // текущий итерат σ
    // Стартовая точка Ньютона: если упругий предиктор выводит какую-то систему за
    // τ/τc > 1, масштабируем приращение так, чтобы max τ/τc = 1. При большом n
    // (≈ rate-independent) это убирает старт "на обрыве" степенного закона и
    // сокращает число итераций в 3–5 раз. На результат не влияет (только старт).
    {
        double lam = 1.0;
        for (int i = 0; i < n; i++) {
            double t_old = Stress_v.dot(params->Ms_v[i]);
            double t_tr  = sigma_trial_v.dot(params->Ms_v[i]);
            if (t_tr > tc[i] && t_tr > t_old) {
                double l = (tc[i] - t_old) / (t_tr - t_old);
                if (l < lam) lam = l;
            }
        }
        if (lam < 1.0) {
            if (lam < 0.0) lam = 0.0;
            Buf_v = Stress_v + lam * trial_inc_v * dt;
        }
    }

    evaluate_Dp(Buf_v, dt);
    res_v = Buf_v - sigma_trial_v + params->P_v * Dp_v * dt;
    double scale = std::max(sigma_trial_v.norm(), params->tc0);
    double rnorm = res_v.norm();

    bool converged = (rnorm < NEWTON_TOL * scale);
    for (int iter = 0; iter < NEWTON_MAX_ITER && !converged; iter++)
    {
        assemble_J(dt);
        Eigen::PartialPivLU<Matrix6d> lu(J_buf);
        delta_sigma_v = lu.solve(-res_v);

        double lambda = 1.0;
        double rnorm_new = rnorm;
        bool accepted = false;
        for (int bt = 0; bt <= NEWTON_MAX_BACKTRACKS; bt++)
        {
            sigma_new_v = Buf_v + lambda * delta_sigma_v;
            evaluate_Dp(sigma_new_v, dt);
            res_new_v = sigma_new_v - sigma_trial_v + params->P_v * Dp_v * dt;
            rnorm_new = res_new_v.norm();
            if (rnorm_new <= (1.0 - 1e-4 * lambda) * rnorm) { accepted = true; break; }
            lambda *= 0.5;
        }
        if (!accepted) break;   // → fallback

        Buf_v = sigma_new_v;
        res_v = res_new_v;
        rnorm = rnorm_new;
        converged = (rnorm < NEWTON_TOL * scale);
    }

    if (!converged)
    {
            if (depth >= MAX_FALLBACK_DEPTH) {
#ifndef ABAQUS_UMAT
                throw std::runtime_error("implicit Newton failed to converge at max fallback depth");
#else
                tangent_valid = false;  // signal: return elastic matrix, stale tangent discarded
                step_failed = true;     // -> CPUMAT::step sets PNEWDT < 1
                return;
#endif
            }
        solve_implicit(0.5 * dt, depth + 1, 0.5 * weight);
        solve_implicit(0.5 * dt, depth + 1, 0.5 * weight);
        return;
    }

    // Обновление состояния по сошедшимся γ̇
    // (последний evaluate_Dp был именно в Buf_v → gamma_dot/tc_eff/G_prev соответствуют точке сходимости)
    Stress_v = Buf_v;
    for (int i = 0; i < n; i++) {
        sp_shifts[i] = gamma_dot(i);
        new_step[i] = (gamma_dot(i) > 0.0) ? 1 : 0;
        shifts[i] += gamma_dot(i) * dt;
        tc[i] = tc_eff(i);
    }
    if (params->hard_mode == 2)
        update_H_KM();

    // Консистентный тангенс: J^{-1}·(P - c·g^T) на сошедшемся состоянии
    // c = P·Dᵖ, g = ∂ln(du)/∂D_v. Поправка учитывает зависимость rate_factor от D.
    assemble_J(dt);
    {
        Eigen::PartialPivLU<Matrix6d> lu(J_buf);
        double du = std::max(computeDeviatoricDu(), 1e-30);
        if (du > 0 && Dp_v.squaredNorm() > 0) {
            double inv_3du2 = 1.0 / (3.0 * du * du);
            Vector6d c = params->P_v * Dp_v;
            Vector6d g;
            double tr3 = (D_v(0) + D_v(1) + D_v(2)) / 3.0;
            g(0) = 2.0 * (D_v(0) - tr3) * inv_3du2;
            g(1) = 2.0 * (D_v(1) - tr3) * inv_3du2;
            g(2) = 2.0 * (D_v(2) - tr3) * inv_3du2;
            g(3) = D_v(3) * inv_3du2;
            g(4) = D_v(4) * inv_3du2;
            g(5) = D_v(5) * inv_3du2;
            Matrix6d P_corrected = params->P_v;
            P_corrected.noalias() -= c * g.transpose();
            // Цепной консистентный тангенс по подшагам: C <- J^{-1}(C_prev + w*P_eff)
            tangent_buf = lu.solve((tangent_buf + weight * P_corrected).eval());
        } else {
            tangent_buf = lu.solve((tangent_buf + weight * params->P_v).eval());
        }
    }
    tangent_valid = true;

    apply_rotation(dt);
}

inline void Grain::computePlasticSpin()
{
    Wp_mat_buf.setZero();
    for (int i = 0; i < 2 * params->num_slip_sys; i++)
        if (new_step[i])
            Wp_mat_buf.noalias() += sp_shifts[i] * params->Wpm[i];
}

inline void Grain::updateOrientation(const Eigen::Matrix3d& Omega, double dt)
{
    omega_buf(0) = Omega(2, 1) - Omega(1, 2);
    omega_buf(1) = Omega(0, 2) - Omega(2, 0);
    omega_buf(2) = Omega(1, 0) - Omega(0, 1);
    omega_buf *= 0.5;

    double angle = omega_buf.norm() * dt;
    if (angle > ROTATION_ANGLE_THRESHOLD)
    {
        dq_buf = Eigen::AngleAxisd(angle, omega_buf / omega_buf.norm());
        orient_q = orient_q * dq_buf;
        orient_q.normalize();
        R_buf = orient_q.toRotationMatrix();
        T_voigt_valid = false;  // инвалидация кэша при повороте
    }
}

inline void Grain::Simulate_Step(double dt, const Eigen::Matrix3d& Lnew)
{
    // Compute equivalent strain increment from Lnew (lab frame, rotation-invariant)
    Buf.noalias() = 0.5 * (Lnew + Lnew.transpose());
    tensorToVoigtInPlace(Buf, D_v);
    double du = computeDeviatoricDu();
    double ds = du * dt;

    // Adaptive sub-stepping: if strain increment exceeds ds_crit, split the step
    double ds_crit = (params->integ_mode == 1) ? DS_CRITICAL_IMPL : DS_CRITICAL;
    int n_substeps = 1;
    if (ds > ds_crit) {
        n_substeps = static_cast<int>(std::ceil(ds / ds_crit));
    }
#ifdef ABAQUS_UMAT
    // Пробная деформация, требующая аномально много подшагов, — признак мусорного
    // глобального итерата. Не перемалываем её, а сигнализируем о срезке шага.
    if (params->integ_mode == 1 && n_substeps > 200) {
        step_failed = true;
        tangent_valid = false;
        return;
    }
#endif
    double sub_dt = dt / n_substeps;

    step_failed = false;
    if (params->integ_mode == 1)
        tangent_buf.setZero();   // аккумулятор цепного тангенса на весь инкремент

    double time_done = 0;
    for (int s = 0; s < n_substeps; s++) {
        double actual_dt = (s == n_substeps - 1) ? (dt - time_done) : sub_dt;

        clear();

        // Transform L to crystal frame (orientation may have changed on previous substep)
        if (params->rot_mode != 0) {
            L.noalias() = R_buf.transpose() * Lnew * R_buf;
        } else {
            L = Lnew;
        }
        Buf.noalias() = 0.5 * (L + L.transpose());
        tensorToVoigtInPlace(Buf, D_v);

        // Elastic-plastic split for this substep
        if (params->integ_mode == 1) {
            solve_implicit(actual_dt, 0, actual_dt / dt);
        } else {
            dT = actual_dt;
            Sum_time = 0;
            elastic_buf = 0;

            while (actual_dt > Sum_time + TIME_EPSILON) {
                if (check_elastic(dT))
                    Sum_time += dT;
                else {
                    double dt_plastic = dT - elastic_buf;
                    Sum_time += dT;
                    if (dt_plastic > TIME_EPSILON)
                        find_SS(dt_plastic);
                }
                dT = actual_dt - Sum_time;
            }
        }

        exp_buf.noalias() = actual_dt * L;
        if (exp_buf.norm() < 1e-2) {
            applyExpTaylor(exp_buf, F, Buf);
        } else {
            F.noalias() = exp_buf.exp() * F;
        }

        time_done += actual_dt;
    }

    strain_gl.noalias() = 0.5 * (F.transpose() * F - Eigen::Matrix3d::Identity());
    sim_time += dt;
    acc_def += ds;
}

inline double Grain::computeDeviatoricDu() const
{
    double trD = D_v(0) + D_v(1) + D_v(2);
    double mean = trD / 3.0;
    double d1 = D_v(0) - mean;
    double d2 = D_v(1) - mean;
    double d3 = D_v(2) - mean;
    double j2 = d1 * d1 + d2 * d2 + d3 * d3
        + 0.5 * (D_v(3) * D_v(3) + D_v(4) * D_v(4) + D_v(5) * D_v(5));
    return std::sqrt(2.0 / 3.0 * j2);
}

inline void Grain::reset_state()
{
    sim_time = 0;
    acc_def = 0;
    F.setIdentity();
    strain_gl.setZero();
    Stress.setZero();
    Stress_v.setZero();
    D_v.setZero();
    Dp_v.setZero();
    Wp.setZero();
    L.setZero();
    T_voigt_valid = false;
    tangent_valid = false;

    int n = 2 * params->num_slip_sys;
    for (int i = 0; i < n; i++)
    {
        last_step[i] = 0;
        new_step[i] = 0;
        tc[i] = params->tc0;
        shifts[i] = 0.0;
        sp_shifts[i] = 0.0;
    }

    if (params->hard_mode == 2)
    {
        H = params->H;
    }
}

inline double Grain::get_stress_intensity() const
{
    double s_mean = (Stress_v(0) + Stress_v(1) + Stress_v(2)) / 3.0;
    double s1 = Stress_v(0) - s_mean;
    double s2 = Stress_v(1) - s_mean;
    double s3 = Stress_v(2) - s_mean;
    double j2s = s1 * s1 + s2 * s2 + s3 * s3
        + 0.5 * (Stress_v(3) * Stress_v(3) + Stress_v(4) * Stress_v(4) + Stress_v(5) * Stress_v(5));
    return std::sqrt(1.5 * j2s);
}

inline size_t Grain::num_slip_systems() const {
    return params->num_slip_sys * 2;
}

inline double Grain::get_tc(size_t i) const {
    return tc[i];
}

inline double Grain::get_shifts(size_t i) const {
    return shifts[i];
}

inline double Grain::get_sp_shifts(size_t i) const {
    return sp_shifts[i];
}

inline uint8_t Grain::get_last_step(size_t i) const {
    return last_step[i];
}

inline uint8_t Grain::get_new_step(size_t i) const {
    return new_step[i];
}

inline const Vector6d& Grain::get_stress() const
{
    return Stress_v;
}

inline Eigen::Quaterniond Grain::get_orientation() const
{
    return orient_q;
}

inline void Grain::set_orientation(const Eigen::Quaterniond& q)
{
    orient_q = q.normalized();
    R_buf = orient_q.toRotationMatrix();
    T_voigt_valid = false;  // инвалидация кэша
}

inline Eigen::Matrix3d Grain::get_stress_lab() const
{
    voigtToTensorInPlace(Stress_v, S_buf);
    Buf.noalias() = R_buf * S_buf;
    S_buf.noalias() = Buf * R_buf.transpose();
    return S_buf;
}

inline void Grain::computeVoigtRotationInPlace(Matrix6d& T, const Eigen::Matrix3d& R) const {
    T.setZero();
    for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q) {
            int P = tensorToVoigtIndex(p, q);
            bool P_shear = (p != q);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    int I = tensorToVoigtIndex(i, j);
                    bool I_shear = (i != j);
                    double v = R(i, p) * R(j, q);
                    if (!P_shear && !I_shear)       T(I, P) += v;          // norm->norm
                    else if (P_shear && I_shear)    T(I, P) += v * 0.5;    // shear->shear
                    else if (!P_shear && I_shear)   T(I, P) += v;          // norm->shear
                    else                            T(I, P) += v * 0.5;    // shear->norm
                }
        }
}

inline void Grain::get_elasticity_lab_into(Matrix6d& result) const {
    if (!T_voigt_valid) {
        computeVoigtRotationInPlace(T_voigt_buf, R_buf);
        computeVoigtRotationInPlace(T_voigt_inv_buf, R_buf.transpose());
        T_voigt_valid = true;
    }
    Buf_P.noalias() = T_voigt_buf * params->P_v;
    result.noalias() = Buf_P * T_voigt_inv_buf;
}

inline void Grain::get_tangent_crystal_into(Matrix6d& result) const
{
    if (!tangent_valid) {
#ifndef ABAQUS_UMAT
        throw std::runtime_error("get_tangent_crystal_into: no implicit step performed yet");
#else
        result = params->P_v;
        return;
#endif
    }
    result = tangent_buf;
}

inline void Grain::get_tangent_lab_into(Matrix6d& result) const
{
    if (!tangent_valid) {
#ifndef ABAQUS_UMAT
        throw std::runtime_error("get_tangent_lab_into: no implicit step performed yet");
#else
        get_elasticity_lab_into(result);
        return;
#endif
    }
    if (!T_voigt_valid) {
        computeVoigtRotationInPlace(T_voigt_buf, R_buf);
        computeVoigtRotationInPlace(T_voigt_inv_buf, R_buf.transpose());
        T_voigt_valid = true;
    }
    Buf_P.noalias() = T_voigt_buf * tangent_buf;
    result.noalias() = Buf_P * T_voigt_inv_buf;
}
