#include "Grain.h"
#include "tensor_utils.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <unsupported/Eigen/MatrixFunctions>

Grain::Grain(std::string info)
    : params(std::make_shared<MaterialParams>(info))
{
    H = params->H;

    init_common();
}

Grain::Grain(MaterialParamsPtr p, const Eigen::Quaterniond& orientation)
    : params(p)
{
    H = params->H;

    init_common();

    orient_q = orientation.normalized();
    R_buf = orient_q.toRotationMatrix();
    T_voigt_valid = false;
}

void Grain::init_common()
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

    T_voigt_valid = false;
    random_state = 0;
}

void Grain::clear()
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

bool Grain::check_elastic(double dt)
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

void Grain::find_SS(double dT)
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

void Grain::integrate(double dt)
{
    int n = 2 * params->num_slip_sys;
    
    for (int j = 0; j < n; j++)
    {
        if (!new_step[j]) continue;
        double shift_j = sp_shifts[j];
        shifts[j] += shift_j * dt;
        for (int i = 0; i < n; i++)
        {
            if (tc[i] < params->ts)
                tc[i] += H(i, j) * shift_j * dt;
        }
    }

    if (params->hard_mode == 2)
    {
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

    if (params->rot_mode != 0)
    {
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
}

void Grain::computePlasticSpin()
{
    Wp_mat_buf.setZero();
    for (int i = 0; i < 2 * params->num_slip_sys; i++)
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
    if (angle > ROTATION_ANGLE_THRESHOLD)
    {
        dq_buf = Eigen::AngleAxisd(angle, omega_buf / omega_buf.norm());
        orient_q = orient_q * dq_buf;
        orient_q.normalize();
        R_buf = orient_q.toRotationMatrix();
        T_voigt_valid = false;  // инвалидация кэша при повороте
    }
}

void Grain::Simulate_Step(double dt, const Eigen::Matrix3d& Lnew)
{
    // Compute equivalent strain increment from Lnew (lab frame, rotation-invariant)
    Buf.noalias() = 0.5 * (Lnew + Lnew.transpose());
    tensorToVoigtInPlace(Buf, D_v);
    double du = computeDeviatoricDu();
    double ds = du * dt;

    // Adaptive sub-stepping: if strain increment exceeds ds_crit, split the step
    int n_substeps = 1;
    if (ds > DS_CRITICAL) {
        n_substeps = static_cast<int>(std::ceil(ds / DS_CRITICAL));
    }
    double sub_dt = dt / n_substeps;

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

double Grain::computeDeviatoricDu() const
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

void Grain::reset_state()
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

double Grain::get_stress_intensity() const
{
    double s_mean = (Stress_v(0) + Stress_v(1) + Stress_v(2)) / 3.0;
    double s1 = Stress_v(0) - s_mean;
    double s2 = Stress_v(1) - s_mean;
    double s3 = Stress_v(2) - s_mean;
    double j2s = s1 * s1 + s2 * s2 + s3 * s3
        + 0.5 * (Stress_v(3) * Stress_v(3) + Stress_v(4) * Stress_v(4) + Stress_v(5) * Stress_v(5));
    return std::sqrt(1.5 * j2s);
}

size_t Grain::num_slip_systems() const {
    return params->num_slip_sys * 2;
}

double Grain::get_tc(size_t i) const {
    return tc[i];
}

double Grain::get_shifts(size_t i) const {
    return shifts[i];
}

double Grain::get_sp_shifts(size_t i) const {
    return sp_shifts[i];
}

uint8_t Grain::get_last_step(size_t i) const {
    return last_step[i];
}

uint8_t Grain::get_new_step(size_t i) const {
    return new_step[i];
}

const Vector6d& Grain::get_stress() const
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
    T_voigt_valid = false;  // инвалидация кэша
}

Eigen::Matrix3d Grain::get_stress_lab() const
{
    voigtToTensorInPlace(Stress_v, S_buf);
    Buf.noalias() = R_buf * S_buf;
    S_buf.noalias() = Buf * R_buf.transpose();
    return S_buf;
}

void Grain::computeVoigtRotationInPlace(Matrix6d& T) const {
    T.setZero();
    for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q) {
            int P = tensorToVoigtIndex(p, q);
            bool P_shear = (p != q);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    int I = tensorToVoigtIndex(i, j);
                    bool I_shear = (i != j);
                    double v = R_buf(i, p) * R_buf(j, q);
                    if (!P_shear && !I_shear)       T(I, P) += v;          // норм→норм
                    else if (P_shear && I_shear)    T(I, P) += v;          // сдв→сдв (суммируется дважды)
                    else if (!P_shear && I_shear)   T(I, P) += v * 0.5;          // норм→сдвиг: 1.0 (v(сдвиг)=2·σ, σ'_pp=R_pi·R_pj·σ_ij)
                    else                            T(I, P) += v * 2.0;    // сдвиг→норм: 2.0 (2·σ'_pq=2·R_pi·R_qi·σ_ii)
                }
        }
}

void Grain::get_elasticity_lab_into(Matrix6d& result) const {
    if (!T_voigt_valid) {
        computeVoigtRotationInPlace(T_voigt_buf);
        T_voigt_valid = true;
    }
    Buf_P.noalias() = T_voigt_buf * params->P_v;
    result.noalias() = Buf_P * T_voigt_buf.transpose();
}
