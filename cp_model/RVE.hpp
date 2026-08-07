#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <unsupported/Eigen/MatrixFunctions>
#include "cp_model/tensor_utils.hpp"
#include "cp_model/Grain.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

class RVE {

    std::vector<MaterialParamsPtr> materials;
    
    std::vector<Grain> grains;

    Eigen::Matrix3d stress_av;
    Eigen::Matrix3d W_av;
    Eigen::Matrix3d Dp_av;
    Eigen::Matrix3d F_av;

    double sim_time;
    double acc_strain;
    double total_volume;
    std::string config_path;
    unsigned master_seed = 0;
    std::vector<int> grains_per_mat;  // saved for reset() without re-reading config

    // Сохранённые значения с прошлого шага
    Vector6d stress_av_old_v = Vector6d::Zero();
    Vector6d Dp_av_old_v = Vector6d::Zero();
    Eigen::Matrix3d W_av_old = Eigen::Matrix3d::Zero();

    // Буферы для избежания аллокаций
    mutable Matrix6d P_av_buf;
    mutable Matrix6d temp_P_buf;
    mutable Vector6d dS_dt_v_buf;
    mutable Vector6d D_v_buf;
    mutable Eigen::Matrix3d D_buf;
    mutable Eigen::Matrix3d L_buf;
    mutable Vector6d stress_v_buf;
    mutable Vector6d Dp_v_buf;
    mutable Eigen::Matrix3d W_sum_buf;
    mutable Eigen::Matrix3d exp_buf;
    mutable Eigen::Matrix3d temp_3x3_buf;
    mutable Eigen::PartialPivLU<Matrix6d> lu_buf;

public:
    explicit RVE(const std::string& config_file, unsigned master_seed = 0);
    
    void reset();
    void Simulate_Step(double dt, const Eigen::Matrix3d& L);
    void Simulate_Step_Soft(double dt, const Eigen::Matrix3d& dS_dt);
    
    Eigen::Matrix3d get_stress() const { return stress_av; }
    Eigen::Matrix3d get_spin() const { return W_av; }
    Eigen::Matrix3d get_plastic_strain_rate() const { return Dp_av; }
    Eigen::Matrix3d get_deformation_gradient() const { return F_av; }
    double get_time() const { return sim_time; }
    double get_accumulated_strain() const { return acc_strain; }
    double get_stress_intensity() const;  // von Mises equivalent stress
    void get_tangent_lab_into(Matrix6d& result) const;  // осреднённый тангенс RVE (Taylor, ЛСК)

    // Previous step values (for soft loading analysis)
    Eigen::Matrix<double, 6, 1> get_stress_voigt_old() const { return stress_av_old_v; }
    Eigen::Matrix<double, 6, 1> get_Dp_voigt_old() const { return Dp_av_old_v; }
    Eigen::Matrix3d get_spin_old() const { return W_av_old; }

    size_t num_grains() const { return grains.size(); }
    const Grain& get_grain(size_t i) const { return grains[i]; }
    Grain& get_grain(size_t i) { return grains[i]; }
    unsigned get_rot_mode() const { return materials.empty() ? 0 : materials[0]->rot_mode; }
    double get_total_volume() const { return total_volume; }
    unsigned get_master_seed() const { return master_seed; }
    void set_master_seed(unsigned seed) { master_seed = seed; }
};

inline void read_rve_config(const std::string& config_file, int& num_materials,
    std::vector<std::string>& mat_files, std::vector<int>& grains_per_mat)
{
    std::ifstream file(config_file);
    if (!file.is_open())
        throw std::runtime_error("Cannot open RVE config file: " + config_file);

    file >> num_materials;
    mat_files.resize(num_materials);
    grains_per_mat.resize(num_materials);
    for (int m = 0; m < num_materials; m++)
        file >> mat_files[m] >> grains_per_mat[m];
}

inline void build_grains_from_config(const std::vector<MaterialParamsPtr>& materials,
    std::vector<Grain>& grains, const std::vector<int>& grains_per_mat, unsigned master_seed)
{
    grains.clear();
    for (size_t m = 0; m < materials.size(); m++) {
        unsigned seed = (master_seed != 0) ? static_cast<unsigned>(std::rand()) : 0;

        for (int g = 0; g < grains_per_mat[m]; g++) {
            Eigen::Quaterniond random_q = Eigen::Quaterniond::UnitRandom();
            grains.emplace_back(materials[m], random_q);
            if (master_seed != 0) {
                unsigned grain_seed = static_cast<unsigned>(std::rand());
                grains.back().set_random_state(grain_seed);
                std::srand(grain_seed);
            }
        }

        if (master_seed != 0)
            std::srand(seed);
    }
}

inline RVE::RVE(const std::string& config_file, unsigned ms)
    : config_path(config_file), master_seed(ms)
{
    if (master_seed != 0)
        std::srand(master_seed);

    int num_materials;
    std::vector<std::string> mat_files;
    read_rve_config(config_file, num_materials, mat_files, grains_per_mat);

    materials.resize(num_materials);
    for (int m = 0; m < num_materials; m++)
        materials[m] = std::make_shared<MaterialParams>(mat_files[m]);

    build_grains_from_config(materials, grains, grains_per_mat, master_seed);

    total_volume = 0;
    for (const auto& grain : grains)
        total_volume += grain.get_volume();

    stress_av.setZero();
    W_av.setZero();
    Dp_av.setZero();
    F_av.setIdentity();
    sim_time = 0;
    acc_strain = 0;

    stress_av_old_v.setZero();
    Dp_av_old_v.setZero();
    W_av_old.setZero();
}

inline void RVE::reset()
{
    if (master_seed != 0)
        std::srand(master_seed);

    build_grains_from_config(materials, grains, grains_per_mat, master_seed);

    total_volume = 0;
    for (const auto& grain : grains)
        total_volume += grain.get_volume();

    stress_av.setZero();
    W_av.setZero();
    Dp_av.setZero();
    F_av.setIdentity();
    sim_time = 0;
    acc_strain = 0;

    stress_av_old_v.setZero();
    Dp_av_old_v.setZero();
    W_av_old.setZero();
}

inline void RVE::Simulate_Step(double dt, const Eigen::Matrix3d& L) {
    stress_av.setZero();
    W_av.setZero();
    Dp_av.setZero();
    Dp_v_buf.setZero();
    
    int n = static_cast<int>(grains.size());
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < n; i++) {
        grains[i].Simulate_Step(dt, L);
    }
    
    for (int i = 0; i < n; i++) {
        double weight = grains[i].get_volume() / total_volume;
        temp_3x3_buf.noalias() = grains[i].get_stress_lab();
        stress_av.noalias() += temp_3x3_buf * weight;

        // Осредняем Dp для диагностики
        Dp_v_buf += weight * grains[i].get_Dp_v();
    }

    Dp_av(0,0) = Dp_v_buf(0);
    Dp_av(1,1) = Dp_v_buf(1);
    Dp_av(2,2) = Dp_v_buf(2);
    Dp_av(1,2) = Dp_av(2,1) = Dp_v_buf(3) / 2;
    Dp_av(0,2) = Dp_av(2,0) = Dp_v_buf(4) / 2;
    Dp_av(0,1) = Dp_av(1,0) = Dp_v_buf(5) / 2;

    // acc_strain из полной деформации (для контроля цикла)
    D_buf = 0.5 * (L + L.transpose());
    double tr_D = D_buf.trace();
    D_buf.diagonal().array() -= tr_D / 3.0;
    double D_eq = std::sqrt(2.0 / 3.0 * D_buf.squaredNorm());
    acc_strain += D_eq * dt;

    sim_time += dt;
    
    exp_buf.noalias() = dt * L;
    if (exp_buf.norm() < 1e-2) {
        applyExpTaylor(exp_buf, F_av, temp_3x3_buf);
    } else {
        F_av.noalias() = exp_buf.exp() * F_av;
    }
}

inline void RVE::Simulate_Step_Soft(double dt, const Eigen::Matrix3d& dS_dt) {
    int n = static_cast<int>(grains.size());
    
    // 1. Вычислить среднюю упругость (в ЛСК, Фойгт)
    P_av_buf.setZero();
    
    for (int i = 0; i < n; i++) {
        double w = grains[i].get_volume() / total_volume;
        grains[i].get_elasticity_lab_into(temp_P_buf);
        P_av_buf.noalias() += w * temp_P_buf;
    }
    
    // 2. Вычислить D из формулы: D = P^-1 : dS/dt + Dp
    tensorToVoigtInPlace(dS_dt, dS_dt_v_buf);
    
    // LU-инверсия без аллокаций
    lu_buf.compute(P_av_buf);
    D_v_buf.noalias() = lu_buf.solve(dS_dt_v_buf);
    D_v_buf.noalias() += Dp_av_old_v;
    
    // 3. Преобразовать в L (спин с прошлого шага)
    voigtToTensorInPlace(D_v_buf, D_buf);
    
    L_buf.noalias() = D_buf + W_av_old;
    
    // 4. Параллельный расчёт зёрен
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < n; i++) {
        grains[i].Simulate_Step(dt, L_buf);
    }
    
    // 5. Осреднить результаты
    stress_v_buf.setZero();
    Dp_v_buf.setZero();

    for (int i = 0; i < n; i++) {
        double w = grains[i].get_volume() / total_volume;

        const Eigen::Matrix3d& s = grains[i].get_stress_lab();
        stress_v_buf(0) += w * s(0,0);
        stress_v_buf(1) += w * s(1,1);
        stress_v_buf(2) += w * s(2,2);
        stress_v_buf(3) += w * (s(1,2) + s(2,1));
        stress_v_buf(4) += w * (s(0,2) + s(2,0));
        stress_v_buf(5) += w * (s(0,1) + s(1,0));

        Dp_v_buf.noalias() += w * grains[i].get_Dp_v();
    }

    // W зависит только от L_buf (одинаков для всех зёрен), Σw = 1
    W_sum_buf.setZero();
    W_sum_buf(0,1) = 0.5 * (L_buf(0,1) - L_buf(1,0));
    W_sum_buf(0,2) = 0.5 * (L_buf(0,2) - L_buf(2,0));
    W_sum_buf(1,2) = 0.5 * (L_buf(1,2) - L_buf(2,1));
    W_sum_buf(1,0) = -W_sum_buf(0,1);
    W_sum_buf(2,0) = -W_sum_buf(0,2);
    W_sum_buf(2,1) = -W_sum_buf(1,2);

    // 6. Сохранить результаты
    stress_av(0,0) = stress_v_buf(0);
    stress_av(1,1) = stress_v_buf(1);
    stress_av(2,2) = stress_v_buf(2);
    stress_av(1,2) = stress_av(2,1) = stress_v_buf(3) / 2;
    stress_av(0,2) = stress_av(2,0) = stress_v_buf(4) / 2;
    stress_av(0,1) = stress_av(1,0) = stress_v_buf(5) / 2;
    
    stress_av_old_v = stress_v_buf;
    Dp_av_old_v = Dp_v_buf;
    W_av_old = W_sum_buf;
    W_av = W_sum_buf;
    
    // Dp_av
    Dp_av(0,0) = Dp_v_buf(0);
    Dp_av(1,1) = Dp_v_buf(1);
    Dp_av(2,2) = Dp_v_buf(2);
    Dp_av(1,2) = Dp_av(2,1) = Dp_v_buf(3) / 2;
    Dp_av(0,2) = Dp_av(2,0) = Dp_v_buf(4) / 2;
    Dp_av(0,1) = Dp_av(1,0) = Dp_v_buf(5) / 2;
    
    sim_time += dt;
    
    exp_buf.noalias() = dt * L_buf;
    if (exp_buf.norm() < 1e-2) {
        applyExpTaylor(exp_buf, F_av, temp_3x3_buf);
    } else {
        F_av.noalias() = exp_buf.exp() * F_av;
    }
}

inline double RVE::get_stress_intensity() const {
    double s_mean = (stress_av(0,0) + stress_av(1,1) + stress_av(2,2)) / 3.0;
    double s00 = stress_av(0,0) - s_mean;
    double s11 = stress_av(1,1) - s_mean;
    double s22 = stress_av(2,2) - s_mean;
    double s12 = stress_av(0,1);
    double s13 = stress_av(0,2);
    double s23 = stress_av(1,2);
    double j2s = s00 * s00 + s11 * s11 + s22 * s22
        + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
    return std::sqrt(1.5 * j2s);
}

inline void RVE::get_tangent_lab_into(Matrix6d& result) const
{
    result.setZero();
    int n = static_cast<int>(grains.size());
    for (int i = 0; i < n; i++) {
        double w = grains[i].get_volume() / total_volume;
        grains[i].get_tangent_lab_into(temp_P_buf);
        result.noalias() += w * temp_P_buf;
    }
}
