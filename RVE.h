#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <string>
#include <memory>
#include "tensor_utils.h"
#include "Grain.h"

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
