#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <string>
#include <memory>
#include "Grain.h"

class RVE {
    using Vector6d = Eigen::Matrix<double, 6, 1>;
    using Matrix6d = Eigen::Matrix<double, 6, 6>;

    std::vector<MaterialParamsPtr> materials;
    
    std::vector<Grain> grains;
    
    unsigned rot_mode;

    Eigen::Matrix3d stress_av;
    Eigen::Matrix3d W_av;
    Eigen::Matrix3d Dp_av;
    Eigen::Matrix3d F_av;
    
    double time;
    double acc_strain;
    double total_volume;
    Eigen::Matrix3d A_affine;
    
    // Параметры итеративного метода
    int max_iter = 20;
    double tol = 1e-4;
    
    // Сохранённые значения с прошлого шага
    Vector6d stress_av_old_v = Vector6d::Zero();
    Vector6d Dp_av_old_v = Vector6d::Zero();
    Eigen::Matrix3d W_av_old = Eigen::Matrix3d::Zero();

public:
    explicit RVE(const std::string& config_file);
    
    void Simulate_Step(double dt, const Eigen::Matrix3d& L);
    void Simulate_Step_Soft(double dt, const Eigen::Matrix3d& dS_dt);
    
    Eigen::Matrix3d get_stress() const { return stress_av; }
    Eigen::Matrix3d get_spin() const { return W_av; }
    Eigen::Matrix3d get_plastic_strain_rate() const { return Dp_av; }
    Eigen::Matrix3d get_deformation_gradient() const { return F_av; }
    double get_time() const { return time; }
    double get_accumulated_strain() const { return acc_strain; }
    
    size_t num_grains() const { return grains.size(); }
    const Grain& get_grain(size_t i) const { return grains[i]; }
};
