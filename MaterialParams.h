#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <memory>

struct MaterialParams {
    using Vector6d = Eigen::Matrix<double, 6, 1>;
    using Matrix6d = Eigen::Matrix<double, 6, 6>;

    std::vector<Eigen::Matrix3d> M;
    std::vector<Eigen::Matrix3d> Wpm;
    std::vector<Vector6d> Ms_v;
    int num_slip_sys;

    Matrix6d P_v;

    // Режимы модели (считываются из файла материала)
    unsigned hard_mode;
    unsigned rot_mode;

    double tc0;
    double strain_rate_infl;
    double n_power = 100.0;

    Eigen::MatrixXd H;
    double A = 0, q = 0, h0 = 0, ts = 1e30, alfa = 0;

    explicit MaterialParams(const std::string& filename);
};

using MaterialParamsPtr = std::shared_ptr<const MaterialParams>;
