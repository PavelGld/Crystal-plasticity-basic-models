#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <memory>
#include "tensor_utils.h"

struct MaterialParams {

    std::vector<Eigen::Matrix3d> M;
    std::vector<Eigen::Matrix3d> Wpm;
    std::vector<Vector6d> Ms_v;
    std::vector<Vector6d> MsE_v;   // симметричная часть в инженерной нотации (сдвиги ×2)
    int num_slip_sys;

    Matrix6d P_v;

    // Режимы модели (считываются из файла материала)
    unsigned hard_mode;
    unsigned rot_mode;
    unsigned integ_mode = 0;       // 0=явная схема, 1=неявная (Newton-Raphson)

    double tc0;
    double strain_rate_infl;
    double n_power = 100.0;

    Eigen::MatrixXd H;
    double A = 0, q = 0, h0 = 0, ts = 1e30, alfa = 0;
    double iso_a = 0, iso_b = 0;

    explicit MaterialParams(const std::string& filename);
};

using MaterialParamsPtr = std::shared_ptr<const MaterialParams>;
