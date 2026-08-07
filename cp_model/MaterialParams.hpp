#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include "cp_model/tensor_utils.hpp"

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

    MaterialParams() = default;
    explicit MaterialParams(const std::string& filename);
};

using MaterialParamsPtr = std::shared_ptr<const MaterialParams>;

inline MaterialParams::MaterialParams(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open material file: " + filename);
    
    unsigned hard_mode_temp, rot_mode_temp;
    file >> hard_mode_temp >> rot_mode_temp;
    hard_mode = hard_mode_temp;
    rot_mode = rot_mode_temp;
    file >> num_slip_sys;

    std::vector<double> b(3), n(3);
    M.resize(num_slip_sys * 2);
    Wpm.resize(num_slip_sys * 2);
    Ms_v.resize(num_slip_sys * 2);
    MsE_v.resize(num_slip_sys * 2);

    for (int i = 0; i < num_slip_sys; i++) {
        for (int j = 0; j < 3; j++) file >> b[j];
        for (int j = 0; j < 3; j++) file >> n[j];

        if (std::abs(n[0]*b[0] + n[1]*b[1] + n[2]*b[2]) > 1e-10)
            std::cout << "Warning: b and n not orthogonal for system " << i << std::endl;

        double norm_n = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        double norm_b = std::sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
        for (int j = 0; j < 3; j++) { n[j] /= norm_n; b[j] /= norm_b; }

        M[i*2].setZero();
        Wpm[i*2].setZero();
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) {
                M[i*2](j, k) = b[k] * n[j];
                Wpm[i*2](j, k) = 0.5 * (b[j]*n[k] - b[k]*n[j]);
            }
        Eigen::Matrix3d Ms = 0.5 * (M[i*2] + M[i*2].transpose());
        Ms_v[i*2](0) = Ms(0, 0);
        Ms_v[i*2](1) = Ms(1, 1);
        Ms_v[i*2](2) = Ms(2, 2);
        Ms_v[i*2](3) = Ms(1, 2);  // тензорная компонента (без ×2)
        Ms_v[i*2](4) = Ms(0, 2);  // чтобы σ_v · Ms_v = σ : Ms
        Ms_v[i*2](5) = Ms(0, 1);
        MsE_v[i*2] = Ms_v[i*2];
        MsE_v[i*2](3) *= 2.0; MsE_v[i*2](4) *= 2.0; MsE_v[i*2](5) *= 2.0;

        for (int j = 0; j < 3; j++) b[j] = -b[j];
        M[i*2+1].setZero();
        Wpm[i*2+1].setZero();
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) {
                M[i*2+1](j, k) = b[k] * n[j];
                Wpm[i*2+1](j, k) = 0.5 * (b[j]*n[k] - b[k]*n[j]);
            }
        Ms = 0.5 * (M[i*2+1] + M[i*2+1].transpose());
        Ms_v[i*2+1](0) = Ms(0, 0);
        Ms_v[i*2+1](1) = Ms(1, 1);
        Ms_v[i*2+1](2) = Ms(2, 2);
        Ms_v[i*2+1](3) = Ms(1, 2);
        Ms_v[i*2+1](4) = Ms(0, 2);
        Ms_v[i*2+1](5) = Ms(0, 1);
        MsE_v[i*2+1] = Ms_v[i*2+1];
        MsE_v[i*2+1](3) *= 2.0; MsE_v[i*2+1](4) *= 2.0; MsE_v[i*2+1](5) *= 2.0;
    }

    double c1, c2, c3;
    file >> c1 >> c2 >> c3;

    P_v.setZero();
    P_v(0,0) = P_v(1,1) = P_v(2,2) = c1;
    P_v(0,1) = P_v(1,0) = P_v(1,2) = P_v(2,1) = P_v(0,2) = P_v(2,0) = c2;
    P_v(3,3) = P_v(4,4) = P_v(5,5) = c3;

    H.resize(num_slip_sys * 2, num_slip_sys * 2);

    if (hard_mode == 2) {
        file >> A >> q >> h0 >> ts >> alfa;
        H.setConstant(A * q);
        H.diagonal().setConstant(A);
    }
    else if (hard_mode == 1) {
        file >> iso_a >> iso_b;
    }

    file >> tc0;
    file >> strain_rate_infl;
    file >> n_power;

    unsigned im;
    if (file >> im) { integ_mode = im; if (integ_mode > 1) throw std::runtime_error("integ_mode must be 0 or 1"); }   // опциональный токен; EOF → остаётся 0
}
