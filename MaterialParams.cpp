#include "MaterialParams.h"
#include "tensor_utils.h"
#include <iostream>
#include <fstream>
#include <cmath>

static Eigen::Matrix<double, 6, 1> tensorToVoigt(const Eigen::Matrix3d& A) {
    Eigen::Matrix<double, 6, 1> v;
    v.setZero();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            v(tensorToVoigtIndex(i, j)) += A(i, j);
    return v;
}

MaterialParams::MaterialParams(const std::string& filename) {
    std::ifstream file(filename);
    
    unsigned hard_mode_temp, rot_mode_temp;
    file >> hard_mode_temp >> rot_mode_temp;
    hard_mode = hard_mode_temp;
    rot_mode = rot_mode_temp;
    file >> num_slyde_sys;

    std::vector<double> b(3), n(3);
    M.resize(num_slyde_sys * 2);
    Wpm.resize(num_slyde_sys * 2);
    Ms_v.resize(num_slyde_sys * 2);

    for (int i = 0; i < num_slyde_sys; i++) {
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
    }

    double c1, c2, c3;
    file >> c1 >> c2 >> c3;

    P_v.setZero();
    P_v(0,0) = P_v(1,1) = P_v(2,2) = c1;
    P_v(0,1) = P_v(1,0) = P_v(1,2) = P_v(2,1) = P_v(0,2) = P_v(2,0) = c2;
    P_v(3,3) = P_v(4,4) = P_v(5,5) = c3;

    H.resize(num_slyde_sys * 2, num_slyde_sys * 2);
    H.setZero();

    if (hard_mode == 2) {
        file >> A >> q >> h0 >> ts >> alfa;
        for (int i = 0; i < num_slyde_sys * 2; i++)
            for (int j = 0; j < num_slyde_sys * 2; j++)
                H(j, i) = (i == j) ? A : A * q;
    }

    file >> tc0;
    file >> strain_rate_infl;
    file >> n_power;
}
