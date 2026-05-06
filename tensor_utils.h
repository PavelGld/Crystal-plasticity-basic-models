#pragma once
#include <Eigen/Dense>
#include <utility>

// Общие типы для Voigt-нотации (6-мерное пространство)
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;

// Voigt mapping: tensor (i,j) -> Voigt index
//   0=00, 1=11, 2=22, 3=23/32, 4=02/20, 5=01/10
inline int tensorToVoigtIndex(int i, int j) {
    static const int voigt_index[3][3] = {
        {0, 5, 4},
        {5, 1, 3},
        {4, 3, 2}
    };
    return voigt_index[i][j];
}

// Voigt index -> tensor (i,j) pair
inline std::pair<int, int> voigtToTensorIndex(int voigt_idx) {
    static const std::pair<int, int> tensor_index[6] = {
        {0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}
    };
    return tensor_index[voigt_idx];
}

// Tensor (3x3) -> Voigt (6x1), symmetric sum for off-diagonal
inline void tensorToVoigtInPlace(const Eigen::Matrix3d& A, Eigen::Matrix<double, 6, 1>& v) {
    v(0) = A(0, 0);
    v(1) = A(1, 1);
    v(2) = A(2, 2);
    v(3) = A(1, 2) + A(2, 1);
    v(4) = A(0, 2) + A(2, 0);
    v(5) = A(0, 1) + A(1, 0);
}

// Voigt (6x1) -> Tensor (3x3), symmetric reconstruction
inline void voigtToTensorInPlace(const Eigen::Matrix<double, 6, 1>& v, Eigen::Matrix3d& A) {
    A(0, 0) = v(0);
    A(1, 1) = v(1);
    A(2, 2) = v(2);
    A(1, 2) = A(2, 1) = v(3) * 0.5;
    A(0, 2) = A(2, 0) = v(4) * 0.5;
    A(0, 1) = A(1, 0) = v(5) * 0.5;
}

// F = exp(A) * F   using Taylor series when ||A|| is small
// Uses scratch as workspace (3x3 buffer, will be overwritten)
inline void applyExpTaylor(const Eigen::Matrix3d& A, Eigen::Matrix3d& F, Eigen::Matrix3d& scratch) {
    scratch.setIdentity();
    scratch.noalias() += A;
    Eigen::Matrix3d A2 = A * A;
    scratch.noalias() += A2 * 0.5;
    F.noalias() = scratch * F;
}


