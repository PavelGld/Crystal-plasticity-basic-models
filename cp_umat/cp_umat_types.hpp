#pragma once
#include <Eigen/Dense>

using Matrix6d = Eigen::Matrix<double, 6, 6>;

// Abaqus NTENS=6: [11,22,33,12,13,23] — tensor components
// Our Voigt:       [00,11,22,23,13,12] — engineering (shear ?2)

inline const Matrix6d& A_to_V()
{
    static const Matrix6d M = (Matrix6d() <<
        1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 2,
        0, 0, 0, 0, 2, 0,
        0, 0, 0, 2, 0, 0
    ).finished();
    return M;
}

inline const Matrix6d& V_to_A()
{
    static const Matrix6d M = (Matrix6d() <<
        1, 0, 0, 0,   0,   0,
        0, 1, 0, 0,   0,   0,
        0, 0, 1, 0,   0,   0,
        0, 0, 0, 0,   0, 0.5,
        0, 0, 0, 0, 0.5,   0,
        0, 0, 0, 0.5, 0,   0
    ).finished();
    return M;
}
