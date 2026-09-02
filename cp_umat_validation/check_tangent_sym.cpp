// Check if lab-frame tangent is symmetric for random orientation
#include "cp_model/Grain.hpp"
#include <iostream>
#include <random>

int main() {
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> ud(-1,1);
    Eigen::Matrix3d R;
    do {
        Eigen::Matrix3d A; for (int i=0;i<3;++i)for(int j=0;j<3;++j)A(i,j)=ud(gen);
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(A, Eigen::ComputeFullU|Eigen::ComputeFullV);
        R = svd.matrixU()*svd.matrixV().transpose();
    } while (R.determinant()<0);
    R.col(2) = R.col(0).cross(R.col(1));

    MaterialParamsPtr params = std::make_shared<MaterialParams>("grain_implicit.txt");
    Grain gr(params, Eigen::Quaterniond(R));

    double eps_dot = 0.0017;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0,0) = -0.5*eps_dot; L(1,1) = -0.5*eps_dot; L(2,2) = eps_dot;
    gr.Simulate_Step(1e-3, L);

    Matrix6d T;
    gr.get_tangent_lab_into(T);

    double asym = 0;
    for (int i=0;i<6;++i) for (int j=0;j<6;++j)
        asym += std::abs(T(i,j) - T(j,i));

    std::cout << "Tangent asymmetry: " << asym << " MPa\n";
    if (asym > 1e-6) {
        std::cout << "TANGENT IS NON-SYMMETRIC — this breaks Abaqus global Newton without UNSYMM\n";
    } else {
        std::cout << "Tangent is symmetric\n";
    }
    return 0;
}