// Standalone check: DDSDDE consistency for random orientation
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
    Eigen::Quaterniond q(R);

    double eps_dot = 0.0017;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0,0) = -0.5*eps_dot; L(1,1) = -0.5*eps_dot; L(2,2) = eps_dot;
    double dt = 1e-3;
    double h = 1e-8;

    Eigen::Matrix3d Lp = L; Lp(2,2) += h/dt;

    // Normal step from fresh state
    Grain gr_n("grain_implicit.txt");
    gr_n.set_orientation(q);
    gr_n.Simulate_Step(dt, L);
    Vector6d s_n = gr_n.get_stress();
    Matrix6d DDSDDE;
    gr_n.get_tangent_lab_into(DDSDDE);

    // Perturbed step from fresh state
    Grain gr_p("grain_implicit.txt");
    gr_p.set_orientation(q);
    gr_p.Simulate_Step(dt, Lp);
    Vector6d s_p = gr_p.get_stress();

    Vector6d dstress_fd = (s_p - s_n);

    // Voigt strain increment difference
    Eigen::Matrix3d dL = Eigen::Matrix3d::Zero(); dL(2,2) = h/dt;
    Eigen::Matrix3d dD = 0.5*(dL + dL.transpose());
    Vector6d dev;
    tensorToVoigtInPlace(dD, dev);
    Vector6d dstress_tangent = DDSDDE * dev;

    double err = (dstress_fd - dstress_tangent).norm();
    std::cout << "FD stress change:   " << dstress_fd.transpose() << "\n";
    std::cout << "Tangent prediction: " << dstress_tangent.transpose() << "\n";
    std::cout << "FD norm: " << dstress_fd.norm() << "  Tangent norm: " << dstress_tangent.norm() << "\n";
    std::cout << "Relative error: " << err / std::max(dstress_fd.norm(), 1e-10) << "\n";
    return 0;
}