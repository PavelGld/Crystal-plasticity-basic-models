#include "cp_model/tensor_utils.hpp"
#include "cp_model/Grain.hpp"
#include "cp_model/MaterialParams.hpp"
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
    Eigen::Quaterniond q(R);
    Grain gr(params, q);

    double eps_dot = 0.0017;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0,0) = -0.5*eps_dot; L(1,1) = -0.5*eps_dot; L(2,2) = eps_dot;
    gr.Simulate_Step(1e-3, L);

    Matrix6d DDSDDE;
    gr.get_elasticity_lab_into(DDSDDE);

    // Reference: rotate cubic elasticity tensor directly
    // C'_ijkl = R_ia R_jb R_kc R_l d C_abcd
    double c11=106750, c12=60410, c44=28340;
    auto Cabcd = [&](int a,int b,int c,int d) -> double {
        if (a==b && b==c && c==d) return c11;
        if (a==b && c==d && a!=c) return c12;
        if (a==c && b==d && a!=b) return c44;
        if (a==d && b==c && a!=b) return c44;
        return 0;
    };
    double Crot[3][3][3][3];
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)for(int l=0;l<3;++l){
        double s=0;
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)for(int c=0;c<3;++c)for(int d=0;d<3;++d)
            s += R(i,a)*R(j,b)*R(k,c)*R(l,d)*Cabcd(a,b,c,d);
        Crot[i][j][k][l]=s;
    }
    // Map to 6x6 Voigt (doubled shears convention)
    auto vidx = [](int i,int j)->int {
        static const int t[3][3]={{0,5,4},{5,1,3},{4,3,2}};
        return t[i][j];
    };
    Matrix6d Pref = Matrix6d::Zero();
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){
        int I=vidx(i,j); double wi=(i==j)?1.0:2.0;
        for(int k=0;k<3;++k)for(int l=0;l<3;++l){
            int J=vidx(k,l); double wj=(k==l)?1.0:2.0;
            Pref(I,J) += wi*Crot[i][j][k][l]/wj;
        }
    }
    double err = (DDSDDE - Pref).norm();
    std::cout << "|DDSDDE - Pref| = " << err << " MPa\n";
    std::cout << "Relative: " << err/Pref.norm() << "\n";
    return 0;
}
