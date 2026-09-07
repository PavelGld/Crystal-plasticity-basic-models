// Equivalence test: standalone Grain with full L (D+W)  vs  UMAT path with Abaqus-style DSTRAN + DROT
#include "cp_umat/cp_umat_interface.hpp"
#include <unsupported/Eigen/MatrixFunctions>
#include <cstdio>
#include <cstring>
const MaterialParams* cp_params();
extern "C" void umat_(double*,double*,double*,double&,double&,double&,double&,double*,double*,double*,double*,double*,double*,double*,double&,double&,double*,double*,char*,int&,int&,int&,int&,double*,int&,double*,double*,double&,double*,double*,double*,int&,int&,int&,int&,int*,int&);
int main(int argc,char**argv){
    double npow = argc>1?atof(argv[1]):20;
    int ninc = argc>2?atoi(argv[2]):300;
    double PROPS[19]={0,106750,60410,28340,2,2,1,6, 0.588235,npow, -0.12523701,0.09397545,0.59307825,0.78977370, 1,2,115,34,2.25};
    CPUMAT::init_params(PROPS,19);
    // --- reference: standalone grain, full velocity gradient of simple shear L = gdot * e1 (x) e2 (D+W)
    Eigen::Quaterniond q0(PROPS[10],PROPS[11],PROPS[12],PROPS[13]); q0.normalize();
    Grain ref(cp_params(), q0);
    double gdot=0.1, T=3.0, dt=T/ninc;   // gamma = 0.3
    Eigen::Matrix3d L=Eigen::Matrix3d::Zero(); L(0,1)=gdot;
    for(int k=0;k<ninc;k++) ref.Simulate_Step(dt, L);
    Eigen::Matrix3d s_ref = ref.get_stress_lab(); Eigen::Quaterniond q_ref = ref.get_orientation();

    // --- UMAT emulation (NLGEOM=YES style): per increment dF = exp(L dt) = dV * dR
    double STATEV[70]; std::memset(STATEV,0,sizeof(STATEV));
    double S[6]={0},D[36],dummy=0,z=0,pn; int i3=3,i6=6,np=19,one=1,jstep[4]={1,0,0,0}; char nm[4]="CP";
    Eigen::Matrix3d dF = (L*dt).exp();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(dF, Eigen::ComputeFullU|Eigen::ComputeFullV);
    Eigen::Matrix3d dR = svd.matrixU()*svd.matrixV().transpose();
    Eigen::Matrix3d dV = dF*dR.transpose();                       // left stretch: dF = dV dR
    Eigen::Matrix3d logV = dV.log();                              // log strain increment in rotated (end) basis
    double DSTRAN[6]={logV(0,0),logV(1,1),logV(2,2),2*logV(0,1),2*logV(0,2),2*logV(1,2)}; // Abaqus order, eng. shear
    double DROT[9]; for(int i=0;i<3;i++)for(int j=0;j<3;j++) DROT[i+j*3]=(getenv("NO_DROT")? (i==j?1.0:0.0) : dR(i,j));
    for(int k=1;k<=ninc;k++)
        umat_(S,STATEV,D,dummy,dummy,dummy,dummy,0,0,0,0,DSTRAN,&z,&dt,dummy,dummy,0,0,nm,i3,i3,i6,np,PROPS,np,0,DROT,pn,0,0,0,one,one,one,one,jstep,k);
    Eigen::Quaterniond q_um(STATEV[6],STATEV[7],STATEV[8],STATEV[9]);
    // Abaqus stress is in the corotational basis rotated by R_total = dR^ninc ; bring back to lab for comparison
    Eigen::Matrix3d Rtot = Eigen::Matrix3d::Identity(); for(int k=0;k<ninc;k++) Rtot = dR*Rtot;
    Eigen::Matrix3d s_um_corot; s_um_corot<<S[0],S[3],S[4], S[3],S[1],S[5], S[4],S[5],S[2];
    // Abaqus: sigma_corot = Rtot^T sigma_lab Rtot? Our returned sigma_lab already uses R_lat rotated by DROTs,
    // i.e. it IS expressed in the lab frame co-moving with Abaqus' basis -> compare directly.
    Eigen::Matrix3d s_um = s_um_corot;
    double dq = std::min((q_ref.coeffs()-q_um.coeffs()).norm(), (q_ref.coeffs()+q_um.coeffs()).norm());
    printf("n=%g ninc=%d gamma=%.2f\n  |sigma_ref - sigma_umat| / |sigma_ref| = %.2e   (|sigma|=%.1f MPa)\n  |q_ref - q_umat| = %.2e\n",
           npow, ninc, gdot*T, (s_ref-s_um).norm()/s_ref.norm(), s_ref.norm(), dq);
    Eigen::AngleAxisd aa(q0.conjugate()*q_ref); printf("  lattice rotated by %.2f deg (ref)\n", aa.angle()*180/M_PI);
    return 0;
}
