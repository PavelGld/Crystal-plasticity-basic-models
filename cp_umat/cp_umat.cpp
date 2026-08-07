#include "cp_umat_interface.hpp"
#include <cmath>

// === init_from_props ===
void CPUMAT::init_from_props(const double* PROPS, int NPROPS)
{
    // Minimal sanity: check NPROPS against hard_mode
    int hard_mode = (int)PROPS[4];
    int expected = (hard_mode == 2) ? 19 : (hard_mode == 1) ? 16 : 14;
    if (NPROPS != expected) return;   // wrong PROPS count — skip init

    double c11 = PROPS[1], c12 = PROPS[2], c44 = PROPS[3];
    int rot_mode = (int)PROPS[5];
    int integ_mode = (int)PROPS[6];
    double tc0 = PROPS[7];
    double srf = PROPS[8];
    double np  = PROPS[9];
    double qw = PROPS[10], qx = PROPS[11], qy = PROPS[12], qz = PROPS[13];

    raw_params_ = new MaterialParams();
    MaterialParams& p = *raw_params_;
    p.hard_mode = hard_mode;
    p.rot_mode = rot_mode;
    p.integ_mode = integ_mode;
    p.tc0 = tc0;
    p.strain_rate_infl = srf;
    p.n_power = np;

    // FCC {111}<110>, 12 slip systems
    static const double b[12][3] = {
        {0,1,1},{0,1,1},{0,-1,1},{0,-1,1},
        {-1,1,0},{-1,1,0},{1,1,0},{1,1,0},
        {1,0,1},{1,0,1},{-1,0,1},{-1,0,1}
    };
    static const double n[12][3] = {
        {1,-1,1},{1,1,-1},{1,1,1},{-1,1,1},
        {1,1,1},{1,1,-1},{-1,1,1},{1,-1,1},
        {1,1,-1},{-1,1,1},{1,1,1},{1,-1,1}
    };

    int n_sys = 12;
    p.num_slip_sys = n_sys;
    int N = n_sys * 2;
    p.M.resize(N);
    p.Wpm.resize(N);
    p.Ms_v.resize(N);
    p.MsE_v.resize(N);

    for (int i = 0; i < n_sys; i++) {
        double bn[3], nn[3];
        double nb = std::sqrt(b[i][0]*b[i][0]+b[i][1]*b[i][1]+b[i][2]*b[i][2]);
        double nd = std::sqrt(n[i][0]*n[i][0]+n[i][1]*n[i][1]+n[i][2]*n[i][2]);
        for (int j=0;j<3;j++){bn[j]=b[i][j]/nb; nn[j]=n[i][j]/nd;}

        // +b system
        p.M[i*2].setZero(); p.Wpm[i*2].setZero();
        for (int j=0;j<3;j++) for(int k=0;k<3;k++){
            p.M[i*2](j,k) = bn[k] * nn[j];
            p.Wpm[i*2](j,k) = 0.5*(bn[j]*nn[k] - bn[k]*nn[j]);
        }
        Eigen::Matrix3d Ms = 0.5*(p.M[i*2] + p.M[i*2].transpose());
        tensorToVoigtInPlace(Ms, p.Ms_v[i*2]);
        p.MsE_v[i*2] = p.Ms_v[i*2];
        p.MsE_v[i*2](3)*=2; p.MsE_v[i*2](4)*=2; p.MsE_v[i*2](5)*=2;

        // -b system
        for(int j=0;j<3;j++) bn[j] = -bn[j];
        p.M[i*2+1].setZero(); p.Wpm[i*2+1].setZero();
        for (int j=0;j<3;j++) for(int k=0;k<3;k++){
            p.M[i*2+1](j,k) = bn[k] * nn[j];
            p.Wpm[i*2+1](j,k) = 0.5*(bn[j]*nn[k] - bn[k]*nn[j]);
        }
        Ms = 0.5*(p.M[i*2+1] + p.M[i*2+1].transpose());
        tensorToVoigtInPlace(Ms, p.Ms_v[i*2+1]);
        p.MsE_v[i*2+1] = p.Ms_v[i*2+1];
        p.MsE_v[i*2+1](3)*=2; p.MsE_v[i*2+1](4)*=2; p.MsE_v[i*2+1](5)*=2;
    }

    // Elasticity matrix
    p.P_v.setZero();
    p.P_v(0,0)=p.P_v(1,1)=p.P_v(2,2)=c11;
    p.P_v(0,1)=p.P_v(1,0)=p.P_v(1,2)=p.P_v(2,1)=p.P_v(0,2)=p.P_v(2,0)=c12;
    p.P_v(3,3)=p.P_v(4,4)=p.P_v(5,5)=c44;

    // Hardening
    p.H.resize(N, N);
    if (hard_mode == 2) {
        p.A = PROPS[14]; p.q = PROPS[15]; p.h0 = PROPS[16]; p.ts = PROPS[17]; p.alfa = PROPS[18];
        p.H.setConstant(p.A * p.q);
        p.H.diagonal().setConstant(p.A);
    } else if (hard_mode == 1) {
        p.iso_a = PROPS[14]; p.iso_b = PROPS[15];
    }

    Eigen::Quaterniond q(qw, qx, qy, qz);
    q.normalize();
    grain_ = new Grain(raw_params_, q);
    initialized_ = true;
}

// === load_state ===
void CPUMAT::load_state(const double* STATEV)
{
    if (!grain_) return;
    Vector6d sv;
    for (int i=0;i<6;i++) sv(i) = STATEV[i];
    grain_->set_stress_v(sv);

    Eigen::Quaterniond q(STATEV[6], STATEV[7], STATEV[8], STATEV[9]);
    grain_->set_orientation(q);

    int N = (int)grain_->num_slip_systems(); // = 24 for FCC
    std::vector<double> tc(N), shifts(N);
    for(int i=0;i<N;i++) tc[i] = STATEV[10+i];
    for(int i=0;i<N;i++) shifts[i] = STATEV[34+i];
    grain_->set_tc_all(tc);
    grain_->set_shifts_all(shifts);

    grain_->set_sim_time(STATEV[67]);
    grain_->set_acc_def(STATEV[68]);
}

// === save_state ===
void CPUMAT::save_state(double* STATEV) const
{
    if (!grain_) return;
    const Vector6d& sv = grain_->get_stress();
    for(int i=0;i<6;i++) STATEV[i] = sv(i);

    Eigen::Quaterniond q = grain_->get_orientation();
    STATEV[6]=q.w(); STATEV[7]=q.x(); STATEV[8]=q.y(); STATEV[9]=q.z();

    int N = (int)grain_->num_slip_systems();
    const auto& tc = grain_->get_tc_all();
    const auto& sh = grain_->get_shifts_all();
    for(int i=0;i<N;i++){ STATEV[10+i] = tc[i]; STATEV[34+i] = sh[i]; }

    STATEV[67] = grain_->get_time();
    STATEV[68] = grain_->get_acc_def();
}

// === step ===
void CPUMAT::step(const double* DSTRAN, double DTIME, double* STRESS, double* DDSDDE, double& PNEWDT)
{
    PNEWDT = 1.0;
    if (!grain_) return;

    // DSTRAN (Abaqus tensor components) -> D_v (Voigt engineering)
    // Manual copy: Abaqus Fortran arrays not Eigen-aligned
    Eigen::Matrix<double,6,1> de_abaqus;
    for (int i = 0; i < 6; i++) de_abaqus(i) = DSTRAN[i];
    D_v_buf_.noalias() = A_to_V() * de_abaqus;

    // D_v is strain INCREMENT (delta-epsilon), not rate. Convert to tensor L:
    // L_tensor = voigtToTensor(D_v) / DTIME (pure strain, no spin -- NLGEOM=NO)
    Eigen::Matrix3d D_tensor;
    voigtToTensorInPlace(D_v_buf_, D_tensor);
    Eigen::Matrix3d L = D_tensor / DTIME;

    grain_->Simulate_Step(DTIME, L);

    // Stress: lab frame -> Voigt -> Abaqus
    Eigen::Matrix3d s_lab = grain_->get_stress_lab();
    Vector6d s_v;
    tensorToVoigtInPlace(s_lab, s_v);
    s_v.noalias() = V_to_A() * s_v;
    for (int i = 0; i < 6; i++) STRESS[i] = s_v(i);

    // DDSDDE: tangent in lab frame, Voigt -> Abaqus
    grain_->get_tangent_lab_into(T_buf_);
    // DDSDDE_abaqus = V_to_A * T * A_to_V, stored column-major (Fortran)
    Matrix6d D = V_to_A() * T_buf_ * A_to_V();
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            DDSDDE[i + j*6] = D(i, j);  // column-major (Fortran)
}

// === extern "C" UMAT entry ===
extern "C" void umat_(
    double* STRESS, double* STATEV, double* DDSDDE,
    double& SSE,  double& SPD,  double& SCD,  double& RPL,
    double* DDSDDT, double* DRPLDE, double* DRPLDT,
    double* STRAN, double* DSTRAN, double* TIME,  double* DTIME,
    double& TEMP,  double& DTEMP, double* PREDEF, double* DPRED,
    char*   CMNAME, int& NDI,   int& NSHR,   int& NTENS, int& NSTATV,
    double* PROPS, int& NPROPS, double* COORDS, double* DROT,
    double& PNEWDT, double* CELENT, double* DFGRD0, double* DFGRD1,
    int& NOEL, int& NPT, int& LAYER, int& KSPT, int& KINC)
{
    static CPUMAT cpumat;

    // Initialize on first call
    if (KINC == 1 && std::abs(STATEV[69] - 123456.0) > 0.5) {
        cpumat.init_from_props(PROPS, NPROPS);
        STATEV[69] = 123456.0;
    }

    cpumat.load_state(STATEV);
    cpumat.step(DSTRAN, *DTIME, STRESS, DDSDDE, PNEWDT);
    cpumat.save_state(STATEV);
}

// Provide undecorated name for Abaqus linker
extern "C" void umat(
    double* STRESS, double* STATEV, double* DDSDDE,
    double& SSE,  double& SPD,  double& SCD,  double& RPL,
    double* DDSDDT, double* DRPLDE, double* DRPLDT,
    double* STRAN, double* DSTRAN, double* TIME,  double* DTIME,
    double& TEMP,  double& DTEMP, double* PREDEF, double* DPRED,
    char*   CMNAME, int& NDI,   int& NSHR,   int& NTENS, int& NSTATV,
    double* PROPS, int& NPROPS, double* COORDS, double* DROT,
    double& PNEWDT, double* CELENT, double* DFGRD0, double* DFGRD1,
    int& NOEL, int& NPT, int& LAYER, int& KSPT, int& KINC)
{
    umat_(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
          DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME,
          TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
          PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,
          NOEL, NPT, LAYER, KSPT, KINC);
}
