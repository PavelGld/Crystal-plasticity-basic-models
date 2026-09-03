#include "cp_umat_interface.hpp"
#include <cmath>

// === Hardcoded FCC {111}<110> slip systems ===
static const double fcc_b[12][3] = {
    {0,1,1},{0,1,1},{0,-1,1},{0,-1,1},
    {-1,1,0},{-1,1,0},{1,1,0},{1,1,0},
    {1,0,1},{1,0,1},{-1,0,1},{-1,0,1}
};
static const double fcc_n[12][3] = {
    {1,-1,1},{1,1,-1},{1,1,1},{-1,1,1},
    {1,1,1},{1,1,-1},{-1,1,1},{1,-1,1},
    {1,1,-1},{-1,1,1},{1,1,1},{1,-1,1}
};

static void build_slip_systems(MaterialParams& p)
{
    int n_sys = 12;
    p.num_slip_sys = n_sys;
    int N = n_sys * 2;
    p.M.resize(N); p.Wpm.resize(N); p.Ms_v.resize(N); p.MsE_v.resize(N);
    for (int i = 0; i < n_sys; i++) {
        double bn[3], nn[3];
        double nb = std::sqrt(fcc_b[i][0]*fcc_b[i][0]+fcc_b[i][1]*fcc_b[i][1]+fcc_b[i][2]*fcc_b[i][2]);
        double nd = std::sqrt(fcc_n[i][0]*fcc_n[i][0]+fcc_n[i][1]*fcc_n[i][1]+fcc_n[i][2]*fcc_n[i][2]);
        for (int j=0;j<3;j++){bn[j]=fcc_b[i][j]/nb; nn[j]=fcc_n[i][j]/nd;}
        p.M[i*2].setZero(); p.Wpm[i*2].setZero();
        for (int j=0;j<3;j++) for(int k=0;k<3;k++){
            p.M[i*2](j,k) = bn[k] * nn[j];
            p.Wpm[i*2](j,k) = 0.5*(bn[j]*nn[k] - bn[k]*nn[j]);
        }
        Eigen::Matrix3d Ms = 0.5*(p.M[i*2] + p.M[i*2].transpose());
        tensorToVoigtInPlace(Ms, p.Ms_v[i*2]);
        p.MsE_v[i*2] = p.Ms_v[i*2]; p.MsE_v[i*2](3)*=2; p.MsE_v[i*2](4)*=2; p.MsE_v[i*2](5)*=2;
        for(int j=0;j<3;j++) bn[j] = -bn[j];
        p.M[i*2+1].setZero(); p.Wpm[i*2+1].setZero();
        for (int j=0;j<3;j++) for(int k=0;k<3;k++){
            p.M[i*2+1](j,k) = bn[k] * nn[j];
            p.Wpm[i*2+1](j,k) = 0.5*(bn[j]*nn[k] - bn[k]*nn[j]);
        }
        Ms = 0.5*(p.M[i*2+1] + p.M[i*2+1].transpose());
        tensorToVoigtInPlace(Ms, p.Ms_v[i*2+1]);
        p.MsE_v[i*2+1] = p.Ms_v[i*2+1]; p.MsE_v[i*2+1](3)*=2; p.MsE_v[i*2+1](4)*=2; p.MsE_v[i*2+1](5)*=2;
    }
}

// === init_params: standard 19-PROPS format, same as single-element ===
void CPUMAT::init_params(const double* PROPS, int NPROPS)
{
    // Parse hard_mode first to know expected NPROPS
    int hard_mode = -1;
    if (NPROPS >= 5) hard_mode = (int)PROPS[4];
    int expected = (hard_mode == 2) ? 19 : (hard_mode == 1) ? 16 : 14;
    if (NPROPS != expected) return;
    if (hard_mode < 0 || hard_mode > 2) return;

    auto p = new MaterialParams();
    p->hard_mode = hard_mode;
    p->rot_mode = (int)PROPS[5];
    p->integ_mode = (int)PROPS[6];
    p->tc0 = PROPS[7];
    p->strain_rate_infl = PROPS[8];
    p->n_power = PROPS[9];

    build_slip_systems(*p);

    double c11 = PROPS[1], c12 = PROPS[2], c44 = PROPS[3];
    p->P_v.setZero();
    p->P_v(0,0)=p->P_v(1,1)=p->P_v(2,2)=c11;
    p->P_v(0,1)=p->P_v(1,0)=p->P_v(1,2)=p->P_v(2,1)=p->P_v(0,2)=p->P_v(2,0)=c12;
    p->P_v(3,3)=p->P_v(4,4)=p->P_v(5,5)=2.0*c44;

    int N = p->num_slip_sys * 2;
    p->H.resize(N, N);
    if (hard_mode == 2) {
        p->A = PROPS[14]; p->q = PROPS[15]; p->h0 = PROPS[16]; p->ts = PROPS[17]; p->alfa = PROPS[18];
        p->H.setConstant(p->A * p->q);
        p->H.diagonal().setConstant(p->A);
    } else if (hard_mode == 1) {
        p->iso_a = PROPS[14]; p->iso_b = PROPS[15];
    }

    params_ = p;
}

Grain* CPUMAT::grain()
{
    if (!grain_ && params_)
        grain_ = new Grain(params_, Eigen::Quaterniond::Identity());
    return grain_;
}

void CPUMAT::load_state(Grain* g, const double* STATEV)
{
    if (!g) return;
    Vector6d sv;
    for (int i=0;i<6;i++) sv(i) = STATEV[i];
    g->set_stress_v(sv);
    g->set_orientation(Eigen::Quaterniond(STATEV[6], STATEV[7], STATEV[8], STATEV[9]));
    int N = (int)g->num_slip_systems();
    std::vector<double> tc(N), shifts(N);
    for (int i=0;i<N;i++) tc[i]=STATEV[10+i];
    for (int i=0;i<N;i++) shifts[i]=STATEV[34+i];
    g->set_tc_all(tc); g->set_shifts_all(shifts);
    g->update_H_KM_if_needed();  // recompute H from restored tc (hard_mode==2 only)
    g->set_sim_time(STATEV[67]); g->set_acc_def(STATEV[68]);
}

void CPUMAT::save_state(const Grain* g, double* STATEV)
{
    if (!g) return;
    const Vector6d& sv = g->get_stress();
    for (int i=0;i<6;i++) STATEV[i]=sv(i);
    Eigen::Quaterniond q = g->get_orientation();
    STATEV[6]=q.w(); STATEV[7]=q.x(); STATEV[8]=q.y(); STATEV[9]=q.z();
    int N = (int)g->num_slip_systems();
    const auto& tc=g->get_tc_all(), &sh=g->get_shifts_all();
    for (int i=0;i<N;i++){ STATEV[10+i]=tc[i]; STATEV[34+i]=sh[i]; }
    STATEV[67]=g->get_time(); STATEV[68]=g->get_acc_def();
}

void CPUMAT::step(Grain* g, const double* DSTRAN, double DTIME, double* STRESS, double* DDSDDE, double& PNEWDT)
{
    PNEWDT = 1.0;
    if (!g) return;
    if (DTIME <= 0.0) {
        // Abaqus may probe the material with a zero time increment (initial stiffness):
        // return current stress in lab frame and the elastic lab stiffness, no integration.
        Eigen::Matrix3d sl0 = g->get_stress_lab();
        Vector6d sv0; tensorToVoigtInPlace(sl0, sv0);
        sv0 = (V_to_A() * sv0).eval();
        for (int i=0;i<6;i++) STRESS[i]=sv0(i);
        Matrix6d E; g->get_elasticity_lab_into(E);
        Matrix6d D0 = V_to_A() * E * A_to_V();
        for (int i=0;i<6;i++) for (int j=0;j<6;j++) DDSDDE[i+j*6]=D0(i,j);
        return;
    }

    Eigen::Matrix<double,6,1> de;
    for (int i=0;i<6;i++) de(i)=DSTRAN[i];
    Vector6d dv; dv.noalias() = A_to_V() * de;
    Eigen::Matrix3d Dt; voigtToTensorInPlace(dv, Dt);
    g->Simulate_Step(DTIME, Dt / DTIME);

    Eigen::Matrix3d sl = g->get_stress_lab();
    Vector6d sv; tensorToVoigtInPlace(sl, sv);
    sv = (V_to_A() * sv).eval();   // FIX: aliasing — noalias() при самоприсваивании портил sigma23
    for (int i=0;i<6;i++) STRESS[i]=sv(i);

    if (g->get_step_failed())
        PNEWDT = 0.25;   // local integrator failed: ask Abaqus to cut the increment

    Matrix6d T; g->get_tangent_lab_into(T);   // on failure: elastic fallback (tangent_valid=false)
    Matrix6d D = V_to_A() * T * A_to_V();
    for (int i=0;i<6;i++) for (int j=0;j<6;j++) DDSDDE[i+j*6]=D(i,j);
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
    int& NOEL, int& NPT, int& LAYER, int& KSPT, int* JSTEP, int& KINC)
{
    static CPUMAT cpumat;
    static bool ready = false;
    if (!ready) { cpumat.init_params(PROPS, NPROPS); ready = true; }

    Grain* g = cpumat.grain();

    bool first = (std::abs(STATEV[69] - 123456.0) > 0.5);
    if (first) {
        // Initial orientation from PROPS[10-13] for this element
        Eigen::Quaterniond q(PROPS[10], PROPS[11], PROPS[12], PROPS[13]);
        q.normalize();
        g->set_orientation(q);
        g->reset_state();  // clean state for this IP (tc=tc0, stress=0, shifts=0)
    CPUMAT::save_state(g, STATEV);
    // Diagnostic: store orientation quaternion in STATEV[60-63]
    Eigen::Quaterniond qdiag = g->get_orientation();
    STATEV[60] = qdiag.w(); STATEV[61] = qdiag.x(); STATEV[62] = qdiag.y(); STATEV[63] = qdiag.z();
    STATEV[64] = (double)NOEL;
        STATEV[69] = 123456.0;
    } else {
        CPUMAT::load_state(g, STATEV);
    }

    cpumat.step(g, DSTRAN, *DTIME, STRESS, DDSDDE, PNEWDT);
    CPUMAT::save_state(g, STATEV);
}

extern "C" void umat(
    double* STRESS, double* STATEV, double* DDSDDE,
    double& SSE,  double& SPD,  double& SCD,  double& RPL,
    double* DDSDDT, double* DRPLDE, double* DRPLDT,
    double* STRAN, double* DSTRAN, double* TIME,  double* DTIME,
    double& TEMP,  double& DTEMP, double* PREDEF, double* DPRED,
    char*   CMNAME, int& NDI,   int& NSHR,   int& NTENS, int& NSTATV,
    double* PROPS, int& NPROPS, double* COORDS, double* DROT,
    double& PNEWDT, double* CELENT, double* DFGRD0, double* DFGRD1,
    int& NOEL, int& NPT, int& LAYER, int& KSPT, int* JSTEP, int& KINC)
{
    umat_(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
          DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME,
          TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
          PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,
          NOEL, NPT, LAYER, KSPT, JSTEP, KINC);
}
