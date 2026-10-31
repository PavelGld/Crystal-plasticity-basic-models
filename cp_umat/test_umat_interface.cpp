// Interface tests for the UMAT: component order, DROT objectivity, energy
// balance, and large-strain consistency with the standalone Grain.
//
// Build/run: cp_umat\run_test_umat_interface.bat <material>
// Exit code 0 = PASS.

#include "cp_model/Grain.hpp"
#include "cp_umat_types.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>

extern "C" void umat_(
    double* STRESS, double* STATEV, double* DDSDDE,
    double& SSE, double& SPD, double& SCD, double& RPL,
    double* DDSDDT, double* DRPLDE, double* DRPLDT,
    double* STRAN, double* DSTRAN, double* TIME, double* DTIME,
    double& TEMP, double& DTEMP, double* PREDEF, double* DPRED,
    char* CMNAME, int& NDI, int& NSHR, int& NTENS, int& NSTATV,
    double* PROPS, int& NPROPS, double* COORDS, double* DROT,
    double& PNEWDT, double* CELENT, double* DFGRD0, double* DFGRD1,
    int& NOEL, int& NPT, int& LAYER, int& KSPT, int* JSTEP, int& KINC);

namespace {

constexpr int NV = 70;
constexpr double PI = 3.14159265358979323846;

struct UmatState {
    double PROPS[19];
    double STATEV[NV];
    double STRESS[6];
    double DDSDDE[36];
    double SSE = 0, SPD = 0, SCD = 0, RPL = 0;
    double STRAN[6] = {0,0,0,0,0,0};
    double TIME = 0.0;
    int nprop = 19;
};

void fill_props(double* P, const MaterialParams& p, const Eigen::Quaterniond& q)
{
    P[0] = 0.0;
    P[1] = p.P_v(0, 0); P[2] = p.P_v(0, 1); P[3] = 0.5 * p.P_v(3, 3);
    P[4] = (double)p.hard_mode; P[5] = (double)p.rot_mode; P[6] = (double)p.integ_mode;
    P[7] = p.tc0; P[8] = p.strain_rate_infl; P[9] = p.n_power;
    P[10] = q.w(); P[11] = q.x(); P[12] = q.y(); P[13] = q.z();
    P[14] = p.A; P[15] = p.q; P[16] = p.h0; P[17] = p.ts; P[18] = p.alfa;
}

void init_state(UmatState& u, const MaterialParams& p, const Eigen::Quaterniond& q)
{
    std::memset(&u, 0, sizeof(u));
    u.nprop = (p.hard_mode == 2) ? 19 : (p.hard_mode == 1) ? 16 : 14;
    fill_props(u.PROPS, p, q);
}

void step_umat(UmatState& u, const double dstr[6], double dt, const double drot[9])
{
    double DSTRAN[6], DROT[9];
    for (int i = 0; i < 6; i++) DSTRAN[i] = dstr[i];
    for (int i = 0; i < 9; i++) DROT[i] = drot[i];

    double DTIME = dt;
    double TIME[2] = {u.TIME, u.TIME + dt};
    double DDSDDT[6] = {0}, DRPLDE[6] = {0}, DRPLDT[1] = {0};
    double TEMP = 0, DTEMP = 0, PREDEF[1] = {0}, DPRED[1] = {0};
    char CMNAME[80] = "CPFCC";
    int NDI = 3, NSHR = 3, NTENS = 6, NSTATV = NV, NPROPS = u.nprop;
    double COORDS[3] = {0,0,0};
    double PNEWDT = 1e10, CELENT = 1.0;
    double DFGRD0[9] = {1,0,0,0,1,0,0,0,1}, DFGRD1[9] = {1,0,0,0,1,0,0,0,1};
    int NOEL = 1, NPT = 1, LAYER = 1, KSPT = 1, KINC = 0;
    int JSTEP[4] = {1,1,1,1};

    umat_(u.STRESS, u.STATEV, u.DDSDDE, u.SSE, u.SPD, u.SCD, u.RPL,
          DDSDDT, DRPLDE, DRPLDT, u.STRAN, DSTRAN, TIME, &DTIME,
          TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
          u.PROPS, NPROPS, COORDS, DROT, PNEWDT, &CELENT, DFGRD0, DFGRD1,
          NOEL, NPT, LAYER, KSPT, JSTEP, KINC);
    u.TIME += dt;
    for (int i = 0; i < 6; i++) u.STRAN[i] += dstr[i];
}

void identity_rot(double R[9]) { for (int i = 0; i < 9; i++) R[i] = (i % 4 == 0) ? 1.0 : 0.0; }

void rot_z(double angle, double R[9])
{
    double c = std::cos(angle), s = std::sin(angle);
    double M[9] = {c,-s,0, s,c,0, 0,0,1};   // column-major
    for (int i = 0; i < 9; i++) R[i] = M[i];
}

Eigen::Quaterniond state_quat(const UmatState& u)
{
    Eigen::Quaterniond q(u.STATEV[6], u.STATEV[7], u.STATEV[8], u.STATEV[9]);
    q.normalize();
    return q;
}

double max_abs6(const double* a, const double* b)
{
    double m = 0;
    for (int i = 0; i < 6; i++) m = std::max(m, std::fabs(a[i] - b[i]));
    return m;
}

// ---------------------------------------------------------------- ordering
bool test_ordering(const MaterialParams& p)
{
    const double c11 = p.P_v(0,0), c12 = p.P_v(0,1), c44 = 0.5 * p.P_v(3,3);
    const double d = 1e-6;
    const double I9[9] = {1,0,0,0,1,0,0,0,1};

    struct Case { int comp; double expect[6]; };
    Case cases[6] = {
        {0, {c11*d, c12*d, c12*d, 0, 0, 0}},
        {1, {c12*d, c11*d, c12*d, 0, 0, 0}},
        {2, {c12*d, c12*d, c11*d, 0, 0, 0}},
        {3, {0, 0, 0, c44*d, 0, 0}},
        {4, {0, 0, 0, 0, c44*d, 0}},
        {5, {0, 0, 0, 0, 0, c44*d}},
    };

    bool ok = true;
    for (const auto& c : cases) {
        UmatState u; init_state(u, p, Eigen::Quaterniond::Identity());
        double dstr[6] = {0,0,0,0,0,0}; dstr[c.comp] = d;
        step_umat(u, dstr, 1.0, I9);
        for (int i = 0; i < 6; i++) {
            double scale = std::max(std::fabs(c.expect[i]), 1e-3);
            if (std::fabs(u.STRESS[i] - c.expect[i]) > 1e-6 * scale) {
                std::printf("  ordering comp %d: STRESS[%d] = %.6e, expected %.6e\n",
                            c.comp, i, u.STRESS[i], c.expect[i]);
                ok = false;
            }
        }
    }

    // DDSDDE in Abaqus order: (0,0)=c11, (1,0)=c12, shear diagonal = c44
    UmatState u; init_state(u, p, Eigen::Quaterniond::Identity());
    double dstr[6] = {d,0,0,0,0,0};
    step_umat(u, dstr, 1.0, I9);
    auto chk = [&](int i, int j, double e) {
        double got = u.DDSDDE[i + j*6];
        double scale = std::max(std::fabs(e), 1e3);
        if (std::fabs(got - e) > 1e-6 * scale) {
            std::printf("  ordering DDSDDE(%d,%d) = %.6e, expected %.6e\n", i, j, got, e);
            ok = false;
        }
    };
    chk(0,0,c11); chk(1,0,c12); chk(2,0,c12);
    chk(3,3,c44); chk(4,4,c44); chk(5,5,c44);

    std::printf("  ordering: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ------------------------------------------------------------- objectivity
bool test_objectivity(const MaterialParams& p)
{
    const double I9[9] = {1,0,0,0,1,0,0,0,1};
    bool ok = true;

    for (double deg : {30.0, 90.0, 179.0}) {
        UmatState u; init_state(u, p, Eigen::Quaterniond::Identity());
        double uniax[6] = {-0.5e-3, -0.5e-3, 1.0e-3, 0, 0, 0};
        for (int i = 0; i < 100; i++) step_umat(u, uniax, 1.0, I9);

        double s_before[6]; for (int i = 0; i < 6; i++) s_before[i] = u.STATEV[i];
        double lab_before[6]; for (int i = 0; i < 6; i++) lab_before[i] = u.STRESS[i];
        Eigen::Quaterniond q_before = state_quat(u);

        double R[9]; rot_z(deg * PI / 180.0, R);
        double zero[6] = {0,0,0,0,0,0};
        step_umat(u, zero, 1.0, R);

        double ds = max_abs6(u.STATEV, s_before);
        if (ds > 1e-3) { std::printf("  objectivity %.0f deg: crystal stress changed by %.3e\n", deg, ds); ok = false; }

        Eigen::Matrix3d Rm;
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) Rm(i,j) = R[i+j*3];
        Eigen::Quaterniond q_rot(Rm);
        Eigen::Quaterniond q_expect = q_rot * q_before;
        double qd = std::fabs(q_expect.dot(state_quat(u)));
        if (qd < 1.0 - 1e-10) { std::printf("  objectivity %.0f deg: quaternion mismatch (dot=%.12f)\n", deg, qd); ok = false; }

        // lab stress must rotate: S' = R S R^T
        Eigen::Matrix3d Sl;
        Sl << lab_before[0], lab_before[3], lab_before[4],
              lab_before[3], lab_before[1], lab_before[5],
              lab_before[4], lab_before[5], lab_before[2];
        Eigen::Matrix3d Sr = Rm * Sl * Rm.transpose();
        double exp6[6] = {Sr(0,0), Sr(1,1), Sr(2,2), Sr(0,1), Sr(0,2), Sr(1,2)};
        double dl = max_abs6(u.STRESS, exp6);
        if (dl > 1e-3) { std::printf("  objectivity %.0f deg: lab stress diff %.3e\n", deg, dl); ok = false; }

        // inverse rotation returns the orientation and keeps the stress
        double Rinv[9]; rot_z(-deg * PI / 180.0, Rinv);
        step_umat(u, zero, 1.0, Rinv);
        double qb = std::fabs(q_before.dot(state_quat(u)));
        if (qb < 1.0 - 1e-10) { std::printf("  objectivity %.0f deg: inverse rotation mismatch (dot=%.12f)\n", deg, qb); ok = false; }
        if (max_abs6(u.STATEV, s_before) > 1e-3) { std::printf("  objectivity %.0f deg: stress not restored\n", deg); ok = false; }
    }

    std::printf("  objectivity: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ----------------------------------------------------------------- energy
bool test_energy(const MaterialParams& p)
{
    const double I9[9] = {1,0,0,0,1,0,0,0,1};
    UmatState u; init_state(u, p, Eigen::Quaterniond::Identity());

    double uniax[6] = {-0.5e-3, -0.5e-3, 1.0e-3, 0, 0, 0};
    double sigma_prev[6] = {0,0,0,0,0,0};
    double W = 0.0;
    for (int i = 0; i < 100; i++) {
        step_umat(u, uniax, 1.0, I9);
        for (int j = 0; j < 6; j++)
            W += 0.5 * (sigma_prev[j] + u.STRESS[j]) * uniax[j];
        for (int j = 0; j < 6; j++) sigma_prev[j] = u.STRESS[j];
    }

    double Eint = u.SSE + u.SPD;
    double rel = std::fabs(Eint - W) / std::max(std::fabs(W), 1e-30);
    bool ok = (u.SPD > 0.0) && (u.SSE > 0.0) && (rel < 1e-6);

    std::printf("  energy: W=%.9e  SSE=%.9e  SPD=%.9e  |SSE+SPD-W|/W=%.2e  %s\n",
                W, u.SSE, u.SPD, rel, ok ? "PASS" : "FAIL");
    return ok;
}

// ------------------------------------------------------------ large strain
bool test_large_strain(const MaterialParams& p)
{
    const Eigen::Quaterniond q0 = Eigen::Quaterniond::Identity();
    Grain g(&p, q0);
    UmatState u; init_state(u, p, q0);

    const double dt = 1e-3;
    const int n = 1000;
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L(0,1) = 1.0;

    Eigen::Matrix3d Finc = (L * dt).exp();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(Finc, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d Rinc = svd.matrixU() * svd.matrixV().transpose();
    double R[9];
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) R[i+j*3] = Rinc(i,j);

    Eigen::Matrix3d Dm = 0.5 * (L + L.transpose());
    Eigen::Matrix3d Dinc = Rinc * Dm * Rinc.transpose();
    Vector6d dvp; tensorToVoigtInPlace(Dinc, dvp);
    Vector6d dva = A_to_V() * dvp;
    double dstr[6]; for (int i = 0; i < 6; i++) dstr[i] = dva(i) * dt;

    for (int i = 0; i < n; i++) {
        g.Simulate_Step(dt, L);
        step_umat(u, dstr, dt, R);
    }

    Eigen::Matrix3d sl = g.get_stress_lab();
    Vector6d sv; tensorToVoigtInPlace(sl, sv);
    Vector6d sa = V_to_A() * sv;
    double ds = 0.0;
    for (int i = 0; i < 6; i++) ds = std::max(ds, std::fabs(u.STRESS[i] - sa(i)));

    Eigen::Quaterniond qs = g.get_orientation();
    double qd = std::fabs(qs.dot(state_quat(u)));

    double scale = 1.0;
    for (int i = 0; i < 6; i++) scale = std::max(scale, std::fabs(sa(i)));
    bool ok = (ds < 1e-6 * scale) && (qd > 1.0 - 1e-8);
    std::printf("  large strain (shear 100%%): stress diff %.3e Pa (rel %.1e), quat dot %.10f  %s\n",
                ds, ds / scale, qd, ok ? "PASS" : "FAIL");
    return ok;
}

} // namespace

int main(int argc, char** argv)
{
    const char* path = (argc > 1) ? argv[1] : "grain_implicit.txt";
    MaterialParams p(path);
    std::printf("UMAT interface tests, material %s\n", path);

    bool ok = true;
    ok = test_ordering(p) && ok;
    ok = test_objectivity(p) && ok;
    ok = test_energy(p) && ok;
    ok = test_large_strain(p) && ok;

    std::printf("%s\n", ok ? "UMAT_INTERFACE_PASS" : "UMAT_INTERFACE_FAIL");
    return ok ? 0 : 1;
}
