// Regression: the UMAT path (cp_umat.cpp) must reproduce the standalone Grain
// exactly for the same strain increments, material and orientation.
//
// This is the test that would have caught the doubled-shear Ms_v construction:
// tau = sigma_v * Ms_v must use tensor shear components (no x2), while Dp uses
// the engineering MsE_v (x2). A wrong Ms_v shows up as a large stress/tangent
// mismatch for any non-diagonal crystal-frame stress (random orientation,
// simple shear), even though identity-orientation uniaxial tension passes.
//
// Build/run: cp_umat\run_test_umat_vs_grain.bat
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

constexpr double TOL_STRESS  = 1e-8;   // relative to max |sigma|
constexpr double TOL_TANGENT = 1e-6;   // relative to max |DDSDDE|

void fill_props(double* P, const MaterialParams& p, const Eigen::Quaterniond& q)
{
    P[0] = 0.0;
    P[1] = p.P_v(0, 0);          // c11
    P[2] = p.P_v(0, 1);          // c12
    P[3] = 0.5 * p.P_v(3, 3);    // c44
    P[4] = (double)p.hard_mode;
    P[5] = (double)p.rot_mode;
    P[6] = (double)p.integ_mode;
    P[7] = p.tc0;
    P[8] = p.strain_rate_infl;
    P[9] = p.n_power;
    P[10] = q.w(); P[11] = q.x(); P[12] = q.y(); P[13] = q.z();
    P[14] = p.A; P[15] = p.q; P[16] = p.h0; P[17] = p.ts; P[18] = p.alfa;
}

struct CaseResult {
    double stress_err;
    double tangent_err;
    double stress_scale;
};

CaseResult run_case(const MaterialParams& p, const Eigen::Quaterniond& q0,
                    const double dstr[6], int ninc)
{
    const int nprop = (p.hard_mode == 2) ? 19 : (p.hard_mode == 1) ? 16 : 14;
    double PROPS[19];
    fill_props(PROPS, p, q0);

    Grain g(&p, q0);

    double STRESS[6] = {0, 0, 0, 0, 0, 0};
    double STATEV[70];
    std::memset(STATEV, 0, sizeof(STATEV));
    double DDSDDE[36] = {0};
    double DSTRAN[6] = {0, 0, 0, 0, 0, 0}, STRAN[6] = {0, 0, 0, 0, 0, 0};
    double TIME[2] = {0.0, 0.0}, DTIME = 1.0;
    double SSE = 0, SPD = 0, SCD = 0, RPL = 0;
    double DDSDDT[6] = {0}, DRPLDE[6] = {0}, DRPLDT[1] = {0};
    double TEMP = 0, DTEMP = 0, PREDEF[1] = {0}, DPRED[1] = {0};
    char CMNAME[80] = "CPFCC";
    int NDI = 3, NSHR = 3, NTENS = 6, NSTATV = 70, NPROPS = nprop;
    double COORDS[3] = {0, 0, 0};
    double DROT[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    double PNEWDT = 1e10, CELENT = 1.0;
    double DFGRD0[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    double DFGRD1[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    int NOEL = 1, NPT = 1, LAYER = 1, KSPT = 1, KINC = 0;
    int JSTEP[4] = {1, 1, 1, 1};

    CaseResult r{0.0, 0.0, 1.0};
    for (int inc = 0; inc < ninc; inc++) {
        for (int i = 0; i < 6; i++) DSTRAN[i] = dstr[i];

        // standalone: same D in the lab frame
        Vector6d dv;
        for (int i = 0; i < 6; i++) dv(i) = DSTRAN[i];
        Vector6d dvv = A_to_V() * dv;
        Eigen::Matrix3d Dt;
        voigtToTensorInPlace(dvv, Dt);
        g.Simulate_Step(DTIME, Dt / DTIME);

        Eigen::Matrix3d sl = g.get_stress_lab();
        Vector6d sv;
        tensorToVoigtInPlace(sl, sv);
        Vector6d sa = V_to_A() * sv;

        TIME[1] = (double)(inc + 1);
        umat_(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
              DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, &DTIME,
              TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
              PROPS, NPROPS, COORDS, DROT, PNEWDT, &CELENT, DFGRD0, DFGRD1,
              NOEL, NPT, LAYER, KSPT, JSTEP, KINC);
        for (int i = 0; i < 6; i++) STRAN[i] += DSTRAN[i];

        for (int i = 0; i < 6; i++) {
            r.stress_scale = std::max(r.stress_scale, std::fabs(sa(i)));
            r.stress_err = std::max(r.stress_err, std::fabs(STRESS[i] - sa(i)));
        }

        Matrix6d Tlab;
        g.get_tangent_lab_into(Tlab);
        Matrix6d Da = V_to_A() * Tlab * A_to_V();
        double tscale = 1.0;
        for (int i = 0; i < 36; i++) tscale = std::max(tscale, std::fabs(Da(i / 6, i % 6)));
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                r.tangent_err = std::max(r.tangent_err,
                                         std::fabs(DDSDDE[i + j * 6] - Da(i, j)) / tscale);
    }
    r.stress_err /= r.stress_scale;
    return r;
}

} // namespace

int main(int argc, char** argv)
{
    const char* path = (argc > 1) ? argv[1] : "grain_implicit.txt";
    MaterialParams p(path);

    Eigen::Quaterniond q_id(1, 0, 0, 0);
    Eigen::Quaterniond q_rnd(0.61, -0.7639, 0.125, -0.1697);
    q_rnd.normalize();

    const double uniax[6] = {-0.5e-3, -0.5e-3, 1.0e-3, 0, 0, 0};
    const double shear[6] = {0, 0, 0, 1.0e-3, 0, 0};   // Abaqus order 12,13,23

    struct { const char* name; Eigen::Quaterniond q; const double* d; } cases[] = {
        {"uniaxial, identity", q_id,  uniax},
        {"uniaxial, random  ", q_rnd, uniax},
        {"shear,    identity", q_id,  shear},
        {"shear,    random  ", q_rnd, shear},
    };

    bool ok = true;
    std::printf("UMAT vs standalone Grain, material %s\n", path);
    for (const auto& c : cases) {
        CaseResult r = run_case(p, c.q, c.d, 100);
        bool pass = (r.stress_err < TOL_STRESS) && (r.tangent_err < TOL_TANGENT);
        ok = ok && pass;
        std::printf("  %s : stress %.2e  tangent %.2e  %s\n",
                    c.name, r.stress_err, r.tangent_err, pass ? "PASS" : "FAIL");
    }

    std::printf("%s\n", ok ? "UMAT_REGRESSION_PASS" : "UMAT_REGRESSION_FAIL");
    return ok ? 0 : 1;
}
