// FD check of DDSDDE against d(STRESS)/d(DSTRAN) through umat_, on a loading path
// usage: fd <deps> [q: i|r]
#include "cp_umat/cp_umat_interface.hpp"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <array>
extern "C" void umat_(double*,double*,double*,double&,double&,double&,double&,double*,double*,double*,
    double*,double*,double*,double*,double&,double&,double*,double*,char*,int&,int&,int&,int&,
    double*,int&,double*,double*,double&,double*,double*,double*,int&,int&,int&,int&,int*,int&);
static double call(double* SV, double* PROPS, const double* DS, double dt, double* S, double* D)
{
    double dummy = 0, pn = 1e30; int i3 = 3, i6 = 6, np = 19, nst = 70, one = 1, jstep[4] = {1,0,0,0};
    char nm[8] = "CP"; double DROT[9] = {1,0,0,0,1,0,0,0,1}, time[2] = {0,0};
    double ds[6]; std::memcpy(ds, DS, sizeof(ds));
    umat_(S, SV, D, dummy, dummy, dummy, dummy, nullptr, nullptr, nullptr, nullptr, ds, time, &dt,
          dummy, dummy, nullptr, nullptr, nm, i3, i3, i6, nst, PROPS, np, nullptr, DROT, pn,
          nullptr, nullptr, nullptr, one, one, one, one, jstep, one);
    return pn;
}
int main(int argc, char** argv)
{
    double de = argc > 1 ? atof(argv[1]) : 5e-4;
    bool rnd = argc > 2 && argv[2][0] == 'r';
    double PROPS[19] = {0, 106750, 60410, 28340, 2, 2, 1, 6, 0.588235, 50, 1, 0, 0, 0, 1, 2, 115, 34, 2.25};
    if (argc > 3) PROPS[5] = atof(argv[3]);
    if (rnd) { PROPS[10] = -0.12523701; PROPS[11] = 0.09397545; PROPS[12] = 0.59307825; PROPS[13] = 0.78977370; }
    std::array<double,70> sv; sv.fill(0);
    double base[6] = {-0.45 * de, -0.45 * de, de, 0.06 * de, -0.04 * de, 0.08 * de};
    double dt = 0.005, S[6], D[36];
    double eps_checks[] = {0.005, 0.02, 0.05, 0.1};
    int nc = 0;
    int ninc = (int)std::lround(0.1 / de);
    double worst = 0;
    for (int inc = 1; inc <= ninc; inc++) {
        double eps = inc * de;
        if (nc < 4 && eps >= eps_checks[nc] - 1e-12) {
            nc++;
            std::array<double,70> a = sv;
            call(a.data(), PROPS, base, dt, S, D);
            double Dn = 0; for (int k = 0; k < 36; k++) Dn = std::max(Dn, std::fabs(D[k]));
            double emax = 0; int wi = 0, wj = 0; double FD[36];
            for (int j = 0; j < 6; j++) {
                double h = 1e-7 * de / 5e-4, dp[6], dm[6], Sp[6], Sm[6], Dd[36];
                std::memcpy(dp, base, 48); std::memcpy(dm, base, 48); dp[j] += h; dm[j] -= h;
                std::array<double,70> b1 = sv, b2 = sv;
                call(b1.data(), PROPS, dp, dt, Sp, Dd); call(b2.data(), PROPS, dm, dt, Sm, Dd);
                for (int i = 0; i < 6; i++) {
                    FD[i + 6 * j] = (Sp[i] - Sm[i]) / (2 * h);
                    double e = std::fabs((Sp[i] - Sm[i]) / (2 * h) - D[i + 6 * j]);
                    if (e > emax) { emax = e; wi = i; wj = j; }
                }
            }
            Eigen::Map<Eigen::Matrix<double,6,6>> Ma(D), Mf(FD);
            double Et_a = 1.0 / Ma.inverse()(2, 2), Et_f = 1.0 / Mf.inverse()(2, 2);
            double relEt = std::fabs(Et_a - Et_f) / std::fabs(Et_f);
            if (relEt > worst) worst = relEt;
            printf("  eps=%.3f  uniaxial E_t: analytic %9.1f  FD %9.1f  rel.err %.2e |", eps, Et_a, Et_f, relEt);
            printf("  eps=%.3f  |DDSDDE|max=%8.0f  max|DDSDDE-FD|=%9.2f MPa  rel=%.2e  at(%d,%d)\n",
                   eps, Dn, emax, emax / Dn, wi, wj);
        }
        call(sv.data(), PROPS, base, dt, S, D);
    }
    printf("  worst E_t rel.err = %.2e\n", worst);
    return 0;
}
