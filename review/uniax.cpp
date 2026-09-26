// Uniaxial-stress material-point driver through the exact UMAT entry (umat_).
// Mimics one C3D8 element under uniaxial tension (tension.inp): prescribe dε33 per
// increment, solve the other 5 DSTRAN components by Newton on DDSDDE (with backtracking
// line search) so that the other 5 stress components vanish. Each equilibrium iteration
// restarts from the start-of-increment STATEV (as Abaqus does). Abaqus-like increment
// control: cutback x0.25 on failure / PNEWDT<1, growth x1.5 after two easy increments,
// capped by the requested increment.
//
// usage: uniax <dt> <deps33> <eps_end> [q: i|r] [n_power]
#include "cp_umat/cp_umat_interface.hpp"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <array>

extern "C" void umat_(double*,double*,double*,double&,double&,double&,double&,double*,double*,double*,
    double*,double*,double*,double*,double&,double&,double*,double*,char*,int&,int&,int&,int&,
    double*,int&,double*,double*,double&,double*,double*,double*,int&,int&,int&,int&,int*,int&);

static long g_calls = 0;
static double call_umat(double* STATEV, double* PROPS, double* DSTRAN, double dt,
                        double* S, double* D, int kinc)
{
    double dummy = 0, pn = 1e30; int i3 = 3, i6 = 6, np = 19, nst = 70, one = 1;
    int jstep[4] = {1, 0, 0, 0}; char nm[8] = "CP";
    double DROT[9] = {1,0,0, 0,1,0, 0,0,1};
    double time[2] = {0, 0};
    umat_(S, STATEV, D, dummy, dummy, dummy, dummy, nullptr, nullptr, nullptr, nullptr, DSTRAN, time, &dt,
          dummy, dummy, nullptr, nullptr, nm, i3, i3, i6, nst, PROPS, np, nullptr, DROT, pn,
          nullptr, nullptr, nullptr, one, one, one, one, jstep, kinc);
    g_calls++;
    return pn;
}

int main(int argc, char** argv)
{
    double dt0   = argc > 1 ? atof(argv[1]) : 0.005;
    double de0   = argc > 2 ? atof(argv[2]) : 5e-4;
    double e_end = argc > 3 ? atof(argv[3]) : 0.1;
    bool rnd     = argc > 4 && argv[4][0] == 'r';
    double npow  = argc > 5 ? atof(argv[5]) : 50.0;
    const double rate = de0 / dt0;

    double PROPS[19] = {0, 106750, 60410, 28340, 2, 2, 1, 6, 0.588235, npow,
                        1, 0, 0, 0, 1, 2, 115, 34, 2.25};
    if (rnd) { PROPS[10] = -0.12523701; PROPS[11] = 0.09397545; PROPS[12] = 0.59307825; PROPS[13] = 0.78977370; }

    std::array<double,70> sv; sv.fill(0.0);
    const int free_idx[5] = {0, 1, 3, 4, 5};
    double dir[6] = {-0.37, -0.37, 1, 0, 0, 0};   // strain-increment "direction" for extrapolation
    double S[6], D[36];
    double eps = 0.0, de = de0;
    int inc = 0, tot_it = 0, max_it = 0, cutbacks = 0, easy = 0;
    double checkpoints[] = {0.01, 0.02, 0.05, 0.08, 0.10, 0.2, 0.3};
    int cp = 0;
    auto rmax_of = [&](const double* s) { double r = 0; for (int k : free_idx) r = std::max(r, std::fabs(s[k])); return r; };
    auto t0 = std::chrono::high_resolution_clock::now();
    while (eps < e_end - 1e-12) {
        if (inc > 100000) { printf("too many increments\n"); return 1; }
        double de_try = std::min(de, e_end - eps);
        double dt = de_try / rate;
        double DS[6]; for (int k = 0; k < 6; k++) DS[k] = dir[k] * de_try;
        std::array<double,70> sv_new;
        bool ok = false; int it;
        double pn = call_umat((sv_new = sv).data(), PROPS, DS, dt, S, D, inc + 1);
        double r = (pn < 1.0) ? 1e300 : rmax_of(S);
        for (it = 1; it <= 16 && pn >= 1.0; it++) {
            if (r < 1e-7 * std::max(1.0, std::fabs(S[2]))) { ok = true; break; }
            Eigen::Matrix<double,5,5> K; Eigen::Matrix<double,5,1> rv;
            for (int a = 0; a < 5; a++) {
                rv(a) = S[free_idx[a]];
                for (int b = 0; b < 5; b++) K(a, b) = D[free_idx[a] + 6 * free_idx[b]];
            }
            Eigen::Matrix<double,5,1> dx = K.partialPivLu().solve(-rv);
            double lam = 1.0; bool acc = false;
            for (int bt = 0; bt < 8; bt++) {
                double DT[6]; std::memcpy(DT, DS, sizeof(DT));
                for (int a = 0; a < 5; a++) DT[free_idx[a]] += lam * dx(a);
                double St[6], Dt[36];
                double pnt = call_umat((sv_new = sv).data(), PROPS, DT, dt, St, Dt, inc + 1);
                double rt = (pnt < 1.0) ? 1e300 : rmax_of(St);
                if (rt < r) { std::memcpy(DS, DT, sizeof(DS)); std::memcpy(S, St, sizeof(S)); std::memcpy(D, Dt, sizeof(D)); r = rt; acc = true; break; }
                lam *= 0.5;
            }
            if (!acc) break;
        }
        if (!ok) { cutbacks++; de *= 0.25; easy = 0; if (de < 1e-9) { printf("cutback limit\n"); return 1; } continue; }
        inc++; tot_it += it; if (it > max_it) max_it = it;
        sv = sv_new; eps += de_try;
        for (int k = 0; k < 6; k++) dir[k] = DS[k] / de_try;                 // extrapolation
        if (it <= 4) { if (++easy >= 2) { de = std::min(de * 1.5, de0); } } else easy = 0;
        if (cp < 7 && eps >= checkpoints[cp] - 1e-12 && eps <= e_end + 1e-12) {
            double tcm = 0; for (int i = 0; i < 24; i++) tcm += sv[10 + i]; tcm /= 24;
            double sh = 0;  for (int i = 0; i < 24; i++) sh += sv[34 + i];
            printf("  eps=%.3f  S33=%9.4f MPa   <tc>=%8.4f   sum_gamma=%.4f\n", eps, S[2], tcm, sh);
            cp++;
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
    printf("  incs=%d cutbacks=%d  equil.iters avg=%.2f max=%d  UMAT calls=%ld  %.2f us/call  total %.1f ms\n",
           inc, cutbacks, (double)tot_it / inc, max_it, g_calls, us / g_calls, us / 1000);
    return 0;
}
