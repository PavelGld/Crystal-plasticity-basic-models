// FD consistency test: DDSDDE vs d(STRESS)/d(DSTRAN) through the exact UMAT path.
#include "cp_umat/cp_umat_interface.hpp"
#include <cstdio>
#include <cstring>
#include <cmath>

static CPUMAT g_cp;
static bool g_ready = false;

// replicate umat_() flow for one integration point
void umat_call(double* STATEV, const double* PROPS, const double* DSTRAN, double DTIME,
               double* STRESS, double* DDSDDE, double& PNEWDT)
{
    if (!g_ready) { int np = ((int)PROPS[4]==2)?19:((int)PROPS[4]==1)?16:14; g_cp.init_params(PROPS, np); g_ready = true; }
    Grain* g = g_cp.grain();
    bool first = (std::fabs(STATEV[69] - 123456.0) > 0.5);
    if (first) {
        Eigen::Quaterniond q(PROPS[10], PROPS[11], PROPS[12], PROPS[13]);
        q.normalize();
        g->set_orientation(q);
        g->reset_state();
        CPUMAT::save_state(g, STATEV);
        STATEV[69] = 123456.0;
    } else {
        CPUMAT::load_state(g, STATEV);
    }
    g_cp.step(g, DSTRAN, DTIME, STRESS, DDSDDE, PNEWDT);
    CPUMAT::save_state(g, STATEV);
}

int main(int argc, char** argv)
{
    // PROPS from tension_8el_random.inp, quaternion of CP1
    double PROPS[19] = {0, 106750, 60410, 28340, 2, 2, 1, 6,
                        0.588235, 50,
                        -0.19326972, -0.23191724, -0.14930902, 0.94157741,
                        1, 2, 115, 34, 2.25};
    if (argc > 1 && argv[1][0]=='i') { PROPS[10]=1; PROPS[11]=PROPS[12]=PROPS[13]=0; printf("== identity orientation =="); }
    else printf("== random orientation ==");
    for (int a=2; a<argc; a++) { if (argv[a][0]>='0'&&argv[a][0]<='9') continue;
        if (argv[a][0]=='R') { PROPS[5]=0; printf(" rot_mode=0"); }          // no lattice rotation
        if (argv[a][0]=='H') { PROPS[4]=0; printf(" hard_mode=0"); }         // no hardening (NPROPS mismatch ok: init reads first 10)
        if (argv[a][0]=='C') { printf(" (clamp probe)"); }
    }
    printf("\n");

    double STATEV[70]; std::memset(STATEV, 0, sizeof(STATEV));
    double DTIME = 0.005;
    // strain path: uniaxial z with elastic-ish lateral contraction + small shears (generic direction)
    double base[6] = {-1.8e-4, -1.8e-4, 5.0e-4, 2e-5, -1e-5, 3e-5};

    double STRESS[6], DDSDDE[36]; double PNEWDT;

    for (int inc = 1; inc <= 60; inc++) {
        // FD check at selected increments, from the same start-of-increment state
        if (inc==1 || inc==10 || inc==30 || inc==60) {
            double sv0[70]; std::memcpy(sv0, STATEV, sizeof(sv0));
            // unperturbed call to get DDSDDE
            double sv[70]; std::memcpy(sv, sv0, sizeof(sv));
            umat_call(sv, PROPS, base, DTIME, STRESS, DDSDDE, PNEWDT);
            double maxrel = 0, maxabs = 0; double Dnorm = 0;
            for (int k=0;k<36;k++) Dnorm = std::max(Dnorm, std::fabs(DDSDDE[k]));
            for (int j = 0; j < 6; j++) {
                double h = 1e-8; { double t_ = atof(argv[argc-1]); if (t_>0 && t_<1e-3) h = t_; }
                double dp[6], dm[6]; std::memcpy(dp, base, sizeof(dp)); std::memcpy(dm, base, sizeof(dm));
                dp[j] += h; dm[j] -= h;
                double Sp[6], Sm[6], Dd[36]; double pn;
                double svp[70]; std::memcpy(svp, sv0, sizeof(svp));
                umat_call(svp, PROPS, dp, DTIME, Sp, Dd, pn);
                double svm[70]; std::memcpy(svm, sv0, sizeof(svm));
                umat_call(svm, PROPS, dm, DTIME, Sm, Dd, pn);
                for (int i = 0; i < 6; i++) {
                    double fd = (Sp[i]-Sm[i])/(2*h);
                    double an = DDSDDE[i + j*6];
                    double err = std::fabs(fd - an);
                    maxabs = std::max(maxabs, err);
                    maxrel = std::max(maxrel, err / Dnorm);
                }
            }
            printf("inc %3d | max|DDSDDE| %10.1f | FD err abs %12.4f  rel %10.3e | S33 %8.4f | PNEWDT %.2f\n",
                   inc, Dnorm, maxabs, maxrel, STRESS[2], PNEWDT);
            std::memcpy(STATEV, sv, sizeof(STATEV)); // accept unperturbed step
        } else {
            umat_call(STATEV, PROPS, base, DTIME, STRESS, DDSDDE, PNEWDT);
        }
    }
    printf("final S33 = %.4f MPa, sim strain e33 = %.4f\n", STRESS[2], 60*base[2]);
    return 0;
}
