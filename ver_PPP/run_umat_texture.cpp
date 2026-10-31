// UMAT-driven texture validation: 1000 grains, 4 loading modes, checkpoints
// at E = 0.3 / 0.5 / 1.0, quaternions written in the same format as
// rotation_validation_test_results/orientations_*.txt (id qw qx qy qz).
// usage: run_umat_texture <mode> [dt]
#include "cp_model/Grain.hpp"
#include "cp_umat/cp_umat_types.hpp"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

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

static const int NV = 70;

struct GrainState {
    double PROPS[19];
    double STATEV[NV];
    double STRESS[6];
};

static void fill_props(double* P, const MaterialParams& p, const Eigen::Quaterniond& q)
{
    P[0] = 0.0;
    P[1] = p.P_v(0, 0); P[2] = p.P_v(0, 1); P[3] = 0.5 * p.P_v(3, 3);
    P[4] = (double)p.hard_mode; P[5] = (double)p.rot_mode; P[6] = (double)p.integ_mode;
    P[7] = p.tc0; P[8] = p.strain_rate_infl; P[9] = p.n_power;
    P[10] = q.w(); P[11] = q.x(); P[12] = q.y(); P[13] = q.z();
    P[14] = p.A; P[15] = p.q; P[16] = p.h0; P[17] = p.ts; P[18] = p.alfa;
}

static void export_quats(const std::vector<GrainState>& gs, const std::string& path)
{
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return;
    std::fprintf(f, "# id qw qx qy qz\n");
    for (size_t i = 0; i < gs.size(); i++) {
        const double* s = gs[i].STATEV;
        double n = std::sqrt(s[6]*s[6] + s[7]*s[7] + s[8]*s[8] + s[9]*s[9]);
        if (n < 1e-30) n = 1.0;
        std::fprintf(f, "%d %.12e %.12e %.12e %.12e\n",
                     (int)i, s[6]/n, s[7]/n, s[8]/n, s[9]/n);
    }
    std::fclose(f);
}

int main(int argc, char** argv)
{
    if (argc < 2) { std::fprintf(stderr, "usage: run_umat_texture <mode> [dt]\n"); return 2; }
    std::string mode = argv[1];
    const double dt = (argc > 2) ? std::atof(argv[2]) : 1e-4;

    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    if (mode == "simple_shear") {
        L(0, 1) = 0.1;
    } else if (mode == "uniaxial_tension") {
        L(0,0) = -0.05; L(1,1) = -0.05; L(2,2) = 0.1;
    } else if (mode == "uniaxial_compression") {
        L(0,0) = 0.05; L(1,1) = 0.05; L(2,2) = -0.1;
    } else if (mode == "plane_strain") {
        L(1,1) = 0.1; L(2,2) = -0.1;
    } else {
        std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
        return 2;
    }

    MaterialParams p("ver_PPP/grain_texture.txt");

    std::ifstream qf("ver_PPP/texture_initial.txt");
    if (!qf) { std::fprintf(stderr, "cannot open texture_initial.txt\n"); return 2; }
    std::string line;
    std::getline(qf, line);
    std::vector<GrainState> gs;
    while (std::getline(qf, line)) {
        if (line.empty() || line[0] == '#') continue;
        int id; double w, x, y, z;
        if (std::sscanf(line.c_str(), "%d %lf %lf %lf %lf", &id, &w, &x, &y, &z) != 5) continue;
        GrainState g;
        std::memset(&g, 0, sizeof(g));
        Eigen::Quaterniond q(w, x, y, z); q.normalize();
        fill_props(g.PROPS, p, q);
        gs.push_back(g);
    }
    const int NG = (int)gs.size();
    std::fprintf(stderr, "%s: %d grains, dt=%g\n", mode.c_str(), NG, dt);

    // rotation increment and strain increment in the co-rotated basis
    Eigen::Matrix3d Finc = (L * dt).exp();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(Finc, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d Rinc = svd.matrixU() * svd.matrixV().transpose();
    double DROT[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) DROT[i + j * 3] = Rinc(i, j);

    Eigen::Matrix3d Dm = 0.5 * (L + L.transpose());
    Eigen::Matrix3d Dinc = Rinc * Dm * Rinc.transpose();
    Vector6d dvp; tensorToVoigtInPlace(Dinc, dvp);
    Vector6d dva = A_to_V() * dvp;
    double DSTRAN[6];
    for (int i = 0; i < 6; i++) DSTRAN[i] = dva(i) * dt;

    // D_eq for the equivalent strain accumulation
    Eigen::Matrix3d dev = Dm - Dm.trace() / 3.0 * Eigen::Matrix3d::Identity();
    double Deq = std::sqrt(2.0 / 3.0 * (dev * dev).trace());

    // warm-up: build shared MaterialParams serially
    {
        double DTIME0 = 0.0, DDSDDE0[36] = {0}, DSTRAN0[6] = {0,0,0,0,0,0};
        double STRAN0[6] = {0,0,0,0,0,0}, TIME0[2] = {0,0};
        double SSE0=0,SPD0=0,SCD0=0,RPL0=0, DDSDDT0[6]={0},DRPLDE0[6]={0},DRPLDT0[1]={0};
        double TEMP0=0,DTEMP0=0,PREDEF0[1]={0},DPRED0[1]={0};
        char CMN0[80]="CPFCC";
        int NDI0=3,NSHR0=3,NTENS0=6,NSTATV0=NV,NPROPS0=19;
        double COORDS0[3]={0,0,0},DROT0[9]={1,0,0,0,1,0,0,0,1};
        double PNEWDT0=1e10,CELENT0=1.0;
        double DFG0[9]={1,0,0,0,1,0,0,0,1},DFG1[9]={1,0,0,0,1,0,0,0,1};
        int NOEL0=1,NPT0=1,LAYER0=1,KSPT0=1,KINC0=0,JSTEP0[4]={1,1,1,1};
        umat_(gs[0].STRESS, gs[0].STATEV, DDSDDE0, SSE0, SPD0, SCD0, RPL0,
              DDSDDT0, DRPLDE0, DRPLDT0, STRAN0, DSTRAN0, TIME0, &DTIME0,
              TEMP0, DTEMP0, PREDEF0, DPRED0, CMN0, NDI0, NSHR0, NTENS0, NSTATV0,
              gs[0].PROPS, NPROPS0, COORDS0, DROT0, PNEWDT0, &CELENT0, DFG0, DFG1,
              NOEL0, NPT0, LAYER0, KSPT0, JSTEP0, KINC0);
        std::memset(gs[0].STATEV, 0, sizeof(gs[0].STATEV));
        std::memset(gs[0].STRESS, 0, sizeof(gs[0].STRESS));
    }

    const double checkpoints[3] = {0.3, 0.5, 1.0};
    int cp_done = 0;
    double E = 0.0;
    int step = 0;
    std::vector<Eigen::Matrix3d> Sgrain(NG);

    const double DSTRANd[6] = {DSTRAN[0],DSTRAN[1],DSTRAN[2],DSTRAN[3],DSTRAN[4],DSTRAN[5]};
    while (E < 1.0 - 1e-12 && cp_done < 3) {
        double TIME[2] = {0.0, (step + 1) * dt};
        double DDSDDT[6] = {0}, DRPLDE[6] = {0}, DRPLDT[1] = {0};
        double STRAN[6] = {0,0,0,0,0,0};
        double TEMP = 0, DTEMP = 0, PREDEF[1] = {0}, DPRED[1] = {0};
        char CMNAME[80] = "CPFCC";
        int NDI = 3, NSHR = 3, NTENS = 6, NSTATV = NV, NPROPS = 19;
        double COORDS[3] = {0,0,0};
        double CELENT = 1.0;
        double DFGRD0[9] = {1,0,0,0,1,0,0,0,1}, DFGRD1[9] = {1,0,0,0,1,0,0,0,1};
        int LAYER = 1, KSPT = 1, KINC = 0;
        int JSTEP[4] = {1,1,1,1};

        #pragma omp parallel for schedule(dynamic)
        for (int gi = 0; gi < NG; gi++) {
            GrainState& g = gs[gi];
            double DTIME = dt;
            double DDSDDE[36] = {0};
            double SSE = 0, SPD = 0, SCD = 0, RPL = 0;
            double PNEWDT = 1e10;
            int NOEL = gi + 1, NPT = 1;
            umat_(g.STRESS, g.STATEV, DDSDDE, SSE, SPD, SCD, RPL,
                  DDSDDT, DRPLDE, DRPLDT, STRAN, (double*)DSTRANd, TIME, &DTIME,
                  TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
                  g.PROPS, NPROPS, COORDS, DROT, PNEWDT, &CELENT, DFGRD0, DFGRD1,
                  NOEL, NPT, LAYER, KSPT, JSTEP, KINC);
        }

        step++;
        E += Deq * dt;

        while (cp_done < 3 && E >= checkpoints[cp_done] - 1e-12) {
            char name[256];
            std::snprintf(name, sizeof(name), "ver_PPP/texture_umat_%s_%dpct.txt",
                          mode.c_str(), (int)(checkpoints[cp_done] * 100 + 0.5));
            export_quats(gs, name);
            std::fprintf(stderr, "  E=%.4f -> %s\n", E, name);
            cp_done++;
        }
    }

    std::fprintf(stderr, "%s done: %d steps, E=%.6f\n", mode.c_str(), step, E);
    return 0;
}
