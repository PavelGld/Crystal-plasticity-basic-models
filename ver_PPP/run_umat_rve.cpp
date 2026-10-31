// RVE-level check: drive the real UMAT entry (cp_umat.cpp) for 500 grains and
// compare the macro stress with the standalone RVE (run_traj.exe).
// usage: run_umat_rve <traj> <nsteps> <quat-file> <out> [cfg]
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

static const int NG = 500;
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

int main(int argc, char** argv)
{
    if (argc < 5) { std::fprintf(stderr, "usage: run_umat_rve <traj> <nsteps> <quat> <out>\n"); return 2; }
    std::string traj = argv[1];
    int nsteps = std::atoi(argv[2]);

    MaterialParams p("ver_PPP/grain_ref_g1.txt");

    std::vector<GrainState> gs(NG);
    std::ifstream qf(argv[3]);
    std::string line;
    std::getline(qf, line);
    for (auto& g : gs) std::memset(&g, 0, sizeof(g));
    while (std::getline(qf, line)) {
        if (line.empty()) continue;
        int id; double w, x, y, z;
        if (std::sscanf(line.c_str(), "%d %lf %lf %lf %lf", &id, &w, &x, &y, &z) != 5) continue;
        if (id < 0 || id >= NG) continue;
        Eigen::Quaterniond q(w, x, y, z); q.normalize();
        fill_props(gs[id].PROPS, p, q);
    }

    const double dt = 5e-3;
    const double ds = 1e-3;
    Eigen::Matrix3d Lt = Eigen::Matrix3d::Zero();
    Lt(0,0) = ds; Lt(1,1) = -0.5*ds; Lt(2,2) = -0.5*ds;
    Eigen::Matrix3d Ls = Eigen::Matrix3d::Zero();
    Ls(0,1) = std::sqrt(3.0)*ds;
    Eigen::Matrix3d Lp = Eigen::Matrix3d::Zero();
    Lp(0,0) = (std::sqrt(3.0)/2.0)*ds; Lp(1,1) = -(std::sqrt(3.0)/2.0)*ds;

    FILE* f = std::fopen(argv[4], "w");
    std::fprintf(f, "E\tS11\tS22\tS33\tS12\tS13\tS23\ttr\n");

    // Warm-up: build the shared MaterialParams serially (avoids a first-call race
    // between OpenMP threads), then reset grain 0 so its real first call re-inits.
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

    double E = 0.0;
    for (int s = 1; s <= nsteps; s++) {
        Eigen::Matrix3d L;
        if (traj == "tension" || (traj == "complex" && s <= nsteps/2)) L = Lt;
        else if (traj == "compression") L = -Lt;
        else if (traj == "simple_shear" || traj == "complex") L = Ls;
        else if (traj == "planar") L = Lp;
        else { std::fprintf(stderr, "unknown traj\n"); return 2; }

        Eigen::Matrix3d Dm = 0.5 * (L + L.transpose());

        // Abaqus passes DROT = rotation increment of the co-rotational basis.
        // For a prescribed L this is the polar rotation of exp(L*dt).
        Eigen::Matrix3d Finc = (L * dt).exp();
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(Finc, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d Rinc = svd.matrixU() * svd.matrixV().transpose();
        double DROT[9];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) DROT[i + j * 3] = Rinc(i, j);

        // DSTRAN is the strain increment in the NEW (co-rotated) basis:
        //   DSTRAN = DROT * D_lab * DROT^T * dt
        Eigen::Matrix3d Dinc = Rinc * Dm * Rinc.transpose();
        Vector6d dvp; tensorToVoigtInPlace(Dinc, dvp);       // project Voigt (x2 shear)
        Vector6d dva = A_to_V() * dvp;                        // Abaqus order, engineering shear (no /2)
        double DSTRAN[6];
        for (int i = 0; i < 6; i++) DSTRAN[i] = dva(i) * dt;

        double DDSDDE[36] = {0};
        double DSTRANd[6]; for (int i = 0; i < 6; i++) DSTRANd[i] = DSTRAN[i];
        double TIME[2] = {0.0, s * dt};
        double DDSDDT[6] = {0}, DRPLDE[6] = {0}, DRPLDT[1] = {0};
        double STRAN[6] = {0,0,0,0,0,0};
        double TEMP = 0, DTEMP = 0, PREDEF[1] = {0}, DPRED[1] = {0};
        char CMNAME[80] = "CPFCC";
        int NDI = 3, NSHR = 3, NTENS = 6, NSTATV = NV, NPROPS = 19;
        double COORDS[3] = {0,0,0};
        double CELENT = 1.0;
        double DFGRD0[9] = {1,0,0, 0,1,0, 0,0,1}, DFGRD1[9] = {1,0,0, 0,1,0, 0,0,1};
        int LAYER = 1, KSPT = 1, KINC = 0;
        int JSTEP[4] = {1, 1, 1, 1};

        static std::vector<Eigen::Matrix3d> Sgrain(NG);
        #pragma omp parallel for schedule(dynamic)
        for (int gi = 0; gi < NG; gi++) {
            GrainState& g = gs[gi];
            double DTIME = dt;
            double DDSDDE[36] = {0};
            double SSE = 0, SPD = 0, SCD = 0, RPL = 0;
            double PNEWDT = 1e10;
            int NOEL = gi + 1, NPT = 1;
            umat_(g.STRESS, g.STATEV, DDSDDE, SSE, SPD, SCD, RPL,
                  DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRANd, TIME, &DTIME,
                  TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
                  g.PROPS, NPROPS, COORDS, DROT, PNEWDT, &CELENT, DFGRD0, DFGRD1,
                  NOEL, NPT, LAYER, KSPT, JSTEP, KINC);
            Eigen::Matrix3d S;
            S(0,0)=g.STRESS[0]; S(1,1)=g.STRESS[1]; S(2,2)=g.STRESS[2];
            S(0,1)=S(1,0)=g.STRESS[3]; S(0,2)=S(2,0)=g.STRESS[4]; S(1,2)=S(2,1)=g.STRESS[5];
            Sgrain[gi] = S;
        }
        Eigen::Matrix3d Savg = Eigen::Matrix3d::Zero();
        for (int gi = 0; gi < NG; gi++) Savg += Sgrain[gi];
        Savg /= (double)NG;
        E += 5e-6;
        if (s % 2000 == 0) {
            std::fprintf(f, "%.6f\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
                         E, Savg(0,0)/1e9, Savg(1,1)/1e9, Savg(2,2)/1e9,
                         Savg(0,1)/1e9, Savg(0,2)/1e9, Savg(1,2)/1e9,
                         (Savg(0,0)+Savg(1,1)+Savg(2,2))/1e9);
        }
    }
    std::fclose(f);
    return 0;
}
