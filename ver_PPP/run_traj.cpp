// Run one verification trajectory with the standalone RVE (== UMAT after the
// Ms_v fix) and dump macro stress on the E = 0.01 grid.
// usage: run_traj <traj> <nsteps> <quat-file> <out> [conj]
#include "cp_model/RVE.hpp"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

int main(int argc, char** argv)
{
    if (argc < 5) { std::fprintf(stderr, "usage: run_traj <traj> <nsteps> <quat> <out> [cfg] [conj]\n"); return 2; }
    std::string traj = argv[1];
    int nsteps = std::atoi(argv[2]);
    const char* cfg = (argc > 5) ? argv[5] : "ver_PPP/rve_ref.txt";
    bool conj = (argc > 6 && std::strcmp(argv[6], "conj") == 0);

    RVE rve(cfg);

    std::ifstream qf(argv[3]);
    if (!qf) { std::fprintf(stderr, "cannot open %s\n", argv[3]); return 2; }
    std::string line;
    std::getline(qf, line);   // header
    int n = 0;
    while (std::getline(qf, line)) {
        if (line.empty()) continue;
        int id; double w, x, y, z;
        if (std::sscanf(line.c_str(), "%d %lf %lf %lf %lf", &id, &w, &x, &y, &z) != 5) continue;
        Eigen::Quaterniond q(w, x, y, z);
        q.normalize();
        if (conj) q = q.conjugate();
        if (id < (int)rve.num_grains()) rve.get_grain(id).set_orientation(q);
        n++;
    }
    std::fprintf(stderr, "orientations: %d grains=%zu\n", n, rve.num_grains());

    const double dt = 5e-3;
    const double ds = 1e-3;
    Eigen::Matrix3d Lt = Eigen::Matrix3d::Zero();
    Lt(0,0) = ds; Lt(1,1) = -0.5*ds; Lt(2,2) = -0.5*ds;
    Eigen::Matrix3d Ls = Eigen::Matrix3d::Zero();
    Ls(0,1) = std::sqrt(3.0)*ds;
    Eigen::Matrix3d Lp = Eigen::Matrix3d::Zero();
    Lp(0,0) = (std::sqrt(3.0)/2.0)*ds; Lp(1,1) = -(std::sqrt(3.0)/2.0)*ds;

    FILE* f = std::fopen(argv[4], "w");
    if (!f) { std::fprintf(stderr, "cannot write %s\n", argv[4]); return 2; }
    std::fprintf(f, "E\tS11\tS22\tS33\tS12\tS13\tS23\ttr\n");

    double E = 0.0;
    for (int s = 1; s <= nsteps; s++) {
        Eigen::Matrix3d L;
        if (traj == "tension" || (traj == "complex" && s <= nsteps/2)) L = Lt;
        else if (traj == "compression") L = -Lt;
        else if (traj == "simple_shear" || traj == "complex") L = Ls;
        else if (traj == "planar") L = Lp;
        else { std::fprintf(stderr, "unknown traj\n"); return 2; }

        rve.Simulate_Step(dt, L);
        E += 5e-6;
        if (s % 2000 == 0) {
            Eigen::Matrix3d S = rve.get_stress();
            std::fprintf(f, "%.6f\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
                         E, S(0,0)/1e9, S(1,1)/1e9, S(2,2)/1e9,
                         S(0,1)/1e9, S(0,2)/1e9, S(1,2)/1e9,
                         (S(0,0)+S(1,1)+S(2,2))/1e9);
            std::fprintf(stderr, "E=%.2f S11=%.6f\n", E, S(0,0)/1e9);
        }
    }
    std::fclose(f);

    std::string qout = std::string(argv[4]) + ".quats";
    FILE* qfout = std::fopen(qout.c_str(), "w");
    if (qfout) {
        std::fprintf(qfout, "grain_id\tw\tx\ty\tz\n");
        for (size_t i = 0; i < rve.num_grains(); i++) {
            Eigen::Quaterniond q = rve.get_grain(i).get_orientation();
            q.normalize();
            std::fprintf(qfout, "%d\t%.12e\t%.12e\t%.12e\t%.12e\n",
                         (int)i, q.w(), q.x(), q.y(), q.z());
        }
        std::fclose(qfout);
    }
    return 0;
}
