// Mini-FE: N×N×N C3D8 cube, one random orientation per element, global Newton with the
// exact UMAT entry (umat_). Abaqus-like controls: residual tol 0.5% of reaction force,
// max 16 iterations, backtracking line search, cutback ×0.25 on failure / PNEWDT<1,
// growth ×1.5 after two increments with ≤ I_G iterations (I_G = 4), dt ≤ dtmax.
// usage: fe <N> <eps_end> <dt0> <dtmax> [n_power]
#include "cp_umat/cp_umat_interface.hpp"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <chrono>

extern "C" void umat_(double*,double*,double*,double&,double&,double&,double&,double*,double*,double*,
    double*,double*,double*,double*,double&,double&,double*,double*,char*,int&,int&,int&,int&,
    double*,int&,double*,double*,double&,double*,double*,double*,int&,int&,int&,int&,int*,int&);

static long g_calls = 0;
static double call_umat(double* SV, double* PROPS, double* DS, double dt, double* S, double* D, int noel, int npt)
{
    double dummy = 0, pn = 1e30; int i3 = 3, i6 = 6, np = 19, nst = 70, one = 1, js[4] = {1,0,0,0};
    char nm[8] = "CP"; double DROT[9] = {1,0,0,0,1,0,0,0,1}, tm[2] = {0,0};
    umat_(S, SV, D, dummy, dummy, dummy, dummy, nullptr, nullptr, nullptr, nullptr, DS, tm, &dt,
          dummy, dummy, nullptr, nullptr, nm, i3, i3, i6, nst, PROPS, np, nullptr, DROT, pn,
          nullptr, nullptr, nullptr, noel, npt, one, one, js, one);
    g_calls++;
    return pn;
}

int main(int argc, char** argv)
{
    int N = argc > 1 ? atoi(argv[1]) : 3;
    double e_end = argc > 2 ? atof(argv[2]) : 0.1;
    double dt0 = argc > 3 ? atof(argv[3]) : 5e-4;
    double dtmax = argc > 4 ? atof(argv[4]) : 5e-3;
    double npow = argc > 5 ? atof(argv[5]) : 50;
    const int NN = N + 1, nnode = NN * NN * NN, ndof = 3 * nnode, nel = N * N * N, nip = 8;
    const double h = 1.0 / N;
    auto nid = [&](int i, int j, int k) { return (k * NN + j) * NN + i; };
    std::vector<std::array<int,8>> conn(nel);
    { int e = 0; for (int k = 0; k < N; k++) for (int j = 0; j < N; j++) for (int i = 0; i < N; i++, e++)
        conn[e] = {nid(i,j,k), nid(i+1,j,k), nid(i+1,j+1,k), nid(i,j+1,k), nid(i,j,k+1), nid(i+1,j,k+1), nid(i+1,j+1,k+1), nid(i,j+1,k+1)}; }
    // orientations
    std::mt19937 rng(42); std::normal_distribution<double> gauss(0, 1);
    std::vector<std::array<double,19>> PROPS(nel);
    for (int e = 0; e < nel; e++) {
        double q[4]; double nq = 0; for (int c = 0; c < 4; c++) { q[c] = gauss(rng); nq += q[c] * q[c]; }
        nq = std::sqrt(nq);
        PROPS[e] = {0, 106750, 60410, 28340, 2, 2, 1, 6, 0.588235, npow, q[0]/nq, q[1]/nq, q[2]/nq, q[3]/nq, 1, 2, 115, 34, 2.25};
    }
    // BCs (tension.inp style): bottom U3=0, top U3=e*L, node(0,0,0) U1=U2=0, node(N,0,0) U2=0
    std::vector<int> fx; std::vector<double> amp;
    for (int j = 0; j < NN; j++) for (int i = 0; i < NN; i++) { fx.push_back(3*nid(i,j,0)+2); amp.push_back(0); }
    for (int j = 0; j < NN; j++) for (int i = 0; i < NN; i++) { fx.push_back(3*nid(i,j,N)+2); amp.push_back(1.0); }
    fx.push_back(3*nid(0,0,0)); amp.push_back(0); fx.push_back(3*nid(0,0,0)+1); amp.push_back(0);
    fx.push_back(3*nid(N,0,0)+1); amp.push_back(0);
    std::vector<char> isfx(ndof, 0); for (int d : fx) isfx[d] = 1;
    // B matrices (regular cube elements)
    double gp = 1.0 / std::sqrt(3.0), lc[8][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
    double xi[8][3]; { int a = 0; for (int c = 0; c < 2; c++) for (int b = 0; b < 2; b++) for (int aa = 0; aa < 2; aa++, a++) { xi[a][0] = aa ? gp : -gp; xi[a][1] = b ? gp : -gp; xi[a][2] = c ? gp : -gp; } }
    double detJ = std::pow(h / 2, 3);
    std::vector<Eigen::Matrix<double,6,24>> B(nip);
    for (int g = 0; g < nip; g++) {
        B[g].setZero();
        for (int a = 0; a < 8; a++) {
            double dN[3];
            for (int d = 0; d < 3; d++) { double p = 0.125 * lc[a][d]; for (int o = 0; o < 3; o++) if (o != d) p *= (1 + lc[a][o] * xi[g][o]); dN[d] = p / (h / 2); }
            int c = 3 * a;
            B[g](0,c) = dN[0]; B[g](1,c+1) = dN[1]; B[g](2,c+2) = dN[2];
            B[g](3,c) = dN[1]; B[g](3,c+1) = dN[0]; B[g](4,c) = dN[2]; B[g](4,c+2) = dN[0]; B[g](5,c+1) = dN[2]; B[g](5,c+2) = dN[1];
        }
    }
    std::vector<std::array<double,70>> svc(nel * nip), svp(nel * nip);
    for (auto& s : svc) s.fill(0);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(ndof), du_prev = Eigen::VectorXd::Zero(ndof);
    double t = 0, dt = dt0, dt_prev = 0;
    int inc = 0, cut = 0, tot_it = 0, easy = 0;
    // assemble residual / stiffness at trial displacement ut
    auto assemble = [&](const Eigen::VectorXd& ut, double dtime, Eigen::VectorXd& R, Eigen::MatrixXd* K, double& pnmin) {
        R.setZero(ndof); if (K) K->setZero(ndof, ndof); pnmin = 1e30;
        for (int e = 0; e < nel; e++) {
            Eigen::Matrix<double,24,1> ue;
            for (int a = 0; a < 8; a++) for (int d = 0; d < 3; d++) ue(3*a+d) = ut(3*conn[e][a]+d) - u(3*conn[e][a]+d);
            for (int g = 0; g < nip; g++) {
                Eigen::Matrix<double,6,1> de = B[g] * ue;
                double DS[6], S[6], D[36]; for (int k = 0; k < 6; k++) DS[k] = de(k);
                int ip = e * nip + g; svp[ip] = svc[ip];
                double pn = call_umat(svp[ip].data(), PROPS[e].data(), DS, dtime, S, D, e + 1, g + 1);
                pnmin = std::min(pnmin, pn);
                Eigen::Map<Eigen::Matrix<double,6,1>> Sv(S); Eigen::Map<Eigen::Matrix<double,6,6>> Dm(D);
                Eigen::Matrix<double,24,1> f = B[g].transpose() * Sv * detJ;
                for (int a = 0; a < 24; a++) R(3*conn[e][a/3] + a%3) += f(a);
                if (K) {
                    Eigen::Matrix<double,24,24> ke = B[g].transpose() * Dm * B[g] * detJ;
                    for (int a = 0; a < 24; a++) for (int b = 0; b < 24; b++) (*K)(3*conn[e][a/3]+a%3, 3*conn[e][b/3]+b%3) += ke(a, b);
                }
            }
        }
    };
    auto rnorms = [&](const Eigen::VectorXd& R, double& rmax, double& fref) {
        rmax = 0; fref = 0; for (int d = 0; d < ndof; d++) { if (isfx[d]) fref += std::fabs(R(d)); else rmax = std::max(rmax, std::fabs(R(d))); }
        fref /= fx.size();
    };
    auto t0 = std::chrono::high_resolution_clock::now();
    while (t < e_end - 1e-12) {
        if (inc > 5000) { printf("too many increments\n"); return 1; }
        double tn = std::min(t + dt, e_end), dtime = tn - t;
        Eigen::VectorXd ut = u;
        if (dt_prev > 0) ut += du_prev * (dtime / dt_prev);        // extrapolation
        for (size_t c = 0; c < fx.size(); c++) ut(fx[c]) = amp[c] * tn;
        Eigen::VectorXd R; Eigen::MatrixXd K; double pn, rmax, fref;
        bool ok = false; int it;
        assemble(ut, dtime, R, &K, pn);
        for (it = 1; it <= 16; it++) {
            if (pn < 1.0) break;
            rnorms(R, rmax, fref);
            if (rmax <= 0.005 * std::max(fref, 1e-10)) { ok = true; break; }
            for (int d = 0; d < ndof; d++) if (isfx[d]) { K.row(d).setZero(); K.col(d).setZero(); K(d, d) = 1; R(d) = 0; }
            Eigen::VectorXd dx = K.partialPivLu().solve(-R);
            double r0 = R.norm(), lam = 1.0; bool acc = false;
            for (int bt = 0; bt < 5; bt++) {
                Eigen::VectorXd uu = ut + lam * dx, Rt; Eigen::MatrixXd Kt; double pnt;
                assemble(uu, dtime, Rt, &Kt, pnt);
                Eigen::VectorXd Rf = Rt; for (int d = 0; d < ndof; d++) if (isfx[d]) Rf(d) = 0;
                if (pnt >= 1.0 && (Rf.norm() <= 0.95 * r0 || bt == 4)) { ut = uu; R = Rt; K = Kt; pn = pnt; acc = true; break; }
                lam *= 0.5;
            }
            if (!acc) break;
        }
        if (!ok) { cut++; dt *= 0.25; easy = 0; if (dt < 1e-9) { printf("cutback limit at t=%g\n", t); return 1; } continue; }
        inc++; tot_it += it;
        du_prev = ut - u; dt_prev = dtime;
        u = ut; t = tn; svc = svp;
        if (it <= 4) { if (++easy >= 2) { dt = std::min(dt * 1.5, dtmax); easy = 0; } } else easy = 0;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    // mean S33 per element from last converged state: probe with zero increment
    double smin = 1e30, smax = -1e30, smean = 0;
    for (int e = 0; e < nel; e++) {
        double s33 = 0;
        for (int g = 0; g < nip; g++) { auto sv = svc[e*nip+g]; double DS[6] = {0}, S[6], D[36]; call_umat(sv.data(), PROPS[e].data(), DS, 0.0, S, D, e+1, g+1); s33 += S[2] / nip; }
        smin = std::min(smin, s33); smax = std::max(smax, s33); smean += s33 / nel;
    }
    printf("  t=%.3f incs=%d cutbacks=%d global iters=%d (avg %.2f/inc)  UMAT calls=%ld  %.1f us/call  wall %.2f s | S33 mean %.3f [%.2f..%.2f]\n",
           t, inc, cut, tot_it, (double)tot_it / inc, g_calls, 1000 * ms / g_calls, ms / 1000, smean, smin, smax);
    return 0;
}
