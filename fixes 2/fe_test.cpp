// Mini-FE: 2x2x2 C3D8 cube, per-element orientations, global Newton with UMAT.
// Replicates tension_8el_random.inp (BCs, loading, increments) without Abaqus.
#include "cp_umat/cp_umat_interface.hpp"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <array>
#include <cstdlib>
#include <algorithm>

static CPUMAT g_cp;
static bool g_ready = false;

void umat_call(double* STATEV, const double* PROPS, const double* DSTRAN, double DTIME,
               double* STRESS, double* DDSDDE, double& PNEWDT)
{
    if (!g_ready) { g_cp.init_params(PROPS, 19); g_ready = true; }
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
    setvbuf(stdout, NULL, _IONBF, 0);
    bool same_q = (argc > 1 && argv[1][0]=='s');   // "s" -> all elements same non-identity quaternion
    bool identity = (argc > 1 && argv[1][0]=='i');
    double t_end = (argc > 2) ? atof(argv[2]) : 1.0;

    // ----- mesh (matches gen_per_material.py, n=2) -----
    const int nx=3;
    std::vector<Eigen::Vector3d> X;
    for (int k=0;k<3;k++) for (int j=0;j<3;j++) for (int i=0;i<3;i++)
        X.push_back(Eigen::Vector3d(i*0.5, j*0.5, k*0.5));
    auto nid=[&](int i,int j,int k){ return k*9+j*3+i; };
    int conn[8][8]; int e=0;
    for (int k=0;k<2;k++) for (int j=0;j<2;j++) for (int i=0;i<2;i++,e++) {
        int c[8]={nid(i,j,k),nid(i+1,j,k),nid(i+1,j+1,k),nid(i,j+1,k),
                  nid(i,j,k+1),nid(i+1,j,k+1),nid(i+1,j+1,k+1),nid(i,j+1,k+1)};
        for(int a=0;a<8;a++) conn[e][a]=c[a];
    }
    const int nnode=27, ndof=81, nelem=8, nip=8;

    // ----- per-element PROPS (19), quaternions from tension_8el_random.inp -----
    double qs[8][4] = {
        {-0.19326972,-0.23191724,-0.14930902, 0.94157741},
        {-0.081675831,-0.95853303,0.21273409,-0.17113649},
        {-0.17904615, 0.095634454,0.19170478, 0.96023216},
        { 0.46223193, 0.077790292,-0.51973324,-0.71426023},
        { 0.18398828, 0.97922124, 0.031112716,-0.079410838},
        { 0.32163641,-0.87915382,-0.18887598, 0.296588},
        { 0.86303465,-0.23777273, 0.37212857, 0.24526647},
        { 0.37011153,-0.52665506, 0.26883411,-0.71650549}};
    double PROPS[8][19];
    for (int el=0; el<8; el++) {
        double base[19] = {0,106750,60410,28340,2,2,1,6, 0.588235,50, 0,0,0,0, 1,2,115,34,2.25};
        if (getenv("NPOW")) base[9]=atof(getenv("NPOW"));
        std::memcpy(PROPS[el], base, sizeof(base));
        if (identity) { PROPS[el][10]=1; }
        else if (same_q) { for(int c=0;c<4;c++) PROPS[el][10+c]=qs[0][c]; }
        else            { for(int c=0;c<4;c++) PROPS[el][10+c]=qs[el][c]; }
    }
    printf("mode: %s, t_end=%g\n", identity?"identity":(same_q?"same non-identity q":"8 random q"), t_end);

    // ----- BCs from the .inp -----
    // bottom face (z=0): U3=0 ; top face (z=1): U3=0.1*t ; node(0,0,0): U1=U2=0 ; node(2,0,0): U2=0
    std::vector<int> fixed_dofs; std::vector<double> fixed_amp; // amp * t (ramp)
    for (int j=0;j<3;j++) for (int i=0;i<3;i++){ fixed_dofs.push_back(nid(i,j,0)*3+2); fixed_amp.push_back(0.0); }
    for (int j=0;j<3;j++) for (int i=0;i<3;i++){ fixed_dofs.push_back(nid(i,j,2)*3+2); fixed_amp.push_back(0.1); }
    fixed_dofs.push_back(nid(0,0,0)*3+0); fixed_amp.push_back(0.0);
    fixed_dofs.push_back(nid(0,0,0)*3+1); fixed_amp.push_back(0.0);
    fixed_dofs.push_back(nid(2,0,0)*3+1); fixed_amp.push_back(0.0);
    std::vector<char> is_fixed(ndof,0);
    for (int d: fixed_dofs) is_fixed[d]=1;

    // ----- Gauss / B-matrices (precomputed: cube elements, constant Jacobian) -----
    double gp = 1.0/std::sqrt(3.0);
    double xi[8][3]; { int a=0; for(int c=0;c<2;c++)for(int b=0;b<2;b++)for(int aa=0;aa<2;aa++,a++){
        xi[a][0]=(aa?gp:-gp); xi[a][1]=(b?gp:-gp); xi[a][2]=(c?gp:-gp);} }
    double lc[8][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
    // dN/dx: J = diag(0.25) (element size 0.5) -> dN/dx = dN/dxi / 0.25 ; detJ = 0.25^3
    double detJ = 0.015625;
    std::vector<Eigen::Matrix<double,6,24>> B(nip);
    for (int g=0; g<nip; g++) {
        B[g].setZero();
        for (int a=0; a<8; a++) {
            double dN[3];
            for (int d=0; d<3; d++) {
                double p = 0.125*lc[a][d];
                for (int o=0;o<3;o++) if(o!=d) p *= (1.0+lc[a][o]*xi[g][o]);
                dN[d] = p/0.25;
            }
            int c = 3*a;
            B[g](0,c+0)=dN[0]; B[g](1,c+1)=dN[1]; B[g](2,c+2)=dN[2];
            B[g](3,c+0)=dN[1]; B[g](3,c+1)=dN[0];             // g12
            B[g](4,c+0)=dN[2]; B[g](4,c+2)=dN[0];             // g13
            B[g](5,c+1)=dN[2]; B[g](5,c+2)=dN[1];             // g23
        }
    }

    // ----- state -----
    std::vector<std::array<double,70>> sv_conv(nelem*nip), sv_pend(nelem*nip);
    for (auto& s: sv_conv) s.fill(0.0);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(ndof), u_old = u;

    double t=0, dt=0.005, dt_max=0.005;
    int inc=0, total_cut=0, consecutive_easy=0;
    while (t < t_end - 1e-12) {
        if (inc > 3000) { printf("increment cap reached at t=%.5f\n", t); break; }
        inc++;
        int attempts=0;
        bool inc_done=false;
        while (!inc_done) {
            double t_new = std::min(t+dt, t_end);
            double dtime = t_new - t;
            Eigen::VectorXd ut = u_old;
            for (size_t c=0;c<fixed_dofs.size();c++) ut(fixed_dofs[c]) = fixed_amp[c]*t_new;

            int it=0; bool conv=false; double pn_min=1.0;
            for (it=1; it<=16; it++) {
                Eigen::VectorXd Rres = Eigen::VectorXd::Zero(ndof);
                Eigen::MatrixXd K = Eigen::MatrixXd::Zero(ndof,ndof);
                pn_min = 1.0;
                for (int el=0; el<nelem; el++) {
                    Eigen::Matrix<double,24,1> ue, ue0;
                    for (int a=0;a<8;a++) for (int d=0;d<3;d++){
                        ue(3*a+d)=ut(conn[el][a]*3+d); ue0(3*a+d)=u_old(conn[el][a]*3+d);}
                    for (int g=0; g<nip; g++) {
                        Eigen::Matrix<double,6,1> de = B[g]*(ue-ue0);
                        double DSTRAN[6]; for(int i=0;i<6;i++) DSTRAN[i]=de(i);
                        double S[6], D[36], pn;
                        int ip = el*nip+g;
                        sv_pend[ip] = sv_conv[ip];                     // start-of-increment STATEV
                        umat_call(sv_pend[ip].data(), PROPS[el], DSTRAN, dtime, S, D, pn);
                        pn_min = std::min(pn_min, pn);
                        Eigen::Map<Eigen::Matrix<double,6,1>> Sv(S);
                        Eigen::Map<Eigen::Matrix<double,6,6>> Dm(D); // column-major == our fill
                        Rres.segment<0>(0);
                        Eigen::Matrix<double,24,1> fint = B[g].transpose()*Sv*detJ;
                        Eigen::Matrix<double,24,24> ke = B[g].transpose()*Dm*B[g]*detJ;
                        for (int a=0;a<24;a++){
                            Rres(conn[el][a/3]*3+a%3) += fint(a);
                            for (int b=0;b<24;b++)
                                K(conn[el][a/3]*3+a%3, conn[el][b/3]*3+b%3) += ke(a,b);
                        }
                    }
                }
                {
                    double demax=0; for (int el=0; el<nelem; el++){
                        Eigen::Matrix<double,24,1> ue, ue0;
                        for (int a=0;a<8;a++) for (int d=0;d<3;d++){ue(3*a+d)=ut(conn[el][a]*3+d); ue0(3*a+d)=u_old(conn[el][a]*3+d);}
                        for (int g=0; g<nip; g++){ Eigen::Matrix<double,6,1> de=B[g]*(ue-ue0); demax=std::max(demax,de.cwiseAbs().maxCoeff());}
                    }
                    double rmax0=0,fref0=0; for (int d=0; d<ndof; d++){ if(is_fixed[d]) fref0=std::max(fref0,std::fabs(Rres(d))); else rmax0=std::max(rmax0,std::fabs(Rres(d)));}
                    printf("    it %2d  |de|max %.3e  rmax %.4e  fref %.4e  pn %.2f\n", it, demax, rmax0, fref0, pn_min);
                }
                if (getenv("FDPROBE") && it>=2 && it<=6) {
                    // FD-check DDSDDE at the IP with the largest strain increment at this iterate
                    int bel=0,bg=0; double bmax=-1;
                    for (int el=0; el<nelem; el++){ Eigen::Matrix<double,24,1> ue, ue0;
                        for (int a=0;a<8;a++) for (int d=0;d<3;d++){ue(3*a+d)=ut(conn[el][a]*3+d); ue0(3*a+d)=u_old(conn[el][a]*3+d);}
                        for (int g=0; g<nip; g++){ double m=(B[g]*(ue-ue0)).cwiseAbs().maxCoeff(); if(m>bmax){bmax=m;bel=el;bg=g;} } }
                    Eigen::Matrix<double,24,1> ue, ue0;
                    for (int a=0;a<8;a++) for (int d=0;d<3;d++){ue(3*a+d)=ut(conn[bel][a]*3+d); ue0(3*a+d)=u_old(conn[bel][a]*3+d);}
                    Eigen::Matrix<double,6,1> de=B[bg]*(ue-ue0);
                    double base[6]; for(int i=0;i<6;i++) base[i]=de(i);
                    double S0[6],D0[36],pn0; std::array<double,70> s0=sv_conv[bel*nip+bg];
                    umat_call(s0.data(),PROPS[bel],base,dtime,S0,D0,pn0);
                    double Dn=0; for(int k=0;k<36;k++) Dn=std::max(Dn,std::fabs(D0[k]));
                    double emax=0; int wi=-1,wj=-1;
                    for(int j=0;j<6;j++){ double h=1e-8,dp[6],dm[6],Sp[6],Sm[6],Dd[36],pp;
                        std::memcpy(dp,base,48); std::memcpy(dm,base,48); dp[j]+=h; dm[j]-=h;
                        std::array<double,70> a=sv_conv[bel*nip+bg], b=a;
                        umat_call(a.data(),PROPS[bel],dp,dtime,Sp,Dd,pp); umat_call(b.data(),PROPS[bel],dm,dtime,Sm,Dd,pp);
                        for(int i=0;i<6;i++){ double e=std::fabs((Sp[i]-Sm[i])/(2*h)-D0[i+6*j]); if(e>emax){emax=e;wi=i;wj=j;} } }
                    printf("        [FD probe el%d ip%d] |de|=%.2e  DDSDDE-FD max err %.1f MPa (rel %.2e) at (%d,%d), |S|max=%.2f\n",
                           bel+1,bg,bmax,emax,emax/Dn,wi,wj, std::max({std::fabs(S0[0]),std::fabs(S0[1]),std::fabs(S0[2])}));
                }
                if (pn_min < 1.0) break;   // local failure -> cutback
                // residual on free DOFs, Abaqus-like tolerance 0.5% of typical force
                double fref=0; for (int d=0; d<ndof; d++) if (is_fixed[d]) fref = std::max(fref, std::fabs(Rres(d)));
                double rmax=0; for (int d=0; d<ndof; d++) if (!is_fixed[d]) rmax = std::max(rmax, std::fabs(Rres(d)));
                if (rmax <= 0.005*std::max(fref,1e-8)) { conv=true; break; }
                // solve K_ff du = -R_f  (row/col elimination)
                for (int d=0; d<ndof; d++) if (is_fixed[d]) {
                    K.row(d).setZero(); K.col(d).setZero(); K(d,d)=1.0; Rres(d)=0.0;
                }
                Eigen::VectorXd du = Eigen::PartialPivLU<Eigen::MatrixXd>(K).solve(-Rres);
                ut += du;
            }
            if (conv) {
                sv_conv = sv_pend; u_old = ut; u = ut; t = t_new; inc_done = true;
                if (it<=5) { if (++consecutive_easy>=2) dt = std::min(dt*1.5, dt_max); }
                else consecutive_easy=0;
                if (inc<=8 || inc%20==0 || dt<1e-3) printf("inc %4d t=%.4f dt=%.5f iters=%d\n", inc, t, dt, it);
            } else {
                attempts++; total_cut++;
                printf("  CUTBACK inc %d attempt %d: iters=%d pn_min=%.3f -> dt=%.3e\n", inc, attempts, it, pn_min, dt*((pn_min<1.0)?std::max(pn_min,0.1):0.25));
                dt *= (pn_min<1.0)? std::max(pn_min,0.1) : 0.25;
                consecutive_easy=0;
                if (attempts>=5 || dt < 1e-10) { printf("TOO MANY ATTEMPTS at inc %d, t=%.5f (iters=%d, pn=%.2f)\n", inc, t, it, pn_min); return 1; }
            }
        }
    }
    // report S33 per element at last increment (IP-average)
    printf("DONE: t=%.3f, increments=%d, cutbacks=%d\n", t, inc, total_cut);
    printf("S33 per element (crystal->lab, from STATEV via one zero-DSTRAN probe):\n");
    for (int el=0; el<nelem; el++) {
        double s33=0;
        for (int g=0; g<nip; g++) {
            double DSTRAN[6]={0,0,0,0,0,0}, S[6], D[36], pn;
            auto sv = sv_conv[el*nip+g];
            umat_call(sv.data(), PROPS[el], DSTRAN, 1e-6, S, D, pn);
            s33 += S[2]/nip;
        }
        printf("  el %d: S33 = %8.3f MPa\n", el+1, s33);
    }
    return 0;
}
