#include "cp_umat/cp_umat_interface.hpp"
#include <cstdio>
#include <cstring>
#include <cmath>
static CPUMAT g_cp; static bool g_ready=false;
void umat_call(double* STATEV, const double* PROPS, const double* DSTRAN, double DTIME,
               double* STRESS, double* DDSDDE, double& PNEWDT){
    if(!g_ready){ int np=((int)PROPS[4]==2)?19:((int)PROPS[4]==1)?16:14; g_cp.init_params(PROPS,np); g_ready=true;}
    Grain* g=g_cp.grain();
    if(std::fabs(STATEV[69]-123456.0)>0.5){
        Eigen::Quaterniond q(PROPS[10],PROPS[11],PROPS[12],PROPS[13]); q.normalize();
        g->set_orientation(q); g->reset_state(); CPUMAT::save_state(g,STATEV); STATEV[69]=123456.0;
    } else CPUMAT::load_state(g,STATEV);
    g_cp.step(g,DSTRAN,DTIME,STRESS,DDSDDE,PNEWDT);
    CPUMAT::save_state(g,STATEV);
}
int main(int argc, char** argv){
    if(argc>1 && argv[1][0]=='e'){ /* set below */ }
    double PROPS[19]={0,106750,60410,28340, 0, 0, 1, 6, 0.588235,50, 1,0,0,0, 1,2,115,34,2.25}; // hard=0, rot=0, identity
    // args: [e]lastic-only (tc0 huge)
    extern char** environ; (void)environ;
    if(argc>1 && argv[1][0]=='e'){ PROPS[7]=1e9; printf("ELASTIC ONLY (tc0=1e9)\n"); }
    double STATEV[70]; std::memset(STATEV,0,sizeof(STATEV));
    double base[6]={-1.8e-4,-1.8e-4,5.0e-4,2e-5,-1e-5,3e-5};
    double S[6],D[36],pn; double DT=0.005;
    for(int inc=1;inc<10;inc++) umat_call(STATEV,PROPS,base,DT,S,D,pn);
    double sv0[70]; std::memcpy(sv0,STATEV,sizeof(sv0));
    double sv[70]; std::memcpy(sv,sv0,sizeof(sv));
    umat_call(sv,PROPS,base,DT,S,D,pn);
    double FD[36];
    for(int j=0;j<6;j++){
        double h=1e-8,dp[6],dm[6],Sp[6],Sm[6],Dd[36],p2;
        std::memcpy(dp,base,48); std::memcpy(dm,base,48); dp[j]+=h; dm[j]-=h;
        double a[70],b[70]; std::memcpy(a,sv0,560); std::memcpy(b,sv0,560);
        umat_call(a,PROPS,dp,DT,Sp,Dd,p2); umat_call(b,PROPS,dm,DT,Sm,Dd,p2);
        for(int i=0;i<6;i++) FD[i+6*j]=(Sp[i]-Sm[i])/(2*h);
    }
    printf("FD:\n");
    for(int i=0;i<6;i++){ for(int j=0;j<6;j++) printf("%9.1f ", FD[i+6*j]); printf("\n");}
    printf("DDSDDE - FD (MPa):\n");
    for(int i=0;i<6;i++){ for(int j=0;j<6;j++) printf("%9.2f ", D[i+6*j]-FD[i+6*j]); printf("\n");}
    printf("\nDDSDDE:\n");
    for(int i=0;i<6;i++){ for(int j=0;j<6;j++) printf("%9.1f ", D[i+6*j]); printf("\n");}
    return 0;
}
