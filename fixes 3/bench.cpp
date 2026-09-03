// Per-IP benchmark: realistic uniaxial path with lateral contraction + shears, random orientation
#include "cp_umat/cp_umat_interface.hpp"
#include <cstdio>
#include <cstring>
#include <chrono>
extern "C" void umat_(double*,double*,double*,double&,double&,double&,double&,double*,double*,double*,double*,double*,double*,double*,double&,double&,double*,double*,char*,int&,int&,int&,int&,double*,int&,double*,double*,double&,double*,double*,double*,int&,int&,int&,int&,int*,int&);

int main(int argc,char**argv){
    double npow = argc>1? atof(argv[1]) : 50;
    int ninc = argc>2? atoi(argv[2]) : 600;
    double PROPS[19]={0,106750,60410,28340,2,2,1,6, 0.588235,npow, -0.12523701,0.09397545,0.59307825,0.78977370, 1,2,115,34,2.25};
    double STATEV[70]; std::memset(STATEV,0,sizeof(STATEV));
    double S[6]={0},D[36],dummy=0,dsum[6]={0},z=0; int i3=3,i6=6,np=19,one=1,jstep[4]={1,0,0,0};
    double dt=0.3/ninc;  // like the 343el step: total time 0.3
    double base[6]={-1.7e-4,-1.7e-4,5.0e-4, 3e-5,-2e-5,4e-5}; for(int i=0;i<6;i++) base[i]*= (0.005/dt)*(dt/0.005); // strain per increment ~ 5e-4 scaled by dt
    for(int i=0;i<6;i++) base[i]*= dt/0.005;
    double tot=0; auto t0=std::chrono::high_resolution_clock::now();
    double pn; char name[8]="CP1";
    for(int inc=1;inc<=ninc;inc++){
        // emulate ~4 equilibrium iterations per increment: same start STATEV, slightly different DSTRAN
        double sv0[70]; std::memcpy(sv0,STATEV,sizeof(sv0));
        for(int it=0;it<4;it++){
            double ds[6]; for(int i=0;i<6;i++) ds[i]=base[i]*(1.0+0.02*(it==0?0.5:(it%2?1:-1)));
            std::memcpy(STATEV,sv0,sizeof(sv0));
            umat_(S,STATEV,D,dummy,dummy,dummy,dummy,nullptr,nullptr,nullptr,nullptr,ds,&z,&dt,dummy,dummy,nullptr,nullptr,name,i3,i3,i6,np,PROPS,np,nullptr,nullptr,pn,nullptr,nullptr,nullptr,one,one,one,one,jstep,inc);
        }
    }
    auto t1=std::chrono::high_resolution_clock::now();
    double us=std::chrono::duration<double,std::micro>(t1-t0).count()/(ninc*4);
    printf("n=%g: %.1f us per UMAT call, S33=%.2f MPa\n", npow, us, S[2]);
}
