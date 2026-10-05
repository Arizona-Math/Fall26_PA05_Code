#ifndef PA05_PROBLEMS_H
#define PA05_PROBLEMS_H
#include "sparse.h"
#include <algorithm>
#include <cmath>
#include <vector>

struct Problem {
    int nx, ny;
    std::vector<double> d,w,e,s,n,b,exact;
    Grid grid() const { return {nx,ny,d.data(),w.data(),e.data(),s.data(),n.data()}; }
};

// material=0: smooth manufactured Poisson test, zero boundary;
// material=1: uniform plate; material=2: conducting channel with inclusion.
// Plate boundary temperature is 1-x; the volume heat source is zero.
inline Problem make_problem(int nx, int ny, int material) {
    const int count=nx*ny;
    Problem p{}; p.nx=nx; p.ny=ny;
    p.d.resize(count); p.w.resize(count); p.e.resize(count);
    p.s.resize(count); p.n.resize(count); p.b.resize(count);
    p.exact.resize(count);
    const double hx=1.0/(nx+1), hy=1.0/(ny+1), pi=std::acos(-1.0);
    auto conductivity=[material](double x,double y) {
        if (material!=2) return 1.0;
        if ((x-0.55)*(x-0.55)+(y-0.5)*(y-0.5)<0.12*0.12) return 0.25;
        return std::abs(y-0.5)<0.08 ? 8.0 : 1.0;
    };
    for(int j=0;j<ny;++j) for(int i=0;i<nx;++i) {
        const int k=j*nx+i;
        const double x=(i+1)*hx,y=(j+1)*hy;
        // Identical midpoint arithmetic at both endpoints of a shared face.
        const double gw=conductivity((i+.5)*hx,y)/(hx*hx);
        const double ge=conductivity((i+1.5)*hx,y)/(hx*hx);
        const double gs=conductivity(x,(j+.5)*hy)/(hy*hy);
        const double gn=conductivity(x,(j+1.5)*hy)/(hy*hy);
        p.d[k]=gw+ge+gs+gn;
        p.w[k]=i ? -gw : 0; p.e[k]=i+1<nx ? -ge : 0;
        p.s[k]=j ? -gs : 0; p.n[k]=j+1<ny ? -gn : 0;
        p.exact[k]=material==0 ? std::sin(pi*x)*std::sin(pi*y) : 1-x;
        if(material==0) p.b[k]=2*pi*pi*p.exact[k];
        else {
            if(i==0) p.b[k]+=gw;
            if(j==0) p.b[k]+=gs*(1-x);
            if(j+1==ny) p.b[k]+=gn*(1-x);
        }
    }
    return p;
}
#endif
