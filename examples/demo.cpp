#include "problems.h"
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char** argv) {
    if(argc<6) {
        std::cerr<<"Usage: demo nx ny material method omega [output.csv]\n"
                 <<"material: 0 smooth Poisson, 1 uniform plate, 2 composite plate\n"
                 <<"method: 0 Jacobi, 1 SOR, 2 fast; nx,ny positive\n";
        return 2;
    }
    int nx=std::stoi(argv[1]),ny=std::stoi(argv[2]);
    int material=std::stoi(argv[3]),method=std::stoi(argv[4]);
    double omega=std::stod(argv[5]);
    if(nx<1||ny<1||material<0||material>2||method<0||method>2||omega<=0||omega>=2)
        return 2;
    auto p=make_problem(nx,ny,material);
    std::vector<double> x(nx*ny,0),ax(nx*ny);
    auto t=std::chrono::steady_clock::now();
    Result r=method==2 ? solve_fast(p.grid(),p.b.data(),x.data(),1e-10,1e-12,200000)
                      : solve(p.grid(),p.b.data(),x.data(),method,omega,1e-10,1e-12,200000);
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
    // Separate diagnostic: do not rely only on the reported residual.
    apply(p.grid(),x.data(),ax.data());
    double residual=0,error=0;
    for(int k=0;k<nx*ny;++k) {
        residual=std::max(residual,std::abs(p.b[k]-ax[k]));
        error=std::max(error,std::abs(x[k]-p.exact[k]));
    }
    std::cout<<std::setprecision(12)<<"seconds="<<seconds<<" iterations="<<r.iterations
             <<" residual="<<residual<<" reported_residual="<<r.residual
             <<" converged="<<r.converged;
    if(material!=2) std::cout<<" error="<<error;
    std::cout<<'\n';
    if(argc>6) {
        std::ofstream out(argv[6]);
        out<<"x,y,u\n"<<std::setprecision(17);
        for(int j=0;j<ny;++j) for(int i=0;i<nx;++i)
            out<<(i+1.0)/(nx+1)<<','<<(j+1.0)/(ny+1)<<','<<x[j*nx+i]<<'\n';
    }
    return r.converged ? 0 : 1;
}
