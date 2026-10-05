#include "problems.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
int main(int argc,char** argv) {
    int n=argc>1 ? std::stoi(argv[1]) : 128;
    int sweeps=argc>2 ? std::stoi(argv[2]) : 100;
    if(n<1||sweeps<1) return 2;
    auto p=make_problem(n,n,0);
    std::vector<double> x(n*n,0),y(n*n);
    auto t=std::chrono::steady_clock::now();
    for(int k=0;k<sweeps;++k) {
        jacobi_sweep(p.grid(),p.b.data(),x.data(),y.data());
        x.swap(y);
    }
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
    std::cout<<std::setprecision(17)<<seconds<<' '
             <<std::accumulate(x.begin(),x.end(),0.0)<<'\n';
}
