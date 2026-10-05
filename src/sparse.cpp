#include "sparse.h"
#include <limits>

// This file compiles before you implement it; the tests should fail.
void apply(const Grid& a, const double* x, double* y) {
    (void)x;
    for (int p=0; p<a.nx*a.ny; ++p) y[p]=0.0; // TODO
}
void jacobi_sweep(const Grid& a, const double* b,
                  const double* old_x, double* new_x) {
    (void)b;
    for (int p=0; p<a.nx*a.ny; ++p) new_x[p]=old_x[p]; // TODO
}
void sor_sweep(const Grid& a, const double* b, double* x, double omega) {
    (void)a; (void)b; (void)x; (void)omega; // TODO
}
Result solve(const Grid& a, const double* b, double* x, int method,
             double omega, double rtol, double atol, int max_iter) {
    (void)a; (void)b; (void)x; (void)method; (void)omega;
    (void)rtol; (void)atol; (void)max_iter; // TODO
    return {0, std::numeric_limits<double>::infinity(), false};
}
Result solve_fast(const Grid& a, const double* b, double* x,
                  double rtol, double atol, int max_iter) {
    // Working fallback once solve is implemented. Optimize after correctness.
    return solve(a,b,x,1,1.5,rtol,atol,max_iter);
}
