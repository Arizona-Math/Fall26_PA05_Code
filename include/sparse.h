#ifndef PA05_SPARSE_H
#define PA05_SPARSE_H

// Five occupied diagonals, p = j*nx+i (x varies fastest).
// Off-diagonal entries are NONPOSITIVE matrix entries, not conductances.
// Missing neighbors have coefficient zero. Arrays have nx*ny entries.
struct Grid {
    int nx, ny;
    const double *d, *w, *e, *s, *n;
};

struct Result {
    int iterations;           // completed sweeps (both colors count as one)
    double residual;          // ||b-Ax||_infinity at the RETURNED x
    bool converged;
};

// Inputs satisfy the assumptions in README.md. Arithmetic is double.
// Preserve all input arrays except the designated outputs. No printing.
void apply(const Grid& a, const double* x, double* y);
void jacobi_sweep(const Grid& a, const double* b,
                  const double* old_x, double* new_x);
// Lexicographic order p=0,...,nx*ny-1. omega=1 is Gauss-Seidel.
void sor_sweep(const Grid& a, const double* b, double* x, double omega);

// method=0: Jacobi; method=1: lexicographic SOR with given omega.
// Test the initial residual, then the residual after EVERY complete sweep.
// Stop at ||b-Ax||_inf <= atol + rtol*||b||_inf, or max_iter sweeps.
Result solve(const Grid& a, const double* b, double* x, int method,
             double omega, double rtol, double atol, int max_iter);

// Your performance entry. Same tolerance, initial-guess and budget contract.
// May check residual less often; always report the true final residual.
// Implement using Jacobi/GS/SOR; reordering and up to four threads allowed.
Result solve_fast(const Grid& a, const double* b, double* x,
                  double rtol, double atol, int max_iter);
#endif
