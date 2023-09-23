#include "solvers.h"

int Jacobi::solve(int M, int N, int x, int b) {
    return M * x + N * b;
}