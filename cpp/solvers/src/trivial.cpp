#include "solvers.h"

int Trivial::solve(int M, int N, int x, int b) {
    return M * x + N * b;
}