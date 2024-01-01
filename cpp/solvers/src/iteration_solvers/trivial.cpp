#include "iteration_solvers.h"

using namespace Solvers;
using namespace Eigen;

void Trivial::prepareSolver() {
    return;
}

VectorXd Trivial::phi(SolverData *data) {
    const int size = data->A.rows();
    VectorXd& x = data->x;
    const VectorXd& b = data->b;

    // Identity matrices
    MatrixXd M = MatrixXd::Identity(size, size) - data->A;
    MatrixXd N = MatrixXd::Identity(size, size);

    // Perform Trivial iteration directly on vector x
    x = M * x + N * b;

    // Calculate the residual term (b - Ax)
    VectorXd residual = b - data->A * x;

    return residual;
}
