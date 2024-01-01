#include "iteration_solvers.h"
#include <iostream>

using namespace Solvers;
using namespace Eigen;

void Richardson::prepareSolver() {
    return;
}

VectorXd Richardson::phi(SolverData *data) {
    const int size = data->A.rows();
    VectorXd& x = data->x;
    const VectorXd& b = data->b;

    // Identity matrix
    MatrixXd I = MatrixXd::Identity(size, size);

    // Calculate matrices M and N
    MatrixXd M = I - this->getOmega() * data->A;
    MatrixXd N = this->getOmega() * I;

    // Perform Richardson iteration directly on vector x
    x = M * x + N * b;

    // Calculate the residual term (b - Ax)
    VectorXd residual = b - data->A * x;

    return residual;
}
