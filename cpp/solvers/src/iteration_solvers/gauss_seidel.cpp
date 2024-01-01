#include "iteration_solvers.h"

using namespace Solvers;
using namespace Eigen;

void GaussSeidel::prepareSolver() {
    return;
}

VectorXd GaussSeidel::phi(SolverData *data) {
    VectorXd& x = data->x;
    const VectorXd& b = data->b;

    // Extract diagonal and strictly lower triangular parts of A
    MatrixXd D = data->A.diagonal().asDiagonal();
    MatrixXd L = data->A.triangularView<StrictlyLower>();

    // Calculate the inverse of (D + omega*L)
    MatrixXd invDL = (D + this->getOmega() * L).inverse();

    // Calculate the matrices M and N
    MatrixXd M = invDL * ((1 - this->getOmega()) * D - this->getOmega() * (data->A - L - D));
    MatrixXd N = this->getOmega() * invDL;

    // Perform Gauss-Seidel iteration directly on vector x
    x = M * x + N * b;

    // Calculate the residual term (b - Ax)
    VectorXd residual = b - data->A * x;

    return residual;
}
