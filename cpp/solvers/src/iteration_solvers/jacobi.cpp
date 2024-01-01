#include "iteration_solvers.h"
#include <iostream>

using namespace Solvers;
using namespace Eigen;

void Jacobi::prepareSolver() {
    // Extract diagonal elements of A
    VectorXd diagA = this->getData()->A.diagonal();

    // Check for division by zero in diagA
    for (int i = 0; i < this->getData()->A.rows(); ++i) {
        if (diagA(i) == 0) {
            // Handle division by zero by setting a small epsilon value
            diagA(i) = 1e-10;
        }
    }

    // Calculate the reciprocal of diagonal elements
    this->m_D_inv = diagA.cwiseInverse();
}

VectorXd Jacobi::phi(SolverData *data) {
    // Calculate the residual term (b - Ax)
    VectorXd residual = data->b - data->A * data->x;

    // Perform Jacobi iteration directly on vector x
    data->x += this->getOmega() * (this->m_D_inv.cwiseProduct(residual));

    return residual;
}

