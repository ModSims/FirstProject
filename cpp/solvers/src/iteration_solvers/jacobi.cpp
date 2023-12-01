#include "iteration_solvers.h"
#include <iostream>

using namespace Solvers;
using namespace Eigen;

VectorXd Jacobi::phi(MatrixXd& A, VectorXd& x, VectorXd& b) {
    MatrixXd I = MatrixXd::Identity(A.rows(), A.cols());
    MatrixXd D = A.diagonal().asDiagonal();
    // Inverse the diagonal matrix D by taking the reciprocal of each element (without method .inverse())
    MatrixXd D_inv = D;
    for (int i = 0; i < D_inv.rows(); ++i) {
        if (D_inv(i, i) != 0) {
            D_inv(i, i) = 1 / D_inv(i, i);
        }
    }

    // omega relaxation
    MatrixXd M = (I - this->getOmega() * D_inv * A);
    MatrixXd N = this->getOmega() * D_inv;
    return M * x + N * b;
}
