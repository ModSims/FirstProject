#include "iteration_solvers.h"
#include <iostream>

using namespace Solvers;
using namespace Eigen;

VectorXd Richardson::phi(MatrixXd& A, VectorXd& x, VectorXd& b) {
    MatrixXd I = MatrixXd::Identity(A.rows(), A.cols());
    MatrixXd M = (I - this->getOmega() * A);
    MatrixXd N = this->getOmega() * I;

    return M * x + N * b;
}