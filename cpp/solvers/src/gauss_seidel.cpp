#include "solvers.h"

using namespace Solvers;
using namespace Eigen;

VectorXd GaussSeidel::phi(MatrixXd& A, VectorXd& x, VectorXd& b) {
    MatrixXd D = A.diagonal().asDiagonal();
    MatrixXd L = A.triangularView<StrictlyLower>();

    MatrixXd M = ((D + this->getOmega()*L).inverse()) * ((1 - this->getOmega())*D - this->getOmega()*(A - L - D));
    MatrixXd N = this->getOmega()*(D + this->getOmega()*L).inverse();

    return M * x + N * b;
}
