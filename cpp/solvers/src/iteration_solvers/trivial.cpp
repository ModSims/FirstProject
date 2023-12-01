#include "iteration_solvers.h"

using namespace Solvers;
using namespace Eigen;

VectorXd Trivial::phi(MatrixXd& A, VectorXd& x, VectorXd& b) {
    MatrixXd M = MatrixXd::Identity(A.rows(), A.cols()) - A;
    MatrixXd N = MatrixXd::Identity(A.rows(), A.cols());
    
    return M * x + N * b;
}
