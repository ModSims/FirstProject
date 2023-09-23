#include "solvers.h"

using namespace Solvers;
using namespace Eigen;

VectorXd Trivial::phi(Ref<MatrixXd> A, Ref<VectorXd> x, Ref<VectorXd> b) {
    MatrixXd M = MatrixXd::Identity(A.rows(), A.cols()) - A;
    MatrixXd N = MatrixXd::Identity(A.rows(), A.cols());
    
    return M * x + N * b;
}
