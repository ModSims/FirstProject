#include "solvers.h"

using namespace Solvers;
using namespace Eigen;

VectorXd Jacobi::phi(Ref<MatrixXd> A, Ref<VectorXd> x, Ref<VectorXd> b) {
    MatrixXd D = A.diagonal().asDiagonal();
    MatrixXd M = D.inverse() * (D - A);
    MatrixXd N = D.inverse(); 

    return M * x + N * b;
}
