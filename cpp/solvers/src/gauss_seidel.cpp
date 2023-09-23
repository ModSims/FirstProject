#include "solvers.h"

using namespace Solvers;
using namespace Eigen;

VectorXd GaussSeidel::phi(Ref<MatrixXd> A, Ref<VectorXd> x, Ref<VectorXd> b) {
    MatrixXd D = A.diagonal().asDiagonal();
    MatrixXd L = A.triangularView<StrictlyLower>();
    MatrixXd M = (D + L).inverse() * (D + L - A);
    MatrixXd N = (D + L).inverse();

    return M * x + N * b;
}
