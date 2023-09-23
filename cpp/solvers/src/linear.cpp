#include "solvers.h"
#include <iostream>
#include <vector>

using namespace Solvers;
using namespace Eigen;

VectorXd Linear::solve(Ref<MatrixXd> A, Ref<VectorXd> x, Ref<VectorXd> b) {
    VectorXd error_list;
    double error_norm = 0;
    
    while (true) {
        VectorXd x_new = phi(A, x, b);
        error_norm = (x_new - x).norm();
        x = x_new;
        if (error_norm > 1e-14 && error_norm < 1e14) {
            error_list.conservativeResize(error_list.size() + 1);
            error_list(error_list.size() - 1) = error_norm;
        } else {
            break;
        }
    }
    return error_list;
}
