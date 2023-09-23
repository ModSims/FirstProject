#include "solvers.h"
#include <iostream>
#include <vector>

using namespace Solvers;
using namespace Eigen;

VectorXd Linear::solve(MatrixXd& A, VectorXd& x, VectorXd& b) {
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
        if (error_list.size() > 5) {
            bool same = true;
            for (int i = 0; i < 5; i++) {
                if (error_list(error_list.size() - 1 - i) != error_list(error_list.size() - 2 - i)) {
                    same = false;
                    break;
                }
            }
            if (same) {
                break;
            }
        }
    }
    return error_list;
}
