#include "solvers.h"
#include <iostream>
#include <vector>

using namespace Solvers;
using namespace Eigen;


/*
A: given matrix
x: initial guess
b: right hand side
*/
VectorXd Linear::solve(MatrixXd& A, VectorXd& x, VectorXd& b) {
    VectorXd resid_list;
    double resid = 0;
    int max_iterations = 100000;
    int iteration = 0;

    // first iteration for while loop to work properly
    VectorXd x_new = phi(A, x, b);
    resid = (x_new - x).norm();
    x = x_new;
    resid_list.conservativeResize(resid_list.size() + 1);
    resid_list(resid_list.size() - 1) = resid;
    while (resid >= 1e-14 && iteration < max_iterations) {
        x_new = phi(A, x, b);
        resid = (x_new - x).norm();
        x = x_new;

        resid_list.conservativeResize(resid_list.size() + 1);
        resid_list(resid_list.size() - 1) = resid;

        ++iteration;
    }
    std::cout << "Number of iterations: " << iteration;
    std::cout << " Residual: " << resid << std::endl;

    return resid_list;
}