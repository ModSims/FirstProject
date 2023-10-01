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
    const double tolerance = 1e-14;
    const int max_iterations = 100000;

    VectorXd resid_list(max_iterations);  // Pre-allocate space for residuals
    int iteration = 0;
    double resid = 0;

    // print initial iteration and residual
    std::cout << "Iteration: " << iteration << " Residual: " << resid << std::endl;

    do {
        VectorXd x_new = phi(A, x, b);
        resid = (b - A * x_new).norm();
        x = x_new;

        resid_list(iteration) = resid;

        // print current iteration and residual
        std::cout << "Iteration: " << iteration << " Residual: " << resid << std::endl;
        
        ++iteration;
    } while (resid >= tolerance && iteration < max_iterations);

    // Resize resid_list to the actual number of iterations
    resid_list.conservativeResize(iteration);

    return resid_list;
}