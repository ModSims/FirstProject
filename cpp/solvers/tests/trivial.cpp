#include "solvers.h"
#include <iostream>

using namespace Eigen;

// main
int main() {
    // Initialize the problem
    MatrixXd A(2, 2);
    A << 0.7, -0.4,
        -0.2,  0.5;
    
    VectorXd x(2);
    x << 21.0, -19.0;
    
    VectorXd b(2);
    b << 0.3, 0.3;
    
    // Solve the problem
    auto solver = Solvers::Trivial();
    VectorXd y = solver.solve(A, x, b);
    
    // Print the length of the error list
    std::cout << y.size() << std::endl;
    
    return 0;
}