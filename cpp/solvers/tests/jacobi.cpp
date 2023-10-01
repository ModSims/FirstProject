#include <catch2/catch_test_macros.hpp>

#include "solvers.h"
#include <iostream>

using namespace Solvers;

TEST_CASE( "Jacobi solver", "[Jacobi]" ) {
    // Initialize the problem
    MatrixXd A(2, 2);
    A << 0.7, -0.4,
        -0.2,  0.5;
    
    VectorXd x(2);
    x << 21.0, -19.0;
    
    VectorXd b(2);
    b << 0.3, 0.3;
    
    // Solve the problem
    auto solver = Solvers::Jacobi();
    VectorXd y = solver.solve(A, x, b);
    
    REQUIRE( y.size() >= 45 );
    REQUIRE( y.size() <= 55 );
}