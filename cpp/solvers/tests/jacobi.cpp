#include <catch2/catch_test_macros.hpp>

#include "iteration_solvers.h"
#include <iostream>

using namespace Kernel;
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
    double dt = 0.0001;
    int max_iterations = 100;
    double omega = 1.0;
    Timer timer = Timer(dt);
    auto solver = Solvers::Jacobi(A, x, b, max_iterations, omega, &timer);

    timer.start();
    while (solver.forward());
    timer.stop();

    std::cout << "Jacobi solver took " << timer.getCurrentTimeStep() << " iterations and " << timer.getDurationInSeconds() << " seconds." << std::endl;
    
    REQUIRE( timer.getCurrentTimeStep() >= 45 );
    REQUIRE( timer.getCurrentTimeStep() <= 55 );
}
