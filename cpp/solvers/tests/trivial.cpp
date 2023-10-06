#include <catch2/catch_test_macros.hpp>

#include "solvers.h"
#include <iostream>

using namespace Kernel;
using namespace Solvers;

TEST_CASE( "Trivial solver", "[trivial]" ) {
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
    Timer timer = Timer(dt);
    auto solver = Solvers::Trivial(A, x, b, &timer);

    timer.start();
    while (solver.forward());
    timer.stop();

    std::cout << "Trivial solver took " << timer.getCurrentTimeStep() << " iterations and " << timer.getDurationInSeconds() << " seconds." << std::endl;
    
    REQUIRE( timer.getCurrentTimeStep() >= 90 );
    REQUIRE( timer.getCurrentTimeStep() <= 110 );
}