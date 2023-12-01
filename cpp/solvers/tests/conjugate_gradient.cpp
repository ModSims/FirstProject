#include <catch2/catch_test_macros.hpp>
#include "maths.h"
#include "krylov_solvers.h"
#include <iostream>

using namespace Kernel;
using namespace Solvers;

TEST_CASE( "Conjugate Gradient solver", "[ConjugateGradient]" ) {
    // Initialize the problem
    MatrixXd A(2, 2);
    A << 0.7, -0.4,
        -0.4,  0.7;
    
    VectorXd x(2);
    x << 21.0, -19.0;
    
    VectorXd b(2);
    b << 0.3, 0.3;
    
    // Solve the problem
    double dt = 0.0001;
    Timer timer = Timer(dt);
    auto solver = Solvers::ConjugateGradient(A, x, b, &timer);
    solver.prepareSolver();

    timer.start();
    while (solver.forward());
    timer.stop();

    std::cout << "Conjugate Gradient solver took " << timer.getCurrentTimeStep() << " iterations and " << timer.getDurationInSeconds() << " seconds." << std::endl;
    
    REQUIRE( timer.getCurrentTimeStep() == 1 );
}