#include <catch2/catch_test_macros.hpp>

#include "iteration_solvers.h"
#include <iostream>

using namespace Kernel;
using namespace Solvers;

TEST_CASE( "Trivial solver", "[trivial]" ) {
    // Initialize the problem
    SolverData data;
    data.A = MatrixXd(2, 2);
    data.A << 0.7, -0.4,
        -0.2,  0.5;
    
    data.x = VectorXd(2);
    data.x << 21.0, -19.0;
    
    data.b = VectorXd(2);
    data.b << 0.3, 0.3;
    
    // Solve the problem
    double dt = 0.0001;
    int max_iterations = 100;
    double omega = 1.0;
    Timer timer = Timer(dt);
    auto solver = Solvers::Trivial(&data, max_iterations, omega, &timer);
    solver.prepareSolver();

    timer.start();
    while (solver.forward());
    timer.stop();

    std::cout << "Trivial solver took " << timer.getCurrentTimeStep() << " iterations and " << timer.getDurationInSeconds() << " seconds." << std::endl;
    
    REQUIRE( timer.getCurrentTimeStep() >= 90 );
    REQUIRE( timer.getCurrentTimeStep() <= 110 );
}