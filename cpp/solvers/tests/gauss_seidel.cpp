#include <catch2/catch_test_macros.hpp>

#include "iteration_solvers.h"
#include <iostream>

using namespace Kernel;
using namespace Solvers;

TEST_CASE( "GaussSeidel solver", "[GaussSeidel]" ) {
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
    auto solver = Solvers::GaussSeidel(&data, max_iterations, omega, &timer);
    solver.prepareSolver();

    timer.start();
    while (solver.forward());
    timer.stop();

    std::cout << "Gauss Seidel solver took " << timer.getCurrentTimeStep() << " iterations and " << timer.getDurationInSeconds() << " seconds." << std::endl;
    
    REQUIRE( timer.getCurrentTimeStep() >= 23 );
    REQUIRE( timer.getCurrentTimeStep() <= 27 );
}

TEST_CASE( "GaussSeidel solver (relaxation)", "[GaussSeidel]" ) {
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
    double omega = 1.0648;
    Timer timer = Timer(dt);
    auto solver = Solvers::GaussSeidel(&data, max_iterations, omega, &timer);

    timer.start();
    while (solver.forward());
    timer.stop();

    std::cout << "Gauss Seidel (Relaxation) solver took " << timer.getCurrentTimeStep() << " iterations and " << timer.getDurationInSeconds() << " seconds." << std::endl;
    
    REQUIRE( timer.getCurrentTimeStep() >= 14 );
    REQUIRE( timer.getCurrentTimeStep() <= 16 );
}