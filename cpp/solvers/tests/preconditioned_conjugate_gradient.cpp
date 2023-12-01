#include <catch2/catch_test_macros.hpp>
#include "maths.h"
#include "krylov_solvers.h"
#include <iostream>

using namespace Kernel;
using namespace Solvers;

TEST_CASE( "Preconditioned Conjugate Gradient solver", "[PreconditionedConjugateGradient]" ) {
    // Initialize the problem
    MatrixXd A(2, 2);
    A << 0.7, -0.4,
        -0.4,  0.7;
    
    VectorXd x(2);
    x << 21.0, -19.0;
    
    VectorXd b(2);
    b << 0.3, 0.3;

    // Preconditioned Conjugate Gradient
    MatrixXd P = Maths::Matrix::IncompleteCholeskyDecomposition(A).getP();
    
    // Solve the problem
    double dt = 0.0001;
    Timer timer = Timer(dt);
    auto solver = Solvers::PreconditionedConjugateGradient(A, x, b, &timer);
    solver.setPreconditioner(P);
    solver.prepareSolver();

    timer.start();
    while (solver.forward());
    timer.stop();

    std::cout << "Preconditioned Conjugate Gradient solver took " << timer.getCurrentTimeStep() << " iterations and " << timer.getDurationInSeconds() << " seconds." << std::endl;
    
    REQUIRE( timer.getCurrentTimeStep() == 2 );
}