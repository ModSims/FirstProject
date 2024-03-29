#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "krylov_solvers.h"
#include "kernel.h"
#include "maths.h"

using namespace Eigen;
using namespace Solvers;
using namespace Kernel;

MatrixXd generateMatrix(int n) {
    double h = 1.0 / (n - 1);

    MatrixXd A = MatrixXd::Zero(n, n);

    for (int i = 0; i < n; ++i) {
        A(i, i) = 2.0;
        if (i > 0) {
            A(i, i - 1) = -1.0;
        }
        if (i < n - 1) {
            A(i, i + 1) = -1.0;
        }
    }

    return A / (h * h);
}

SolverData prepareProblem(int n) {
    SolverData data;
    data.A = generateMatrix(n);
    data.x = VectorXd::Zero(n);
    data.b = VectorXd::Zero(n);
    data.b(0) = 0.0;  // Dirichlet boundary condition
    data.b(n - 1) = 1.0;  // Dirichlet boundary condition

    return data;
}

int main(int argc, char** argv) {
    auto P256 = prepareProblem(256);
    auto P512 = prepareProblem(512);
    auto P1024 = prepareProblem(1024);
    auto P2048 = prepareProblem(2048);

    double dt = 0.0001;

    // Conjugate Gradient

    Timer timer256 = Timer(dt);
    auto solver256 = Solvers::ConjugateGradient(&P256, &timer256);
    solver256.prepareSolver();

    Timer timer512 = Timer(dt);
    auto solver512 = Solvers::ConjugateGradient(&P512, &timer512);
    solver512.prepareSolver();

    Timer timer1024 = Timer(dt);
    auto solver1024 = Solvers::ConjugateGradient(&P1024, &timer1024);
    solver1024.prepareSolver();

    Timer timer2048 = Timer(dt);
    auto solver2048 = Solvers::ConjugateGradient(&P2048, &timer2048);
    solver2048.prepareSolver();

    timer256.start();
    while (solver256.forward());
    timer256.stop();

    timer512.start();
    while (solver512.forward());
    timer512.stop();

    timer1024.start();
    while (solver1024.forward());
    timer1024.stop();

    timer2048.start();
    while (solver2048.forward());
    timer2048.stop();

    saveVector("alphas256.dat", solver256.getAlphas());
    saveVector("alphas512.dat", solver512.getAlphas());
    saveVector("alphas1024.dat", solver1024.getAlphas());
    saveVector("alphas2048.dat", solver2048.getAlphas());

    saveVector("x256.dat", solver256.getX());
    saveVector("x512.dat", solver512.getX());
    saveVector("x1024.dat", solver1024.getX());
    saveVector("x2048.dat", solver2048.getX());

    std::cout << "Conjugate Gradient solver took " << timer256.getCurrentTimeStep() << " iterations and " << timer256.getDurationInSeconds() << " seconds for n = 256." << std::endl;
    std::cout << "Conjugate Gradient solver took " << timer512.getCurrentTimeStep() << " iterations and " << timer512.getDurationInSeconds() << " seconds for n = 512." << std::endl;
    std::cout << "Conjugate Gradient solver took " << timer1024.getCurrentTimeStep() << " iterations and " << timer1024.getDurationInSeconds() << " seconds for n = 1024." << std::endl;
    std::cout << "Conjugate Gradient solver took " << timer2048.getCurrentTimeStep() << " iterations and " << timer2048.getDurationInSeconds() << " seconds for n = 2048." << std::endl;

    // Preconditioned Conjugate Gradient

    Timer timer_p256 = Timer(dt);
    auto solver_p256 = Solvers::PreconditionedConjugateGradient(&P256, &timer_p256);
    MatrixXd C256 = Maths::Matrix::IncompleteCholeskyDecomposition(P256.A).getP();
    solver_p256.setPreconditioner(C256);
    solver_p256.prepareSolver();

    Timer timer_p512 = Timer(dt);
    auto solver_p512 = Solvers::PreconditionedConjugateGradient(&P512, &timer_p512);
    MatrixXd C512 = Maths::Matrix::IncompleteCholeskyDecomposition(P512.A).getP();
    solver_p512.setPreconditioner(C512);
    solver_p512.prepareSolver();

    Timer timer_p1024 = Timer(dt);
    auto solver_p1024 = Solvers::PreconditionedConjugateGradient(&P1024, &timer_p1024);
    MatrixXd C1024 = Maths::Matrix::IncompleteCholeskyDecomposition(P1024.A).getP();
    solver_p1024.setPreconditioner(C1024);
    solver_p1024.prepareSolver();

    Timer timer_p2048 = Timer(dt);
    auto solver_p2048 = Solvers::PreconditionedConjugateGradient(&P2048, &timer_p2048);
    MatrixXd C2048 = Maths::Matrix::IncompleteCholeskyDecomposition(P2048.A).getP();
    solver_p2048.setPreconditioner(C2048);
    solver_p2048.prepareSolver();

    timer_p256.start();
    while (solver_p256.forward());
    timer_p256.stop();

    timer_p512.start();
    while (solver_p512.forward());
    timer_p512.stop();

    timer_p1024.start();
    while (solver_p1024.forward());
    timer_p1024.stop();

    timer_p2048.start();
    while (solver_p2048.forward());
    timer_p2048.stop();


    saveMatrix("Preconditioner256.dat", &C256);
    saveMatrix("Preconditioner512.dat", &C512);
    saveMatrix("Preconditioner1024.dat", &C1024);
    saveMatrix("Preconditioner2048.dat", &C2048);


    saveVector("palphas256.dat", solver_p256.getAlphas());
    saveVector("palphas512.dat", solver_p512.getAlphas());
    saveVector("palphas1024.dat", solver_p1024.getAlphas());
    saveVector("palphas2048.dat", solver_p2048.getAlphas());

    saveVector("px256.dat", solver_p256.getX());
    saveVector("px512.dat", solver_p512.getX());
    saveVector("px1024.dat", solver_p1024.getX());
    saveVector("px2048.dat", solver_p2048.getX());

    std::cout << "Preconditioned Conjugate Gradient solver took " << timer_p256.getCurrentTimeStep() << " iterations and " << timer_p256.getDurationInSeconds() << " seconds for n = 256." << std::endl;
    std::cout << "Preconditioned Conjugate Gradient solver took " << timer_p512.getCurrentTimeStep() << " iterations and " << timer_p512.getDurationInSeconds() << " seconds for n = 512." << std::endl;
    std::cout << "Preconditioned Conjugate Gradient solver took " << timer_p1024.getCurrentTimeStep() << " iterations and " << timer_p1024.getDurationInSeconds() << " seconds for n = 1024." << std::endl;
    std::cout << "Preconditioned Conjugate Gradient solver took " << timer_p2048.getCurrentTimeStep() << " iterations and " << timer_p2048.getDurationInSeconds() << " seconds for n = 2048." << std::endl;

    return 0;
}