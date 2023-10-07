#include <iostream>
#include <Eigen/Dense>
#include "maths.h"

double Maths::Spectral::getSpectralRadiusViaEigen(Eigen::MatrixXd& A) {
    Eigen::EigenSolver<Eigen::MatrixXd> solver(A);
    Eigen::VectorXd eigenvalues = solver.eigenvalues().real(); // Get real part of eigenvalues

    double spectralRadius = eigenvalues.cwiseAbs().maxCoeff(); // Calculate spectral radius

    return spectralRadius;
}

double Maths::Spectral::getGelfandsSpectralApproximation(Eigen::MatrixXd& matrix, int max_iterations) {
    Eigen::MatrixXd matrix_power = matrix;
    double prev_term = 0.0;

    for (int k = 1; k <= max_iterations; ++k) {
        matrix_power = matrix_power * matrix;
        double f_norm = matrix_power.norm();
        double term = pow(f_norm, 1.0 / k);

        // Check for convergence
        if (k > 1 && abs(term - prev_term) < 1e-6) {
            return term;
        }

        prev_term = term;
    }

    // Return the last computed term if it didn't converge
    return prev_term;
}