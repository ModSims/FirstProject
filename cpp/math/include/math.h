#pragma once
#include <Eigen/Dense>

namespace Math {
    using namespace Eigen;

    class Radians {
    public:
        double calculateSpectralRadius(Eigen::MatrixXd& A);
        double gelfands_spectral_approximation(Eigen::MatrixXd& matrix, int max_iterations = 100);
    };
}