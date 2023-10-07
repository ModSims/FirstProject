#pragma once
#include <Eigen/Dense>

namespace Maths {
    using namespace Eigen;

    namespace Spectral {
        double getSpectralRadiusViaEigen(Eigen::MatrixXd& A);
        double getGelfandsSpectralApproximation(Eigen::MatrixXd& matrix, int max_iterations = 100);
    };

    namespace Matrix {
        double getConditionNumber(Eigen::MatrixXd& A);
        bool isSymmetric(Eigen::MatrixXd& A);
        bool isDiagonallyDominant(Eigen::MatrixXd& A);
        bool isSPD(Eigen::MatrixXd& A);
        class IncompleteCholeskyDecomposition {
        public:
            IncompleteCholeskyDecomposition(Eigen::MatrixXd& A);
            IncompleteCholeskyDecomposition(const IncompleteCholeskyDecomposition& other) = delete;
            IncompleteCholeskyDecomposition(IncompleteCholeskyDecomposition&& other) : L(other.L), F(other.F) {}
            Eigen::MatrixXd getL();
            Eigen::MatrixXd getF();
        private:
            Eigen::MatrixXd L;
            Eigen::MatrixXd F;
        };
    };
}