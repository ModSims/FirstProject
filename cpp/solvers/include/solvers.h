#pragma once
#include <Eigen/Dense>

namespace Solvers {
    using namespace Eigen;
    class Linear {
    public:
        virtual VectorXd phi(MatrixXd& A, VectorXd& x, VectorXd& b) = 0;
        VectorXd solve(MatrixXd& A, VectorXd& x, VectorXd& b);
    };

    class Trivial : public Linear {
    public:
        VectorXd phi(MatrixXd& A, VectorXd& x, VectorXd& b) override;
    };

    class GaussSeidel : public Linear {
    public:
        VectorXd phi(MatrixXd& A, VectorXd& x, VectorXd& b) override;
    };

    class Jacobi : public Linear {
    public:
        VectorXd phi(MatrixXd& A, VectorXd& x, VectorXd& b) override;
    };
}