#pragma once
#include <Eigen/Dense>

namespace Solvers {
    using namespace Eigen;
    class Linear {
    public:
        virtual VectorXd phi(
            Ref<MatrixXd> A,
            Ref<VectorXd> x,
            Ref<VectorXd> b
        ) = 0;
        VectorXd solve(
            Ref<MatrixXd> A,
            Ref<VectorXd> x,
            Ref<VectorXd> b
        );
    };

    class Trivial : public Linear {
    public:
        VectorXd phi(
            Ref<MatrixXd> A,
            Ref<VectorXd> x,
            Ref<VectorXd> b
        ) override;
    };

    class GaussSeidel : public Linear {
    public:
        VectorXd phi(
            Ref<MatrixXd> A,
            Ref<VectorXd> x,
            Ref<VectorXd> b
        ) override;
    };

    class Jacobi : public Linear {
    public:
        VectorXd phi(
            Ref<MatrixXd> A,
            Ref<VectorXd> x,
            Ref<VectorXd> b
        ) override;
    };
}