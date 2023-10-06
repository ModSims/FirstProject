#pragma once
#include <Eigen/Dense>
#include "kernel.h"
#include "progressbar.h"

namespace Solvers {
    using namespace Eigen;
    class IterationSolver {
    public:
        IterationSolver(MatrixXd A, VectorXd x, VectorXd b, Kernel::Timer* timer);
        VectorXd solve(MatrixXd A, VectorXd x, VectorXd b);
        virtual VectorXd phi(MatrixXd& A, VectorXd& x, VectorXd& b) = 0;
        bool forward();
        void setX(VectorXd x);
        VectorXd* getResidualList();
        double * getLastResidual();
        int getMaxIterations() const;
        static double m_tolerance;
        static int m_max_iterations;
    private:
        Kernel::Timer* m_timer;
        VectorXd* m_residual_list;
        MatrixXd m_A;
        VectorXd m_x;
        VectorXd m_b;
    };

    class Trivial : public IterationSolver {
    public:
        Trivial(MatrixXd A, VectorXd x, VectorXd b, Kernel::Timer* timer) : IterationSolver(A, x, b, timer) {}
        VectorXd phi(MatrixXd& A, VectorXd& x, VectorXd& b) override;
    };

    class GaussSeidel : public IterationSolver {
    public:
        GaussSeidel(MatrixXd A, VectorXd x, VectorXd b, Kernel::Timer* timer) : IterationSolver(A, x, b, timer) {}
        VectorXd phi(MatrixXd& A, VectorXd& x, VectorXd& b) override;
    };

    class Jacobi : public IterationSolver {
    public:
        Jacobi(MatrixXd A, VectorXd x, VectorXd b, Kernel::Timer* timer) : IterationSolver(A, x, b, timer) {}
        VectorXd phi(MatrixXd& A, VectorXd& x, VectorXd& b) override;
    };
}