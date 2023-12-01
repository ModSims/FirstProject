#pragma once
#include <Eigen/Dense>
#include "kernel.h"

namespace Solvers
{
    using namespace Eigen;
    class IterationSolver
    {
    public:
        IterationSolver(MatrixXd A, VectorXd x, VectorXd b, int max_iterations, double omega, Kernel::Timer *timer);
        IterationSolver(IterationSolver &solver) = delete;
        IterationSolver &operator=(IterationSolver &solver) = delete;
        IterationSolver(IterationSolver&& other) noexcept
            : m_timer(other.m_timer),
            m_residuals(other.m_residuals),
            m_A(std::move(other.m_A)),
            m_x(std::move(other.m_x)),
            m_b(std::move(other.m_b)),
            m_max_iterations(other.m_max_iterations),
            m_omega(other.m_omega)
        {
            // Nullify other's pointers
            other.m_timer = nullptr;
            other.m_residuals = nullptr;
            other.m_max_iterations = 0;
            other.m_omega = 0;
        }
        virtual VectorXd phi(MatrixXd &A, VectorXd &x, VectorXd &b) = 0;
        bool forward();
        void setX(VectorXd x);
        VectorXd *getResiduals();
        double *getLastResidual();
        int getMaxIterations() const;
        double getTolerance();
        double getOmega();
        VectorXd *getX();
        const double m_tolerance = 1e-14;

    private:
        Kernel::Timer *m_timer;
        VectorXd *m_residuals;
        MatrixXd m_A;
        VectorXd m_x;
        VectorXd m_b;
        int m_max_iterations;
        double m_omega;
    };

    class Trivial : public IterationSolver
    {
    public:
        Trivial(MatrixXd A, VectorXd x, VectorXd b, int max_iterations, double omega, Kernel::Timer *timer) : IterationSolver(A, x, b, max_iterations, omega, timer) {}
        VectorXd phi(MatrixXd &A, VectorXd &x, VectorXd &b) override;
    };

    class GaussSeidel : public IterationSolver
    {
    public:
        GaussSeidel(MatrixXd A, VectorXd x, VectorXd b, int max_iterations, double omega, Kernel::Timer *timer) : IterationSolver(A, x, b, max_iterations, omega, timer) {}
        VectorXd phi(MatrixXd &A, VectorXd &x, VectorXd &b) override;
    };

    class Jacobi : public IterationSolver
    {
    public:
        Jacobi(MatrixXd A, VectorXd x, VectorXd b, int max_iterations, double omega, Kernel::Timer *timer) : IterationSolver(A, x, b, max_iterations, omega, timer) {}
        VectorXd phi(MatrixXd &A, VectorXd &x, VectorXd &b) override;
    };

    class Richardson : public IterationSolver
    {
    public:
        Richardson(MatrixXd A, VectorXd x, VectorXd b, int max_iterations, double omega, Kernel::Timer *timer) : IterationSolver(A, x, b, max_iterations, omega, timer) {}
        VectorXd phi(MatrixXd &A, VectorXd &x, VectorXd &b) override;
    };
}