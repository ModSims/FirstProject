#pragma once
#include <Eigen/Dense>
#include "kernel.h"
#include "solvers.h"

namespace Solvers
{
    using namespace Eigen;
    class IterationSolver
    {
    public:
        IterationSolver(SolverData *data, int max_iterations, double omega, Kernel::Timer *timer);
        IterationSolver(IterationSolver &solver) = delete;
        IterationSolver &operator=(IterationSolver &solver) = delete;
        IterationSolver(IterationSolver&& other) noexcept
            : m_timer(other.m_timer),
            m_residuals(other.m_residuals),
            m_data(other.m_data),
            m_max_iterations(other.m_max_iterations),
            m_omega(other.m_omega)
        {
            // Nullify other's pointers
            other.m_timer = nullptr;
            other.m_residuals = nullptr;
            other.m_data = nullptr;
            other.m_max_iterations = 0;
            other.m_omega = 0;
        }
        virtual VectorXd phi(SolverData *data) = 0;
        virtual void prepareSolver() = 0;
        bool forward();
        VectorXd *getResiduals();
        double &getLastResidual();
        int getMaxIterations() const;
        double getTolerance();
        double getOmega();
        SolverData *getData() const { return m_data; }
        const double m_tolerance = 1e-14;

    private:
        Kernel::Timer *m_timer;
        VectorXd *m_residuals;
        SolverData *m_data;
        int m_max_iterations;
        double m_omega;
    };

    class Trivial : public IterationSolver
    {
    public:
        Trivial(SolverData *data, int max_iterations, double omega, Kernel::Timer *timer) : IterationSolver(data, max_iterations, omega, timer) {}
        VectorXd phi(SolverData *data) override;
        void prepareSolver();
    };

    class GaussSeidel : public IterationSolver
    {
    public:
        GaussSeidel(SolverData *data, int max_iterations, double omega, Kernel::Timer *timer) : IterationSolver(data, max_iterations, omega, timer) {}
        VectorXd phi(SolverData *data) override;
        void prepareSolver();
    };

    class Jacobi : public IterationSolver
    {
    public:
        Jacobi(SolverData *data, int max_iterations, double omega, Kernel::Timer *timer) : IterationSolver(data, max_iterations, omega, timer) {}
        VectorXd phi(SolverData *data) override;
        void prepareSolver();
    private:
        VectorXd m_D_inv;
    };

    class Richardson : public IterationSolver
    {
    public:
        Richardson(SolverData *data, int max_iterations, double omega, Kernel::Timer *timer) : IterationSolver(data, max_iterations, omega, timer) {}
        VectorXd phi(SolverData *data) override;
        void prepareSolver();
    };
}