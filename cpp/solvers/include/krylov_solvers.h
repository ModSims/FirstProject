#pragma once
#include <Eigen/Dense>
#include "kernel.h"
#include "solvers.h"

namespace Solvers
{
    using namespace Eigen;
    class KrylovSolver
    {
    public:
        KrylovSolver(SolverData *data, Kernel::Timer *timer);
        KrylovSolver(KrylovSolver &solver) = delete;
        KrylovSolver &operator=(KrylovSolver &solver) = delete;
        KrylovSolver(KrylovSolver&& other) noexcept
            : m_timer(other.m_timer),
            m_alphas(other.m_alphas),
            m_data(other.m_data)
        {
            // Nullify other's pointers
            other.m_timer = nullptr;
            other.m_alphas = nullptr;
            other.m_data = nullptr;
        }
        virtual VectorXd phi() = 0;
        virtual void prepareSolver() = 0;
        void setBreakSolver(bool break_solver);
        bool forward();
        void innerLoop();
        void setX(VectorXd x);
        VectorXd *getAlphas();
        double *getLastAlpha();
        double getTolerance();
        VectorXd *getX();
        const double m_tolerance = 1e-30;

    protected:
        SolverData *m_data;
        VectorXd *m_alphas;
        VectorXd m_r;
        VectorXd m_p;
        VectorXd m_v;
        double m_current_alpha;
        double m_previous_alpha;
        double m_lambda;
        int m_m = 0;

        Kernel::Timer *m_timer;
        bool m_break_solver = false;
    };

    class ConjugateGradient : public KrylovSolver
    {
    public:
        ConjugateGradient(SolverData *data, Kernel::Timer *timer) : KrylovSolver(data, timer) {}
        VectorXd phi();
        void prepareSolver();
    };

    class PreconditionedConjugateGradient : public KrylovSolver
    {
    public:
        PreconditionedConjugateGradient(SolverData *data, Kernel::Timer *timer) : KrylovSolver(data, timer) {}
        VectorXd phi();
        void prepareSolver();
        void setPreconditioner(MatrixXd P);
    private:
        MatrixXd m_P;
        VectorXd m_z;
    };
}