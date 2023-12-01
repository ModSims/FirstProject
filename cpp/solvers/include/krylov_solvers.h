#pragma once
#include <Eigen/Dense>
#include "kernel.h"

namespace Solvers
{
    using namespace Eigen;
    class KrylovSolver
    {
    public:
        KrylovSolver(MatrixXd A, VectorXd x, VectorXd b, Kernel::Timer *timer);
        KrylovSolver(KrylovSolver &solver) = delete;
        KrylovSolver &operator=(KrylovSolver &solver) = delete;
        KrylovSolver(KrylovSolver&& other) noexcept
            : m_timer(other.m_timer),
            m_alphas(other.m_alphas),
            m_A(std::move(other.m_A)),
            m_x(std::move(other.m_x)),
            m_b(std::move(other.m_b))
        {
            // Nullify other's pointers
            other.m_timer = nullptr;
            other.m_alphas = nullptr;
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
        MatrixXd m_A;
        VectorXd m_x;
        VectorXd m_b;
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
        ConjugateGradient(MatrixXd A, VectorXd x, VectorXd b, Kernel::Timer *timer) : KrylovSolver(A, x, b, timer) {}
        VectorXd phi();
        void prepareSolver();
    };

    class PreconditionedConjugateGradient : public KrylovSolver
    {
    public:
        PreconditionedConjugateGradient(MatrixXd A, VectorXd x, VectorXd b, Kernel::Timer *timer) : KrylovSolver(A, x, b, timer) {}
        VectorXd phi();
        void prepareSolver();
        void setPreconditioner(MatrixXd P);
    private:
        MatrixXd m_P;
        VectorXd m_z;
    };
}