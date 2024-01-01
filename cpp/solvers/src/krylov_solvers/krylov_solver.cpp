#include "krylov_solvers.h"

#include <iostream>
#include <vector>

using namespace Solvers;
using namespace Eigen;

KrylovSolver::KrylovSolver(SolverData *data, Kernel::Timer *timer)
    : m_timer(timer), m_data(data)
{
    m_alphas = new VectorXd(100000);
}

void KrylovSolver::setBreakSolver(bool break_solver)
{
    m_break_solver = break_solver;
}

void KrylovSolver::innerLoop()
{
    while (true)
    {
        phi();

        if (m_timer->getCurrentTimeStep() + m_m < m_alphas->size()) {
            m_alphas->data()[m_timer->getCurrentTimeStep() + m_m] = m_current_alpha;
            m_timer->update();
        }

        if (m_current_alpha < m_tolerance)
        {
            m_break_solver = true;
            break;
        }
    }
}

bool KrylovSolver::forward()
{
    if (m_break_solver || m_m >= m_data->A.rows())
    {
        m_break_solver = true;
        m_alphas->conservativeResize(m_timer->getCurrentTimeStep());
        return false;
    }

    innerLoop();

    m_m += 1;
    return true;
}

VectorXd *KrylovSolver::getAlphas()
{
    return m_alphas;
}

double *KrylovSolver::getLastAlpha()
{
    return m_alphas->data() + m_timer->getCurrentTimeStep();
}

double KrylovSolver::getTolerance()
{
    return m_tolerance;
}

void KrylovSolver::setX(VectorXd x)
{
    m_data->x = x;
}

VectorXd *KrylovSolver::getX()
{
    return &m_data->x;
}
