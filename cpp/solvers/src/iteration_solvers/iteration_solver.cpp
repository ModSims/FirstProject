#include "iteration_solvers.h"
#include <iostream>
#include <vector>

using namespace Solvers;
using namespace Eigen;

IterationSolver::IterationSolver(SolverData *data, int max_iterations, double omega, Kernel::Timer *timer)
    : m_timer(timer), m_data(data), m_max_iterations(max_iterations), m_omega(omega)
{
    m_residuals = new VectorXd(m_max_iterations);
}

bool IterationSolver::forward()
{
    if (
        (getLastResidual() <= m_tolerance && m_timer->getCurrentTimeStep() > 0)
        || m_timer->isStopped() 
        || m_timer->getCurrentTimeStep() + 1 >= m_max_iterations
    )
    {
        m_residuals->conservativeResize(m_timer->getCurrentTimeStep());
        return false;
    }

    (*m_residuals)(m_timer->getCurrentTimeStep() + 1) = phi(m_data).norm();

    m_timer->update();
    return true;
}

VectorXd *IterationSolver::getResiduals()
{
    return m_residuals;
}

// Return the last residual by reference
double& IterationSolver::getLastResidual()
{
    return (*m_residuals)(m_timer->getCurrentTimeStep());
}

double IterationSolver::getTolerance()
{
    return m_tolerance;
}

double IterationSolver::getOmega()
{
    return m_omega;
}

int IterationSolver::getMaxIterations() const
{
    return m_max_iterations;
}
