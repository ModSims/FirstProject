#include "solvers.h"

#include <iostream>
#include <vector>

using namespace Solvers;
using namespace Eigen;

double IterationSolver::m_tolerance = 1e-14;

IterationSolver::IterationSolver(MatrixXd A, VectorXd x, VectorXd b, int max_iterations, double omega, Kernel::Timer *timer)
    : m_timer(timer), m_A(A), m_x(x), m_b(b), m_max_iterations(max_iterations), m_omega(omega), m_progressbar(max_iterations)
{
    m_residuals = new VectorXd(m_max_iterations);
}

bool IterationSolver::forward()
{
    if (
        (*this->getLastResidual() <= m_tolerance && m_timer->getCurrentTimeStep() > 0)
        || m_timer->isStopped() 
        || m_timer->getCurrentTimeStep() + 1 >= m_max_iterations
    )
    {
        m_residuals->conservativeResize(m_timer->getCurrentTimeStep());
        return false;
    }

    this->setX(phi(m_A, m_x, m_b));
    double residual = (m_b - m_A * m_x).norm();

    int newTimeStep = m_timer->getCurrentTimeStep() + 1;
    (*m_residuals)(newTimeStep) = residual;

    m_timer->update();
    m_progressbar.update();
    return true;
}

VectorXd *IterationSolver::getResiduals()
{
    return m_residuals;
}

double *IterationSolver::getLastResidual()
{
    return m_residuals->data() + m_timer->getCurrentTimeStep();
}

double IterationSolver::getTolerance()
{
    return m_tolerance;
}

double IterationSolver::getOmega()
{
    return m_omega;
}

void IterationSolver::setX(VectorXd x)
{
    m_x = x;
}

int IterationSolver::getMaxIterations() const
{
    return m_max_iterations;
}