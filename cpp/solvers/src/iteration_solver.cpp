#include "solvers.h"

#include <iostream>
#include <vector>

using namespace Solvers;
using namespace Eigen;

double IterationSolver::m_tolerance = 1e-14;
int IterationSolver::m_max_iterations = 100000;

IterationSolver::IterationSolver(MatrixXd A, VectorXd x, VectorXd b, Kernel::Timer *timer)
    : m_timer(timer), m_A(A), m_x(x), m_b(b)
{
    m_residual_list = new VectorXd(m_max_iterations);
}

bool IterationSolver::forward()
{
    if (
        (*this->getLastResidual() <= m_tolerance && m_timer->getCurrentTimeStep() > 0)
        || m_timer->isStopped() 
        || m_timer->getCurrentTimeStep() + 1 >= m_max_iterations
    )
    {
        m_residual_list->conservativeResize(m_timer->getCurrentTimeStep());
        return false;
    }

    this->setX(phi(m_A, m_x, m_b));
    double residual = (m_b - m_A * m_x).norm();

    int newTimeStep = m_timer->getCurrentTimeStep() + 1;
    (*m_residual_list)(newTimeStep) = residual;

    m_timer->update();
    return true;
}

VectorXd *IterationSolver::getResidualList()
{
    return m_residual_list;
}

double *IterationSolver::getLastResidual()
{
    return m_residual_list->data() + m_timer->getCurrentTimeStep();
}

void IterationSolver::setX(VectorXd x)
{
    m_x = x;
}

int IterationSolver::getMaxIterations() const
{
    return m_max_iterations;
}