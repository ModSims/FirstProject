 #include "krylov_solvers.h"
#include "maths.h"

using namespace Solvers;
using namespace Eigen;

void ConjugateGradient::prepareSolver() {
    m_r = m_b - m_A * m_x;
    m_p = m_r;
    m_current_alpha = m_r.dot(m_r);
}

VectorXd ConjugateGradient::phi() {
    m_v = m_A * m_p;
    m_lambda = (m_current_alpha / m_p.dot(m_v));
    m_x = m_x + m_lambda * m_p;
    m_r = m_r - m_lambda * m_v;
    m_previous_alpha = m_current_alpha;
    m_current_alpha = m_r.dot(m_r);
    m_p = m_r + (m_current_alpha / m_previous_alpha) * m_p;

    return m_x;
}

 