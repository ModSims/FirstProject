#include "krylov_solvers.h"
#include "maths.h"
#include <iostream>

using namespace Solvers;
using namespace Eigen;

void PreconditionedConjugateGradient::setPreconditioner(MatrixXd P) {
    m_P = P;
}

void PreconditionedConjugateGradient::prepareSolver() {
    // break if m_P is not set
    if (m_P.rows() == 0 || m_P.cols() == 0) {
        throw std::runtime_error("Preconditioner not set!");
    }
    m_r = m_data->b - m_data->A * m_data->x;
    m_z = m_P * m_r;
    m_p = m_z;
    m_current_alpha = m_r.dot(m_z);
}

VectorXd PreconditionedConjugateGradient::phi() {
    m_v = m_data->A * m_p;
    m_lambda = (m_current_alpha / m_p.dot(m_v));
    m_data->x = m_data->x + m_lambda * m_p;
    m_r = m_r - m_lambda * m_v;
    m_z = m_P * m_r;
    m_previous_alpha = m_current_alpha;
    m_current_alpha = m_r.dot(m_z);
    m_p = m_z + (m_current_alpha / m_previous_alpha) * m_p;

    return m_data->x;
}

 