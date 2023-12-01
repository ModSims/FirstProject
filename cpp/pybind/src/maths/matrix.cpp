#include "pybind_maths.h"

using namespace Maths;

void AddMatrix(py::module &m)
{
    py::module m_matrix = m.def_submodule("Matrix");
    m_matrix.def("getConditionNumber", &Maths::Matrix::getConditionNumber);
    m_matrix.def("isSymmetric", &Maths::Matrix::isSymmetric);
    m_matrix.def("isDiagonallyDominant", &Maths::Matrix::isDiagonallyDominant);
    m_matrix.def("isSPD", &Maths::Matrix::isSPD);
    py::class_<Maths::Matrix::IncompleteCholeskyDecomposition>(m_matrix, "IncompleteCholeskyDecomposition")
        .def(py::init<Eigen::MatrixXd&>())
        .def("getL", &Maths::Matrix::IncompleteCholeskyDecomposition::getL)
        .def("getF", &Maths::Matrix::IncompleteCholeskyDecomposition::getF)
        .def("getP", &Maths::Matrix::IncompleteCholeskyDecomposition::getP);
}