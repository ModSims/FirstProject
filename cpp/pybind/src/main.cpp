#include "pybind_kernel.h"
#include "pybind_solvers.h"
#include "pybind_maths.h"

PYBIND11_MODULE(PyModSims, m) {
    m.doc() = "Python bindings for the C++ simulation code";
    AddTimer(m);
    AddJacobiSolver(m);
    AddGaussSeidelSolver(m);
    AddTrivialSolver(m);
    AddRichardsonSolver(m);
    AddSpectral(m);
    AddMatrix(m);
}
