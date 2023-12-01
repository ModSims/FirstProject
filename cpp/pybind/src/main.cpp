#include "pybind_kernel.h"
#include "pybind_iteration_solvers.h"
#include "pybind_krylov_solvers.h"
#include "pybind_maths.h"

PYBIND11_MODULE(PyModSims, m) {
    m.doc() = "Python bindings for the C++ simulation code";
    // Kernel
    AddTimer(m);

    // Maths
    AddSpectral(m);
    AddMatrix(m);

    // Iteration solvers
    AddJacobiSolver(m);
    AddGaussSeidelSolver(m);
    AddTrivialSolver(m);
    AddRichardsonSolver(m);

    // Krylov solvers
    AddConjugateGradientSolver(m);
    AddPreconditionedConjugateGradientSolver(m);
}
