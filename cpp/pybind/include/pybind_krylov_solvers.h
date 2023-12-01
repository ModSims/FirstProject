#pragma once

#include <pybind11/eigen.h>

#include "maths.h"
#include "krylov_solvers.h"

namespace py = pybind11;

void AddConjugateGradientSolver(py::module &m);
void AddPreconditionedConjugateGradientSolver(py::module &m);