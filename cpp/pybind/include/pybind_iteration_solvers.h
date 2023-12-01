#pragma once

#include <pybind11/eigen.h>
#include "iteration_solvers.h"

namespace py = pybind11;

void AddTrivialSolver(py::module &m);
void AddJacobiSolver(py::module &m);
void AddGaussSeidelSolver(py::module &m);
void AddRichardsonSolver(py::module &m);
