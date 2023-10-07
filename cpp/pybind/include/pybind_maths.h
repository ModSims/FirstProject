#pragma once

#include <pybind11/eigen.h>
#include "maths.h"

namespace py = pybind11;

void AddSpectral(py::module &m);
void AddMatrix(py::module &m);