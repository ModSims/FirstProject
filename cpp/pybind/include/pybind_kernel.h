#pragma once

#include <pybind11/pybind11.h>
#include "kernel.h"

namespace py = pybind11;

void AddTimer(py::module &m);
