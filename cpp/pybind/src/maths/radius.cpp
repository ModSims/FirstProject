#include "pybind_maths.h"

using namespace Maths;

void AddRadians(py::module &m)
{
    py::class_<Maths::Radius>(m, "Radius")
        .def(py::init<>())
        .def("calculateSpectralRadius", &Radius::calculateSpectralRadius)
        .def("gelfands_spectral_approximation", &Radius::gelfands_spectral_approximation);
}