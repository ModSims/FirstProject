#include "pybind_maths.h"

using namespace Maths;

void AddSpectral(py::module &m)
{
    py::module m_spectral = m.def_submodule("Spectral");
    m_spectral.def("getSpectralRadiusViaEigen", &Spectral::getSpectralRadiusViaEigen);
    m_spectral.def("getGelfandsSpectralApproximation", &Spectral::getGelfandsSpectralApproximation);
}