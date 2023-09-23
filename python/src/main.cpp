#include <pybind11/eigen.h>
#include "solvers.h"

namespace py = pybind11;

PYBIND11_MODULE(ModSims, m) {

    py::class_<Solvers::Trivial>(m, "Trivial")
        .def(py::init<>())
        .def("solve", &Solvers::Trivial::solve);

    py::class_<Solvers::Jacobi>(m, "Jacobi")
        .def(py::init<>())
        .def("solve", &Solvers::Jacobi::solve);

    py::class_<Solvers::GaussSeidel>(m, "GaussSeidel")
        .def(py::init<>())
        .def("solve", &Solvers::GaussSeidel::solve);
}
