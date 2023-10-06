#include "pybind_kernel.h"

namespace py = pybind11;

void AddTimer(py::module &m) {
    py::class_<Kernel::Timer>(m, "Timer")
        .def(py::init<double>())
        .def("start", &Kernel::Timer::start)
        .def("update", &Kernel::Timer::update)
        .def("stop", &Kernel::Timer::stop)
        .def("isStarted", &Kernel::Timer::isStarted)
        .def("isStopped", &Kernel::Timer::isStopped)
        .def("getCurrentTimeStep", &Kernel::Timer::getCurrentTimeStep)
        .def("getDurationInSeconds", &Kernel::Timer::getDurationInSeconds);
}