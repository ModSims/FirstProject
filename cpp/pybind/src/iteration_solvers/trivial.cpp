#include "pybind_iteration_solvers.h"

using namespace Kernel;
using namespace Solvers;

void AddTrivialSolver(py::module &m) {
    py::class_<Trivial>(m, "Trivial")
        .def(py::init<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, int, double, Kernel::Timer*>())
        .def("phi", &Trivial::phi)
        .def("forward", &Trivial::forward)
        .def("setX", &Trivial::setX)
        .def("getResiduals", &Trivial::getResiduals)
        .def("getLastResidual", &Trivial::getLastResidual)
        .def("simulate", [](Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b, int max_iterations, double omega, double dt) {
            Timer timer = Timer(dt);
            Trivial solver = Trivial(A, x, b, max_iterations, omega, &timer);
            timer.start();
            while (solver.forward()) {
                // Check for Python signals
                if (PyErr_CheckSignals() != 0) {
                    throw py::error_already_set();
                }
            }
            timer.stop();
            return solver;
        });
}
