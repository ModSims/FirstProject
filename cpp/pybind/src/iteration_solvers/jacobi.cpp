#include "pybind_iteration_solvers.h"

using namespace Kernel;
using namespace Solvers;

void AddJacobiSolver(py::module &m) {
    py::class_<Jacobi>(m, "Jacobi")
        .def(py::init<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, int, double, Kernel::Timer*>())
        .def("phi", &Jacobi::phi)
        .def("forward", &Jacobi::forward)
        .def("setX", &Jacobi::setX)
        .def("getResiduals", &Jacobi::getResiduals)
        .def("getLastResidual", &Jacobi::getLastResidual)
        .def("simulate", [](Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b, int max_iterations, double omega, double dt) {
            Timer timer = Timer(dt);
            Jacobi solver = Jacobi(A, x, b, max_iterations, omega, &timer);
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