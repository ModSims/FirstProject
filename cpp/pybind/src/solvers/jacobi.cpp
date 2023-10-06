#include "pybind_solvers.h"

using namespace Kernel;
using namespace Solvers;

void AddJacobiSolver(py::module &m) {
    py::class_<Jacobi>(m, "Jacobi")
        .def(py::init<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, Kernel::Timer*>())
        .def("phi", &Jacobi::phi)
        .def("forward", &Jacobi::forward)
        .def("setX", &Jacobi::setX)
        .def("getResidualList", &Jacobi::getResidualList)
        .def("getLastResidual", &Jacobi::getLastResidual)
        .def("simulate", [](Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b, double dt) {
            Timer timer = Timer(dt);
            Jacobi solver = Jacobi(A, x, b, &timer);
            ProgressBar progressBar = ProgressBar(solver.getMaxIterations());
            timer.start();
            while (solver.forward()) {
                // Check for Python signals
                if (PyErr_CheckSignals() != 0) {
                    throw py::error_already_set();
                }
                progressBar.update();
            }
            timer.stop();
            return solver;
        });
}