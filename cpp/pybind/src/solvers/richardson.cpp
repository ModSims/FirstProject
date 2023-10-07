#include "pybind_solvers.h"

using namespace Kernel;
using namespace Solvers;

void AddRichardsonSolver(py::module &m) {
    py::class_<Richardson>(m, "Richardson")
        .def(py::init<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, int, double, Kernel::Timer*>())
        .def("phi", &Richardson::phi)
        .def("forward", &Richardson::forward)
        .def("setX", &Richardson::setX)
        .def("getResiduals", &Richardson::getResiduals)
        .def("getLastResidual", &Richardson::getLastResidual)
        .def("simulate", [](Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b, int max_iterations, double omega, double dt) {
            Timer timer = Timer(dt);
            Richardson solver = Richardson(A, x, b, max_iterations, omega, &timer);
            ProgressBar progressBar = ProgressBar(solver.getMaxIterations());
            timer.start();
            while (solver.forward()) {
                // Check for Python signals
                if (PyErr_CheckSignals() != 0) {
                    throw py::error_already_set();
                }
                //progressBar.update();
            }
            timer.stop();
            return solver;
        });
}