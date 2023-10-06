#include "pybind_solvers.h"

using namespace Kernel;
using namespace Solvers;

void AddTrivialSolver(py::module &m) {
    py::class_<Trivial>(m, "Trivial")
        .def(py::init<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, Kernel::Timer*>())
        .def("phi", &Trivial::phi)
        .def("forward", &Trivial::forward)
        .def("setX", &Trivial::setX)
        .def("getResidualList", &Trivial::getResidualList)
        .def("getLastResidual", &Trivial::getLastResidual)
        .def("simulate", [](Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b, double dt) {
            Timer timer = Timer(dt);
            Trivial solver = Trivial(A, x, b, &timer);
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
