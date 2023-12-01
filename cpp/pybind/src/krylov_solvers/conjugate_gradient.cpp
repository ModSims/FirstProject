#include "pybind_krylov_solvers.h"

using namespace Kernel;
using namespace Solvers;

void AddConjugateGradientSolver(py::module &m) {
    py::class_<ConjugateGradient>(m, "ConjugateGradient")
        .def(py::init<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, Kernel::Timer*>())
        .def("phi", &ConjugateGradient::phi)
        .def("forward", &ConjugateGradient::forward)
        .def("setX", &ConjugateGradient::setX)
        .def("getAlphas", &ConjugateGradient::getAlphas)
        .def("getLastAlpha", &ConjugateGradient::getLastAlpha)
        .def("simulate", [](Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b, double dt) {
            Timer timer = Timer(dt);
            ConjugateGradient solver = ConjugateGradient(A, x, b, &timer);
            solver.prepareSolver();
            timer.start();
            while (solver.forward()) {
                // Check for Python signals
                if (PyErr_CheckSignals() != 0) {
                    solver.setBreakSolver(true);
                    throw py::error_already_set();
                }
            }
            timer.stop();
            return solver;
        });
}