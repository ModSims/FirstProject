#include "pybind_krylov_solvers.h"

using namespace Kernel;
using namespace Solvers;

void AddPreconditionedConjugateGradientSolver(py::module &m) {
    py::class_<PreconditionedConjugateGradient>(m, "PreconditionedConjugateGradient")
        .def(py::init<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, Kernel::Timer*>())
        .def("phi", &PreconditionedConjugateGradient::phi)
        .def("forward", &PreconditionedConjugateGradient::forward)
        .def("setX", &PreconditionedConjugateGradient::setX)
        .def("getAlphas", &PreconditionedConjugateGradient::getAlphas)
        .def("getLastAlpha", &PreconditionedConjugateGradient::getLastAlpha)
        .def("simulate", [](Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b, double dt) {
            Timer timer = Timer(dt);
            PreconditionedConjugateGradient solver = PreconditionedConjugateGradient(A, x, b, &timer);
            solver.setPreconditioner(Maths::Matrix::IncompleteCholeskyDecomposition(A).getP());
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