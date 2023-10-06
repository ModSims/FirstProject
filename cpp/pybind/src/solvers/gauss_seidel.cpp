#include "pybind_solvers.h"

using namespace Kernel;
using namespace Solvers;

void AddGaussSeidelSolver(py::module &m) {
    py::class_<GaussSeidel>(m, "GaussSeidel")
        .def(py::init<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, Kernel::Timer*>())
        .def("phi", &GaussSeidel::phi)
        .def("forward", &GaussSeidel::forward)
        .def("setX", &GaussSeidel::setX)
        .def("getResidualList", &GaussSeidel::getResidualList)
        .def("getLastResidual", &GaussSeidel::getLastResidual)
        .def("simulate", [](Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b, double dt) {
            Timer timer = Timer(dt);
            GaussSeidel solver = GaussSeidel(A, x, b, &timer);
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