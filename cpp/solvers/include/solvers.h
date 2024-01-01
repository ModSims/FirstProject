#pragma once
#include <Eigen/Dense>

namespace Solvers
{
    using namespace Eigen;
    class SolverData
    {
        public:
            SolverData() = default;
            SolverData(const SolverData &other) = default;
            SolverData(SolverData &&other) = default;
            SolverData &operator=(const SolverData &other) = default;
            SolverData &operator=(SolverData &&other) = default;
            virtual ~SolverData() = default;

            // Attributes
            MatrixXd A;
            VectorXd x;
            VectorXd b;
    };
}