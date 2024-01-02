#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "kernel.h"
#include "argparse.h"
#include "staggered_grid.h"
#include "multigrid.h"

namespace CFD {
    using namespace Eigen;

    enum FlagFieldMask {
        CELL_FLUID = 0b11111,
        FLUID_WEST = 0b00100,
        FLUID_EAST = 0b01000,
        FLUID_SOUTH = 0b00010,
        FLUID_NORTH = 0b00001,
        CELL_OBSTACLE = 0b00000,
        OBSTACLE_WEST = 0b11011,
        OBSTACLE_EAST = 0b11101,
        OBSTACLE_SOUTH = 0b11110,
        OBSTACLE_NORTH = 0b11100,
        MASK_CELL_TYPE = 0b10000,
        MASK_WEST = 0b00100,
        MASK_EAST = 0b01000,
        MASK_SOUTH = 0b00010,
        MASK_NORTH = 0b00001
    };

    enum SolverType {
        JACOBI,
        MULTIGRID_JACOBI,
        CONJUGATED_GRADIENT
    };
    SolverType convertSolverType(const std::string& solver);

    class FluidParams {
        public:
            FluidParams(const std::string name, int argc, char* argv[]);
            int imax = 100;
            int jmax = 100;
            double xlength = 1.0;
            double ylength = 1.0;
            double t_end = 5.0;
            double tau = 0.5;
            double eps = 1e-3;
            double omg = 1.7;
            int itermax = 100;
            double alpha = 0.9;
            double Re = 100.0;
            double t = 0;
            double dt = 0.05;
            double save_interval = 0.5;
            SolverType solver_type = SolverType::JACOBI;

            argparse::ArgumentParser argument_parser;
    };

    class FluidSimulation {
        public:
            FluidSimulation(const FluidParams& params) {
                imax = params.imax;
                jmax = params.jmax;
                xlength = params.xlength;
                ylength = params.ylength;
                t_end = params.t_end;
                tau = params.tau;
                eps = params.eps;
                omg = params.omg;
                itermax = params.itermax;
                alpha = params.alpha;
                Re = params.Re;
                t = params.t;
                dt = params.dt;
                solver_type = params.solver_type;
                save_interval = params.save_interval;
                res_norm = 0.0;
                multigrid_hierarchy = nullptr;
                res_norm_over_time = VectorXd::Zero(1e7);
            }
            int imax;
            int jmax;
            double xlength;
            double ylength;
            double t_end;
            double tau;
            double eps;
            double omg;
            int itermax;
            double alpha;
            double Re;
            double t;
            double dt;
            double res_norm;
            StaggeredGrid grid;
            Kernel::Timer timer;
            SolverType solver_type;
            double save_interval;
            VectorXd res_norm_over_time;

            // Multigrid components
            MultigridHierarchy *multigrid_hierarchy;

            void selectDtAccordingToStabilityCondition();
            void computeF();
            void computeG();
            void computeRHS();
            void solveWithJacobi();
            void solveWithMultigridJacobi();
            void solveWithConjugatedGradient();
            void computeResidual();
            void computeU();
            void computeV();
            void run();
            virtual void setBoundaryConditionsU() = 0;
            virtual void setBoundaryConditionsV() = 0;
            virtual void setBoundaryConditionsP() = 0;
            virtual void setBoundaryConditionsPGeometry();
            virtual void setBoundaryConditionsVelocityGeometry();
            virtual void setBoundaryConditionsInterpolatedVelocityGeometry();
            void saveData();

            virtual ~FluidSimulation() = default;
    };

    // Functions
    void saveVTK(FluidSimulation* sim);
}
