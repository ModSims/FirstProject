#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "iteration_solvers.h"
#include "krylov_solvers.h"

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

    class StaggeredGrid {
    public:
        StaggeredGrid(
            const int p_imax = 50,
            const int p_jmax = 50,
            const double p_xlength = 1.0,
            const double p_ylength = 1.0
        ) {
            imax = p_imax;
            jmax = p_jmax;
            xlength = p_xlength;
            ylength = p_ylength;
            p = MatrixXd::Zero(imax + 2, jmax + 2);
            po = MatrixXd::Zero(imax + 2, jmax + 2);
            RHS = MatrixXd::Zero(imax + 2, jmax + 2);
            res = MatrixXd::Zero(imax + 2, jmax + 2);
            u = MatrixXd::Zero(imax + 2, jmax + 3);
            F = MatrixXd::Zero(imax + 2, jmax + 3);
            v = MatrixXd::Zero(imax + 3, jmax + 2);
            G = MatrixXd::Zero(imax + 3, jmax + 2);
            u_interpolated = MatrixXd::Zero(imax + 2, jmax + 2);
            v_interpolated = MatrixXd::Zero(imax + 2, jmax + 2);
            flag_field = MatrixXi::Zero(imax + 2, jmax + 2);
            // Conjugated Gradient components
            Ares = MatrixXd::Zero(imax + 2, jmax + 2);
        }
        double dx() const { return xlength / imax; }
        double dy() const { return ylength / jmax; }
        double findMaxAbsoluteU() const;
        double findMaxAbsoluteV() const;
        void interpolateVelocity();
        int imax;
        int jmax;
        double xlength;
        double ylength;
        MatrixXd p;
        MatrixXd po;
        MatrixXd RHS;
        MatrixXd res;
        MatrixXd u;
        MatrixXd v;
        MatrixXd F;
        MatrixXd G;
        MatrixXd u_interpolated;
        MatrixXd v_interpolated;
        MatrixXi flag_field;
        // Conjugated Gradient components
        MatrixXd Ares;
    };

    class FluidSimulation {
        public:
            FluidSimulation(
                const int p_imax = 50,
                const int p_jmax = 50,
                const double p_xlength = 1.0,
                const double p_ylength = 1.0,
                const double p_t_end = 50.0,
                const double p_tau = 0.5,
                const double p_eps = 1e-3,
                const double p_omg = 1.7,
                const int p_itermax = 100,
                const double p_alpha = 0.9,
                const double p_Re = 100.0,
                double p_t = 0,
                double p_dt = 0.05,
                SolverType p_solver_type = SolverType::JACOBI
            ) {
                imax = p_imax;
                jmax = p_jmax;
                xlength = p_xlength;
                ylength = p_ylength;
                t_end = p_t_end;
                tau = p_tau;
                eps = p_eps;
                omg = p_omg;
                itermax = p_itermax;
                alpha = p_alpha;
                Re = p_Re;
                t = p_t;
                dt = p_dt;
                res_norm = 0.0;
                solver_type = p_solver_type;

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
            void selectDtAccordingToStabilityCondition();
            void computeF();
            void computeG();
            void computeRHS();
            void solveWithJacobi();
            void solveWithMultiGrid();
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
            void saveMatrices();

            virtual ~FluidSimulation() = default;
    };

    class LidDrivenCavity2D : public FluidSimulation {
        public:
            LidDrivenCavity2D(
                const int p_imax = 50,
                const int p_jmax = 50,
                const double p_xlength = 1.0,
                const double p_ylength = 1.0,
                const double p_t_end = 50.0,
                const double p_tau = 0.5,
                const double p_eps = 1e-3,
                const double p_omg = 1.7,
                const int p_itermax = 100,
                const double p_alpha = 0.9,
                const double p_Re = 100.0,
                double p_t = 0.0,
                double p_dt = 0.05,
                SolverType p_solver_type = SolverType::JACOBI
            ) : FluidSimulation(
                p_imax,
                p_jmax,
                p_xlength,
                p_ylength,
                p_t_end,
                p_tau,
                p_eps,
                p_omg,
                p_itermax,
                p_alpha,
                p_Re,
                p_t,
                p_dt,
                p_solver_type
            ) {
                grid = StaggeredGrid(imax, jmax, xlength, ylength);
            }
            void setBoundaryConditionsU() override;
            void setBoundaryConditionsV() override;
            void setBoundaryConditionsP() override;
            void setBoundaryConditionsVelocityGeometry() override {};
            void setBoundaryConditionsInterpolatedVelocityGeometry() override {};
            void setBoundaryConditionsPGeometry() override {};
    };

    class FlowOverStep2D : public FluidSimulation {
        public:
            FlowOverStep2D(
                const int p_imax = 50,
                const int p_jmax = 50,
                const double p_xlength = 1.0,
                const double p_ylength = 1.0,
                const double p_t_end = 50.0,
                const double p_tau = 0.5,
                const double p_eps = 1e-3,
                const double p_omg = 1.7,
                const int p_itermax = 100,
                const double p_alpha = 0.9,
                const double p_Re = 100.0,
                double p_t = 0.0,
                double p_dt = 0.05,
                SolverType p_solver_type = SolverType::JACOBI
            ) : FluidSimulation(
                p_imax,
                p_jmax,
                p_xlength,
                p_ylength,
                p_t_end,
                p_tau,
                p_eps,
                p_omg,
                p_itermax,
                p_alpha,
                p_Re,
                p_t,
                p_dt,
                p_solver_type
            ) {
                grid = StaggeredGrid(imax, jmax, xlength, ylength);
            }
            void setBoundaryConditionsU() override;
            void setBoundaryConditionsV() override;
            void setBoundaryConditionsP() override;
            void run();
    };

    class KarmanVortexStreet2D : public FluidSimulation {
        public:
            KarmanVortexStreet2D(
                const int p_imax = 50,
                const int p_jmax = 50,
                const double p_xlength = 1.0,
                const double p_ylength = 1.0,
                const double p_t_end = 50.0,
                const double p_tau = 0.5,
                const double p_eps = 1e-3,
                const double p_omg = 1.7,
                const int p_itermax = 100,
                const double p_alpha = 0.9,
                const double p_Re = 100.0,
                double p_t = 0.0,
                double p_dt = 0.05,
                SolverType p_solver_type = SolverType::JACOBI
            ) : FluidSimulation(
                p_imax,
                p_jmax,
                p_xlength,
                p_ylength,
                p_t_end,
                p_tau,
                p_eps,
                p_omg,
                p_itermax,
                p_alpha,
                p_Re,
                p_t,
                p_dt,
                p_solver_type
            ) {
                grid = StaggeredGrid(imax, jmax, xlength, ylength);
            }
            void setBoundaryConditionsU() override;
            void setBoundaryConditionsV() override;
            void setBoundaryConditionsP() override;
            void run();
    };

    class PlaneShearFlow2D : public FluidSimulation {
        public:
            PlaneShearFlow2D(
                const int p_imax = 50,
                const int p_jmax = 50,
                const double p_xlength = 1.0,
                const double p_ylength = 1.0,
                const double p_t_end = 50.0,
                const double p_tau = 0.5,
                const double p_eps = 1e-3,
                const double p_omg = 1.7,
                const int p_itermax = 100,
                const double p_alpha = 0.9,
                const double p_Re = 100.0,
                double p_t = 0.0,
                double p_dt = 0.05,
                SolverType p_solver_type = SolverType::JACOBI
            ) : FluidSimulation(
                p_imax,
                p_jmax,
                p_xlength,
                p_ylength,
                p_t_end,
                p_tau,
                p_eps,
                p_omg,
                p_itermax,
                p_alpha,
                p_Re,
                p_t,
                p_dt,
                p_solver_type
            ) {
                grid = StaggeredGrid(imax, jmax, xlength, ylength);
            }
            void setBoundaryConditionsU() override;
            void setBoundaryConditionsV() override;
            void setBoundaryConditionsP() override;
            void setBoundaryConditionsVelocityGeometry() override {};
            void setBoundaryConditionsInterpolatedVelocityGeometry() override {};
            void setBoundaryConditionsPGeometry() override {};
    };

    // Functions
    void saveVTK(FluidSimulation* sim);
}
