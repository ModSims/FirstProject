#pragma once

#include <Eigen/Dense>

namespace CFD {
    using namespace Eigen;
    class Simulation {
        public:
            Simulation(
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
                double p_res = 99999
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
                res = p_res;
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
            double res;
    };
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
            u = MatrixXd::Zero(imax + 2, jmax + 3);
            F = MatrixXd::Zero(imax + 2, jmax + 3);
            v = MatrixXd::Zero(imax + 3, jmax + 2);
            G = MatrixXd::Zero(imax + 3, jmax + 2);
        }
        double dx() const { return xlength / imax; }
        double dy() const { return ylength / jmax; }
        double findMaxAbsoluteU() const;
        double findMaxAbsoluteV() const;
        int imax;
        int jmax;
        double xlength;
        double ylength;
        MatrixXd p;
        MatrixXd po;
        MatrixXd RHS;
        MatrixXd u;
        MatrixXd v;
        MatrixXd F;
        MatrixXd G;
    };
    void selectDtAccordingToStabilityCondition(StaggeredGrid *grid, Simulation *sim);
    void computeF(StaggeredGrid *grid, Simulation *sim);
    void computeG(StaggeredGrid *grid, Simulation *sim);
    void computeRHS(StaggeredGrid *grid, Simulation *sim);
    void updateStepLGLS(StaggeredGrid *grid, Simulation *sim);
    void computeResidual(StaggeredGrid *grid, Simulation *sim);
    void computeU(StaggeredGrid *grid, Simulation *sim);
    void computeV(StaggeredGrid *grid, Simulation *sim);
}
