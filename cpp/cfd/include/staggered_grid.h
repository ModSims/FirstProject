#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace CFD {
    using namespace Eigen;

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
            search_vector = MatrixXd::Zero(imax + 2, jmax + 2);
            Asearch_vector = MatrixXd::Zero(imax + 2, jmax + 2);

            dx = xlength / imax;
            dy = ylength / jmax;
            dx2 = dx * dx;
            dy2 = dy * dy;
            dxdy = dx * dy;
        }
        double findMaxAbsoluteU() const;
        double findMaxAbsoluteV() const;
        void interpolateVelocity();
        double dx;
        double dy;
        double dx2;
        double dy2;
        double dxdy;
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
        MatrixXd search_vector;
        MatrixXd Asearch_vector;
    };
}
