#include <iostream>
#include <chrono>
#include "cfd.h"

namespace CFD {
    namespace ME_X {
        // Slide 15
        double uu_x(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (
                    (1/grid->dx)*(pow((grid->u(i,j) + grid->u(i+1,j))/2,2) - pow((grid->u(i-1,j)+grid->u(i,j))/2,2))
                )
                +
                (
                    (sim->alpha/grid->dx)*(std::abs(grid->u(i,j) + grid->u(i+1,j))*(grid->u(i,j) - grid->u(i+1,j))/4 - std::abs(grid->u(i-1,j) + grid->u(i,j))*(grid->u(i-1,j) - grid->u(i,j))/4)
                )
            );
        }

        double uv_y(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (
                    (1/grid->dy)*((grid->v(i,j) + grid->v(i+1,j))*(grid->u(i,j) + grid->u(i,j+1))/4 - (grid->v(i,j-1) + grid->v(i+1,j-1))*(grid->u(i,j-1) + grid->u(i,j))/4)
                )
                +
                (
                    (sim->alpha/grid->dy)*(std::abs(grid->v(i,j) + grid->v(i+1,j))*(grid->u(i,j) - grid->u(i,j+1))/4 - std::abs(grid->v(i,j-1) + grid->v(i+1,j-1))*(grid->u(i,j-1) - grid->u(i,j))/4)
                )
            );
        }

        double uu_xx(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (grid->u(i+1,j) - 2*grid->u(i,j) + grid->u(i-1,j))/pow(grid->dx, 2)
            );
        }

        double uu_yy(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (grid->u(i,j+1) - 2*grid->u(i,j) + grid->u(i,j-1))/pow(grid->dy, 2)
            );
        }

        double p_x(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (grid->p(i+1,j) - grid->p(i,j))/grid->dx
            );
        }
    }

    namespace ME_Y {
        double uv_x(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (
                    (1/grid->dx)*((grid->u(i,j) + grid->u(i,j+1))*(grid->v(i,j) + grid->v(i+1,j))/4 - (grid->u(i-1,j) + grid->u(i-1,j+1))*(grid->v(i-1,j) + grid->v(i,j))/4)
                )
                +
                (
                    (sim->alpha/grid->dx)*(std::abs(grid->u(i,j) + grid->u(i,j+1))*(grid->v(i,j) - grid->v(i+1,j))/4 - std::abs(grid->u(i-1,j) + grid->u(i-1,j+1))*(grid->v(i-1,j) - grid->v(i,j))/4)
                )
            );
        }

        double vv_y(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (
                    (1/grid->dy)*(pow((grid->v(i,j) + grid->v(i,j+1))/2,2) - pow((grid->v(i,j-1)+grid->v(i,j))/2,2))
                )
                +
                (
                    (sim->alpha/grid->dy)*(std::abs(grid->v(i,j) + grid->v(i,j+1))*(grid->v(i,j) - grid->v(i,j+1))/4 - std::abs(grid->v(i,j-1) + grid->v(i,j))*(grid->v(i,j-1) - grid->v(i,j))/4)
                )
            );
        }

        double vv_xx(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (grid->v(i+1,j) - 2*grid->v(i,j) + grid->v(i-1,j))/pow(grid->dx, 2)
            );
        }

        double vv_yy(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (grid->v(i,j+1) - 2*grid->v(i,j) + grid->v(i,j-1))/pow(grid->dy, 2)
            );
        }

        double p_y(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (grid->p(i,j+1) - grid->p(i,j))/grid->dy
            );
        }
    }

    namespace CE {
        double u_x(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (grid->u(i,j) - grid->u(i-1,j))/(grid->dx)
            );
        }

        double v_y(StaggeredGrid *grid, FluidSimulation *sim, int i, int j) {
            return (
                (grid->v(i,j) - grid->v(i,j-1))/(grid->dy)
            );
        }
    }

    double StaggeredGrid::findMaxAbsoluteU() const {
        double max = 0.0;
        for (int i = 1; i < imax + 1; i++) {
            for (int j = 1; j < jmax + 2; j++) {
                if (std::abs(u(i, j)) > max) {
                    max = std::abs(u(i, j));
                }
            }
        }
        return max;
    }
    double StaggeredGrid::findMaxAbsoluteV() const {
        double max = 0.0;
        for (int i = 1; i < imax + 2; i++) {
            for (int j = 1; j < jmax + 1; j++) {
                if (std::abs(v(i, j)) > max) {
                    max = std::abs(v(i, j));
                }
            }
        }
        return max;
    }
    void StaggeredGrid::interpolateVelocity() {
        for (int i = 0; i < imax + 2; i++) {
            for (int j = 0; j < jmax + 2; j++) {
                u_interpolated(i, j) = (u(i, j) + u(i, j+1)) / 2;
                v_interpolated(i, j) = (v(i, j) + v(i+1, j)) / 2;
            }
        }
    }
    void FluidSimulation::selectDtAccordingToStabilityCondition() {
        double left = (this->Re/2) * pow(((1/pow(this->grid.dx, 2))+(1/pow(this->grid.dy, 2))), -1);
        double middle = this->grid.dx / this->grid.findMaxAbsoluteU();
        double right = this->grid.dy / this->grid.findMaxAbsoluteV();
        this->dt = this->tau * std::min(std::min(left, middle), right);
    }

    void FluidSimulation::computeF() {
        // Boundary conditions
        for (int j = 0; j < this->grid.jmax + 3; j++) {
            this->grid.F(0, j) = this->grid.u(0, j);
            this->grid.F(this->grid.imax, j) = this->grid.u(this->grid.imax, j);
        }

        for (int i = 1; i < this->grid.imax + 1; i++) {
            for (int j = 1; j < this->grid.jmax + 2; j++) {
                this->grid.F(i,j) = this->grid.u(i,j) + this->dt * (
                    (1/this->Re) * (ME_X::uu_xx(&this->grid, this, i, j) + ME_X::uu_yy(&this->grid, this, i, j)) 
                    - 
                    ME_X::uu_x(&this->grid, this, i, j)
                    -
                    ME_X::uv_y(&this->grid, this, i, j)
                );
            }
        }
    }

    void FluidSimulation::computeG() {
        // Boundary conditions
        for (int i = 0; i < this->grid.imax + 3; i++) {
            this->grid.G(i, 0) = this->grid.v(i, 0);
            this->grid.G(i, this->grid.jmax) = this->grid.v(i, this->grid.jmax);
        }

        for (int i = 1; i < this->grid.imax + 2; i++) {
            for (int j = 1; j < this->grid.jmax + 1; j++) {
                this->grid.G(i,j) = this->grid.v(i,j) + this->dt * (
                    (1/this->Re) * (ME_Y::vv_xx(&this->grid, this, i, j) + ME_Y::vv_yy(&this->grid, this, i, j)) 
                    - 
                    ME_Y::uv_x(&this->grid, this, i, j)
                    -
                    ME_Y::vv_y(&this->grid, this, i, j)
                );
            }
        }
    }

    void FluidSimulation::computeRHS() {
        for (int i = 1; i < this->grid.imax + 1; i++) {
            for (int j = 1; j < this->grid.jmax + 1; j++) {
                this->grid.RHS(i,j) = (1/this->dt)*((this->grid.F(i,j) - this->grid.F(i-1,j))/this->grid.dx+(this->grid.G(i,j) - this->grid.G(i,j-1))/this->grid.dy);
            }
        }
    }

    void FluidSimulation::computeResidual() {
        this->grid.res = this->grid.p - this->grid.po;
        this->res_norm = this->grid.res.norm();
    }

    void FluidSimulation::computeU() {
        for (int i = 1; i < this->grid.imax + 1; i++) {
            for (int j = 1; j < this->grid.jmax + 2; j++) {
                this->grid.u(i,j) = this->grid.F(i,j) - (this->dt/this->grid.dx) * (this->grid.p(i+1,j) - this->grid.p(i,j));
            }
        }
    }

    void FluidSimulation::computeV() {
        for (int i = 1; i < this->grid.imax + 2; i++) {
            for (int j = 1; j < this->grid.jmax + 1; j++) {
                this->grid.v(i,j) = this->grid.G(i,j) - (this->dt/this->grid.dy) * (this->grid.p(i,j+1) - this->grid.p(i,j));
            }
        }
    }

    void FluidSimulation::setBoundaryConditionsInterpolatedVelocityGeometry() {
        // Geometry boundaries
        for (int i = 1; i < this->grid.imax + 1; i++) {
            for (int j = 1; j < this->grid.jmax + 1; j++) {
                if ((this->grid.flag_field(i, j) & FlagFieldMask::MASK_CELL_TYPE) == (FlagFieldMask::CELL_OBSTACLE & FlagFieldMask::MASK_CELL_TYPE)) {
                    this->grid.u_interpolated(i,j) = 0.0;
                    this->grid.v_interpolated(i,j) = 0.0;
                }
            }
        }
    }

    void FluidSimulation::setBoundaryConditionsPGeometry() {
        // Geometry boundaries
        double tmp_p = 0.0;
        int counter = 0;
        for (int i = 1; i < this->grid.imax + 1; i++) {
            for (int j = 1; j < this->grid.jmax + 1; j++) {
                tmp_p = 0.0;
                counter = 0;
                // check if is obstacle
                if ((this->grid.flag_field(i, j) & FlagFieldMask::MASK_CELL_TYPE) == (FlagFieldMask::CELL_OBSTACLE & FlagFieldMask::MASK_CELL_TYPE)) {
                    // obstacle cell
                    if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_NORTH) {
                        tmp_p += this->grid.p(i, j+1);
                        counter++;
                    }
                    if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_SOUTH) {
                        tmp_p += this->grid.p(i, j-1);
                        counter++;
                    }
                    if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_WEST) {
                        tmp_p += this->grid.p(i-1, j);
                        counter++;
                    }
                    if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_EAST) {
                        tmp_p += this->grid.p(i+1, j);
                        counter++;
                    }
                    if (counter > 0) {
                        this->grid.p(i, j) = tmp_p / counter;
                    }
                }
            }
        }
    }

    void FluidSimulation::setBoundaryConditionsVelocityGeometry() {
        // Geometry boundaries
        for (int i = 1; i < this->grid.imax + 1; i++) {
            for (int j = 1; j < this->grid.jmax + 1; j++) {
                // check if is obstacle
                if ((this->grid.flag_field(i, j) & FlagFieldMask::MASK_CELL_TYPE) == (FlagFieldMask::CELL_OBSTACLE & FlagFieldMask::MASK_CELL_TYPE)) {
                    // corner cell
                    if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_NORTH && this->grid.flag_field(i, j) & FlagFieldMask::FLUID_EAST) {
                        this->grid.u(i, j) = 0.0;
                        this->grid.v(i, j) = 0.0;
                        this->grid.u(i-1, j) = -this->grid.u(i-1, j+1);
                        this->grid.v(i, j-1) = -this->grid.v(i+1, j-1);
                        this->grid.G(i, j) = this->grid.v(i, j);
                        this->grid.F(i, j) = this->grid.u(i, j);
                    }
                    else if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_NORTH && this->grid.flag_field(i, j) & FlagFieldMask::FLUID_WEST) {
                        this->grid.u(i, j) = 0.0;
                        this->grid.v(i, j) = 0.0;
                        this->grid.u(i-1, j) = -this->grid.u(i-1, j+1);
                        this->grid.v(i, j-1) = -this->grid.v(i-1, j-1);
                        this->grid.G(i, j) = this->grid.v(i, j);
                        this->grid.F(i, j) = this->grid.u(i, j);
                    }
                    else if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_SOUTH && this->grid.flag_field(i, j) & FlagFieldMask::FLUID_EAST) {
                        this->grid.u(i, j) = 0.0;
                        this->grid.v(i, j) = 0.0;
                        this->grid.u(i-1, j) = -this->grid.u(i-1, j-1);
                        this->grid.v(i, j-1) = -this->grid.v(i+1, j-1);
                        this->grid.G(i, j) = this->grid.v(i, j);
                        this->grid.F(i, j) = this->grid.u(i, j);
                    }
                    else if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_SOUTH && this->grid.flag_field(i, j) & FlagFieldMask::FLUID_WEST) {
                        this->grid.u(i, j) = 0.0;
                        this->grid.v(i, j) = 0.0;
                        this->grid.u(i-1, j) = -this->grid.u(i-1, j-1);
                        this->grid.v(i, j-1) = -this->grid.v(i-1, j-1);
                        this->grid.G(i, j) = this->grid.v(i, j);
                        this->grid.F(i, j) = this->grid.u(i, j);
                    }
                    // obstacle cell
                    else if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_NORTH) {
                        this->grid.u(i, j) = -this->grid.u(i, j+1);
                        this->grid.u(i-1, j) = -this->grid.u(i-1, j+1);
                        this->grid.v(i, j) = 0.0;
                        this->grid.G(i, j) = this->grid.v(i, j);
                    }
                    else if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_SOUTH) {
                        this->grid.u(i, j) = -this->grid.u(i, j-1);
                        this->grid.u(i-1, j) = -this->grid.u(i-1, j-1);
                        this->grid.v(i, j) = 0.0;
                        this->grid.G(i, j) = this->grid.v(i, j);
                    }
                    else if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_WEST) {
                        this->grid.v(i, j-1) = -this->grid.v(i-1, j-1);
                        this->grid.v(i, j) = -this->grid.v(i-1, j);
                        this->grid.u(i-1, j) = 0.0;
                        this->grid.F(i-1, j) = this->grid.u(i-1, j);
                    }
                    else if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_EAST) {
                        this->grid.v(i, j-1) = -this->grid.v(i+1, j-1);
                        this->grid.v(i, j) = -this->grid.v(i+1, j);
                        this->grid.u(i+1, j) = 0.0;
                        this->grid.F(i+1, j) = this->grid.u(i+1, j);
                    }
                    else if (this->grid.flag_field(i, j) == FlagFieldMask::CELL_OBSTACLE) {
                        // interior obstacle cell, so no-slip
                        this->grid.u(i,j) = 0.0;
                        this->grid.v(i,j) = 0.0;
                    }
                }
            }
        }
    }

    void FluidSimulation::solveWithJacobi() {
        // reset norm check
        this->res_norm = 0.0;
        int n = 0;

        while ((this->res_norm > this->eps || this->res_norm == 0) && n < this->itermax) {
            this->setBoundaryConditionsP();
            this->setBoundaryConditionsPGeometry();
            this->grid.po = this->grid.p;

            // Jacobi smoother with relaxation factor (omega)
            for (int i = 1; i <= this->grid.imax; i++) {
                for (int j = 1; j <= this->grid.jmax; j++) {
                    this->grid.p(i, j) = (1 - this->omg) * this->grid.p(i, j) +
                                        this->omg * 0.25 * (
                                            this->grid.p(i - 1, j) + this->grid.p(i + 1, j) +
                                            this->grid.p(i, j - 1) + this->grid.p(i, j + 1)
                                            - this->grid.dxdy * (this->grid.RHS(i, j))
                                        ) / (1 + 2 * (this->grid.dx2 + this->grid.dy2));
                }
            }

            this->computeResidual();
            this->res_norm_over_it_with_pressure_solver(this->it) = this->res_norm;
            this->it++;
            n++;
        }
    }


    void FluidSimulation::solveWithMultigridJacobi() {
        // reset norm check
        this->res_norm = 0.0;
        int n = 0;

        while ((this->res_norm > this->eps || this->res_norm == 0) && n < this->itermax) {
            this->setBoundaryConditionsP();
            this->setBoundaryConditionsPGeometry();
            this->grid.po = this->grid.p;

            Multigrid::vcycle(this->multigrid_hierarchy, this->multigrid_hierarchy->numLevels() - 1, this->omg);

            this->computeResidual();
            this->res_norm_over_it_with_pressure_solver(this->it) = this->res_norm;
            this->it++;
            n++;
        }
    }

    void FluidSimulation::solveWithConjugatedGradient() {
        this->setBoundaryConditionsP();
        this->setBoundaryConditionsPGeometry();

        // reset norm check
        this->res_norm = 0.0;
        int n = 0;
        double lambda = 0.0;
        double alpha_top = 0.0;
        double alpha_bottom = 0.0;
        double alpha_top_new = 0.0;
        int maxiterations = std::max(this->grid.imax, this->grid.jmax);

        // Initial residual vector of Ax=b
        for (int i = 1; i < this->grid.imax + 1; i++) {
            for (int j = 1; j < this->grid.jmax + 1; j++) {
                this->grid.res(i,j) = this->grid.RHS(i,j) - (
                    (1/this->grid.dx2)*(this->grid.p(i+1,j) - 2*this->grid.p(i,j) + this->grid.p(i-1,j)) +
                    (1/this->grid.dy2)*(this->grid.p(i,j+1) - 2*this->grid.p(i,j) + this->grid.p(i,j-1))
                );
                // copy residual to search_vector
                this->grid.search_vector(i,j) = this->grid.res(i,j);
                alpha_top += this->grid.res(i,j)*this->grid.res(i,j);
            }
        }

        while ((this->res_norm > this->eps || this->res_norm == 0) && n < maxiterations) {
            this->setBoundaryConditionsP();
            this->setBoundaryConditionsPGeometry();
            this->grid.po = this->grid.p;

            alpha_bottom = 0.0;
            // Laplacian operator of grid.res, because of dot product of <res, Asearch_vector>, A-Matrix is the laplacian operator
            for (int i = 1; i < this->grid.imax + 1; i++) {
                for (int j = 1; j < this->grid.jmax + 1; j++) {
                    this->grid.Asearch_vector(i,j) = (
                        (1/this->grid.dx2)*(this->grid.search_vector(i+1,j) - 2*this->grid.search_vector(i,j) + this->grid.search_vector(i-1,j)) +
                        (1/this->grid.dy2)*(this->grid.search_vector(i,j+1) - 2*this->grid.search_vector(i,j) + this->grid.search_vector(i,j-1))
                    );
                    alpha_bottom += this->grid.search_vector(i,j)*this->grid.Asearch_vector(i,j);
                }
            }
            // Update pressure and new residual
            lambda = alpha_top/alpha_bottom;
            alpha_top_new = 0.0;
            for (int i = 1; i < this->grid.imax + 1; i++) {
                for (int j = 1; j < this->grid.jmax + 1; j++) {
                    this->grid.p(i,j) += lambda*this->grid.search_vector(i,j);
                    this->grid.res(i,j) -= lambda*this->grid.Asearch_vector(i,j);
                    alpha_top_new += this->grid.res(i,j)*this->grid.res(i,j);
                }
            }

            // Calculate new search vector
            lambda = alpha_top_new/alpha_top;
            for (int i = 1; i < this->grid.imax + 1; i++) {
                for (int j = 1; j < this->grid.jmax + 1; j++) {
                    this->grid.search_vector(i,j) = this->grid.res(i,j) + lambda*this->grid.search_vector(i,j);
                }
            }
            alpha_top = alpha_top_new;

            this->res_norm = (this->grid.po - this->grid.p).norm();
            this->res_norm_over_it_with_pressure_solver(this->it) = this->res_norm;
            this->it++;
            n++;
        }
    }

    void FluidSimulation::solveWithMultigridPCG() {
        this->setBoundaryConditionsP();
        this->setBoundaryConditionsPGeometry();

        // reset norm check
        this->res_norm = 0.0;
        int n = 0;
        double alpha = 0.0;
        double alpha_top = 0.0;
        double alpha_bottom = 0.0;
        double beta = 0.0;
        double beta_top = 0.0;
        double beta_bottom = 0.0;

        int maxiterations = std::max(this->grid.imax, this->grid.jmax);

        // Initial residual vector of Ax=b
        for (int i = 1; i < this->grid.imax + 1; i++) {
            for (int j = 1; j < this->grid.jmax + 1; j++) {
                this->grid.res(i,j) = this->grid.RHS(i,j) - (
                    (1/this->grid.dx2)*(this->grid.p(i+1,j) - 2*this->grid.p(i,j) + this->grid.p(i-1,j)) +
                    (1/this->grid.dy2)*(this->grid.p(i,j+1) - 2*this->grid.p(i,j) + this->grid.p(i,j-1))
                );
                // copy residual to preconditioner and reset grids
                this->preconditioner.RHS(i,j) = this->grid.res(i,j);
                this->preconditioner.p(i,j) = 0.0;
            }
        }

        // Initial guess for error vector
        Multigrid::vcycle(this->multigrid_hierarchy, this->multigrid_hierarchy->numLevels() - 1, this->omg);

        // Initial search vector
        this->grid.search_vector = this->preconditioner.p;

        while ((this->res_norm > this->eps || this->res_norm == 0) && n < maxiterations) {
            this->setBoundaryConditionsP();
            this->setBoundaryConditionsPGeometry();
            this->grid.po = this->grid.p;

            alpha_top = 0.0;
            alpha_bottom = 0.0;

            // Calculate alpha
            // Laplacian operator of error_vector from multigrid, because of dot product of <A, Pi>, A-Matrix is the laplacian operator
            for (int i = 1; i < this->grid.imax + 1; i++) {
                for (int j = 1; j < this->grid.jmax + 1; j++) {
                    this->grid.Asearch_vector(i,j) = (
                        (1/this->grid.dx2)*(this->grid.search_vector(i+1,j) - 2*this->grid.search_vector(i,j) + this->grid.search_vector(i-1,j)) +
                        (1/this->grid.dy2)*(this->grid.search_vector(i,j+1) - 2*this->grid.search_vector(i,j) + this->grid.search_vector(i,j-1))
                    );
                    alpha_top += this->preconditioner.p(i,j)*this->grid.res(i,j);
                    alpha_bottom += this->grid.search_vector(i,j)*this->grid.Asearch_vector(i,j);
                }
            }
            alpha = alpha_top/alpha_bottom;

            // Update pressure and new residual
            for (int i = 1; i < this->grid.imax + 1; i++) {
                for (int j = 1; j < this->grid.jmax + 1; j++) {
                    this->grid.p(i,j) += alpha*this->grid.search_vector(i,j);
                    this->grid.res(i,j) -= alpha*this->grid.Asearch_vector(i,j);
                    this->preconditioner.RHS(i,j) = this->grid.res(i,j);
                }
            }
            
            // New guess for error vector
            Multigrid::vcycle(this->multigrid_hierarchy, this->multigrid_hierarchy->numLevels() - 1, this->omg);

            // Calculate beta
            for (int i = 1; i < this->grid.imax + 1; i++) {
                for (int j = 1; j < this->grid.jmax + 1; j++) {
                    beta_top += this->preconditioner.p(i,j)*this->grid.res(i,j);
                }
            }
            beta = beta_top/alpha_top;

            // Calculate new search vector
            for (int i = 1; i < this->grid.imax + 1; i++) {
                for (int j = 1; j < this->grid.jmax + 1; j++) {
                    this->grid.search_vector(i,j) = this->preconditioner.p(i,j) + beta*this->grid.search_vector(i,j);
                }
            }

            this->res_norm = (this->grid.po - this->grid.p).norm();
            this->res_norm_over_it_with_pressure_solver(this->it) = this->res_norm;
            this->it++;
            n++;
        }
    }

    void FluidSimulation::run() {
        double last_saved = 0.0;

        // Function pointer to solver
        void (CFD::FluidSimulation::*pressure_solver)();

        if (this->solver_type == SolverType::JACOBI) {
            pressure_solver = &FluidSimulation::solveWithJacobi;
            std::cout << "Solver: Jacobi" << std::endl;
        }
        else if (this->solver_type == SolverType::MULTIGRID_JACOBI) {
            pressure_solver = &FluidSimulation::solveWithMultigridJacobi;

            // check if imax and jmax are powers of 2, if not throw exception
            if ((this->grid.imax & (this->grid.imax - 1)) != 0 || (this->grid.jmax & (this->grid.jmax - 1)) != 0) {
                throw std::invalid_argument("imax and jmax must be powers of 2");
            }

            int imax_levels = std::log2(this->grid.imax);
            int jmax_levels = std::log2(this->grid.jmax);
            int levels = std::min(imax_levels, jmax_levels);

            this->multigrid_hierarchy = new MultigridHierarchy(levels, &this->grid);

            std::cout << "Solver: Multigrid Jacobi" << std::endl;
        }
        else if (this->solver_type == SolverType::CONJUGATED_GRADIENT) {
            pressure_solver = &FluidSimulation::solveWithConjugatedGradient;
            std::cout << "Solver: Conjugated Gradient" << std::endl;
        }
        else if (this->solver_type == SolverType::MULTIGRID_PCG) {
            pressure_solver = &FluidSimulation::solveWithMultigridPCG;

            // check if imax and jmax are powers of 2, if not throw exception
            if ((this->grid.imax & (this->grid.imax - 1)) != 0 || (this->grid.jmax & (this->grid.jmax - 1)) != 0) {
                throw std::invalid_argument("imax and jmax must be powers of 2");
            }

            int imax_levels = std::log2(this->grid.imax);
            int jmax_levels = std::log2(this->grid.jmax);
            int levels = std::min(imax_levels, jmax_levels);

            this->preconditioner = StaggeredGrid(this->grid.imax, this->grid.jmax, this->grid.xlength, this->grid.ylength);

            this->multigrid_hierarchy = new MultigridHierarchy(levels, &this->preconditioner);

            std::cout << "Solver: Multigrid PCG" << std::endl;
        }
        else {
            throw std::invalid_argument("Invalid solver type");
        }


        while(this->t < this->t_end) {
            this->lastTimestamp = static_cast<int>(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count());
            this->selectDtAccordingToStabilityCondition();
            this->setBoundaryConditionsU();
            this->setBoundaryConditionsV();
            this->setBoundaryConditionsVelocityGeometry();
            this->computeF();
            this->computeG();
            this->setBoundaryConditionsVelocityGeometry();
            this->computeRHS();
                
            (this->*pressure_solver)();

            this->computeResidual();
            this->res_norm_over_it_without_pressure_solver(this->it_wo_pressure_solver) = this->res_norm;

            this->setBoundaryConditionsP();
            this->computeU();
            this->computeV();
            this->setBoundaryConditionsU();
            this->setBoundaryConditionsV();
            this->setBoundaryConditionsVelocityGeometry();
            this->setBoundaryConditionsPGeometry();
            // print dt and residual
            std::cout << "t: " << this->t << " dt: " << this->dt << " res: " << this->res_norm;
            if (this->t - last_saved >= this->save_interval) {
                this->grid.interpolateVelocity();
                this->setBoundaryConditionsInterpolatedVelocityGeometry();
                saveVTK(this);
                last_saved += this->save_interval;
                std::cout << " VTK saved!!!";
            }
            std::cout << std::endl;
            this->duration += static_cast<int>(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count()) - this->lastTimestamp;
            this->res_norm_over_time(this->duration) = this->res_norm;
            this->t = this->t + this->dt;
            this->it_wo_pressure_solver++;
        }

        this->setBoundaryConditionsU();
        this->setBoundaryConditionsV();
        this->setBoundaryConditionsVelocityGeometry();

        this->grid.interpolateVelocity();

        this->setBoundaryConditionsInterpolatedVelocityGeometry();
        this->setBoundaryConditionsP();
        this->setBoundaryConditionsPGeometry();

        this->res_norm_over_it_with_pressure_solver.conservativeResize(this->it);
        this->res_norm_over_it_without_pressure_solver.conservativeResize(this->it_wo_pressure_solver);
        this->res_norm_over_time.conservativeResize(this->duration);

        return;
    }

    void FluidSimulation::saveData() {
        Kernel::saveMatrix("u.dat", &this->grid.u_interpolated);
        Kernel::saveMatrix("v.dat", &this->grid.v_interpolated);
        Kernel::saveMatrix("p.dat", &this->grid.p);
        Kernel::saveVector("residuals_with_pressure_solver.dat", &this->res_norm_over_it_with_pressure_solver);
        Kernel::saveVector("residuals_without_pressure_solver.dat", &this->res_norm_over_it_without_pressure_solver);
        Kernel::saveVector("residuals_over_time.dat", &this->res_norm_over_time);
    }
}