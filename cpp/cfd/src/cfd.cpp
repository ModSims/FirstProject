#include <iostream>
#include "cfd.h"

namespace CFD {
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
    void selectDtAccordingToStabilityCondition(StaggeredGrid *grid, Simulation *sim) {
        double left = (sim->Re/2) * pow(((1/pow(grid->dx(), 2))+(1/pow(grid->dy(), 2))), -1);
        double middle = grid->dx() / grid->findMaxAbsoluteU();
        double right = grid->dy() / grid->findMaxAbsoluteV();
        sim->dt = sim->tau * std::min(std::min(left, middle), right);
    }

    namespace ME_X {
        // Slide 15
        double uu_x(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (
                    (1/grid->dx())*(pow((grid->u(i,j) + grid->u(i+1,j))/2,2) - pow((grid->u(i-1,j)+grid->u(i,j))/2,2))
                )
                +
                (
                    (sim->alpha/grid->dx())*(std::abs(grid->u(i,j) + grid->u(i+1,j))*(grid->u(i,j) - grid->u(i+1,j))/4 - std::abs(grid->u(i-1,j) + grid->u(i,j))*(grid->u(i-1,j) - grid->u(i,j))/4)
                )
            );
        }

        double uv_y(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (
                    (1/grid->dy())*((grid->v(i,j) + grid->v(i+1,j))*(grid->u(i,j) + grid->u(i,j+1))/4 - (grid->v(i,j-1) + grid->v(i+1,j-1))*(grid->u(i,j-1) + grid->u(i,j))/4)
                )
                +
                (
                    (sim->alpha/grid->dy())*(std::abs(grid->v(i,j) + grid->v(i+1,j))*(grid->u(i,j) - grid->u(i,j+1))/4 - std::abs(grid->v(i,j-1) + grid->v(i+1,j-1))*(grid->u(i,j-1) - grid->u(i,j))/4)
                )
            );
        }

        double uu_xx(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (grid->u(i+1,j) - 2*grid->u(i,j) + grid->u(i-1,j))/pow(grid->dx(), 2)
            );
        }

        double uu_yy(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (grid->u(i,j+1) - 2*grid->u(i,j) + grid->u(i,j-1))/pow(grid->dy(), 2)
            );
        }

        double p_x(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (grid->p(i+1,j) - grid->p(i,j))/grid->dx()
            );
        }
    }

    namespace ME_Y {
        double uv_x(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (
                    (1/grid->dx())*((grid->u(i,j) + grid->u(i,j+1))*(grid->v(i,j) + grid->v(i+1,j))/4 - (grid->u(i-1,j) + grid->u(i-1,j+1))*(grid->v(i-1,j) + grid->v(i,j))/4)
                )
                +
                (
                    (sim->alpha/grid->dx())*(std::abs(grid->u(i,j) + grid->u(i,j+1))*(grid->v(i,j) - grid->v(i+1,j))/4 - std::abs(grid->u(i-1,j) + grid->u(i-1,j+1))*(grid->v(i-1,j) - grid->v(i,j))/4)
                )
            );
        }

        double vv_y(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (
                    (1/grid->dy())*(pow((grid->v(i,j) + grid->v(i,j+1))/2,2) - pow((grid->v(i,j-1)+grid->v(i,j))/2,2))
                )
                +
                (
                    (sim->alpha/grid->dy())*(std::abs(grid->v(i,j) + grid->v(i,j+1))*(grid->v(i,j) - grid->v(i,j+1))/4 - std::abs(grid->v(i,j-1) + grid->v(i,j))*(grid->v(i,j-1) - grid->v(i,j))/4)
                )
            );
        }

        double vv_xx(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (grid->v(i+1,j) - 2*grid->v(i,j) + grid->v(i-1,j))/pow(grid->dx(), 2)
            );
        }

        double vv_yy(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (grid->v(i,j+1) - 2*grid->v(i,j) + grid->v(i,j-1))/pow(grid->dy(), 2)
            );
        }

        double p_y(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (grid->p(i,j+1) - grid->p(i,j))/grid->dy()
            );
        }
    }

    namespace CE {
        double u_x(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (grid->u(i,j) - grid->u(i-1,j))/(grid->dx())
            );
        }

        double v_y(StaggeredGrid *grid, Simulation *sim, int i, int j) {
            return (
                (grid->v(i,j) - grid->v(i,j-1))/(grid->dy())
            );
        }
    }

    void computeF(StaggeredGrid *grid, Simulation *sim) {
        // Boundary conditions
        for (int j = 0; j < grid->jmax + 3; j++) {
            grid->F(0, j) = grid->u(0, j);
            grid->F(grid->imax, j) = grid->u(grid->imax, j);
        }

        for (int i = 1; i < grid->imax + 1; i++) {
            for (int j = 1; j < grid->jmax + 2; j++) {
                grid->F(i,j) = grid->u(i,j) + sim->dt * (
                    (1/sim->Re) * (ME_X::uu_xx(grid, sim, i, j) + ME_X::uu_yy(grid, sim, i, j)) 
                    - 
                    ME_X::uu_x(grid, sim, i, j)
                    -
                    ME_X::uv_y(grid, sim, i, j)
                );
            }
        }
    }

    void computeG(StaggeredGrid *grid, Simulation *sim) {
        // Boundary conditions
        for (int i = 0; i < grid->imax + 3; i++) {
            grid->G(i, 0) = grid->v(i, 0);
            grid->G(i, grid->jmax) = grid->v(i, grid->jmax);
        }

        for (int i = 1; i < grid->imax + 2; i++) {
            for (int j = 1; j < grid->jmax + 1; j++) {
                grid->G(i,j) = grid->v(i,j) + sim->dt * (
                    (1/sim->Re) * (ME_Y::vv_xx(grid, sim, i, j) + ME_Y::vv_yy(grid, sim, i, j)) 
                    - 
                    ME_Y::uv_x(grid, sim, i, j)
                    -
                    ME_Y::vv_y(grid, sim, i, j)
                );
            }
        }
    }

    void computeRHS(StaggeredGrid *grid, Simulation *sim) {
        for (int i = 1; i < grid->imax + 1; i++) {
            for (int j = 1; j < grid->jmax + 1; j++) {
                grid->RHS(i,j) = (1/sim->dt)*((grid->F(i,j) - grid->F(i-1,j))/grid->dx()+(grid->G(i,j) - grid->G(i,j-1))/grid->dy());
            }
        }
    }

    void updateStepLGLS(StaggeredGrid *grid, Simulation *sim) {
        // jacobi smoother
        for (int i = 1; i < grid->imax + 1; i++) {
            for (int j = 1; j < grid->jmax + 1; j++) {
                grid->p(i, j) = (
                    (1/(-2*pow(grid->dx(), 2) - 2*pow(grid->dy(), 2)))
                    *
                    (
                        grid->RHS(i,j)*pow(grid->dx(), 2)*pow(grid->dy(), 2) - pow(grid->dy(), 2)*(grid->p(i+1,j) + grid->p(i-1,j)) - pow(grid->dx(), 2)*(grid->p(i,j+1) + grid->p(i,j-1))
                    )
                );
            }
        }
    }

    void computeResidual(StaggeredGrid *grid, Simulation *sim) {
        sim->res = (grid->p - grid->po).norm();
    }

    void computeU(StaggeredGrid *grid, Simulation *sim) {
        for (int i = 1; i < grid->imax + 1; i++) {
            for (int j = 1; j < grid->jmax + 2; j++) {
                grid->u(i,j) = grid->F(i,j) - (sim->dt/grid->dx()) * (grid->p(i+1,j) - grid->p(i,j));
            }
        }
    }

    void computeV(StaggeredGrid *grid, Simulation *sim) {
        for (int i = 1; i < grid->imax + 2; i++) {
            for (int j = 1; j < grid->jmax + 1; j++) {
                grid->v(i,j) = grid->G(i,j) - (sim->dt/grid->dy()) * (grid->p(i,j+1) - grid->p(i,j));
            }
        }
    }
}