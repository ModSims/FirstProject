#include <bitset>
#include <iostream>
#include "cfd.h"
#include "kernel.h"

using namespace CFD;
using namespace Kernel; 

void setBoundaryConditionsU(StaggeredGrid *grid) {
    // Inflow at left boundary
    for (int j = 0; j < grid->jmax + 3; j++) {
        grid->u(0, j) = 1.0;
    }

    // Outflow at right boundary
    for (int j = 0; j < grid->jmax + 3; j++) {
        grid->u(grid->imax, j) = grid->u(grid->imax-1, j);
    }

    // no-slip at top and bottom
    for (int i = 0; i < grid->imax + 2; i++) {
        grid->u(i, 0) = -grid->u(i, 1);
        grid->u(i, grid->jmax + 1) = -grid->u(i, grid->jmax);
    }
}

void setBoundaryConditionsV(StaggeredGrid *grid) {
    // Inflow at left boundary
    for (int j = 0; j < grid->jmax + 2; j++) {
        grid->v(0, j) = 0.0;
    }

    // Outflow at right boundary
    for (int j = 0; j < grid->jmax + 2; j++) {
        grid->v(grid->imax + 1, j) = grid->v(grid->imax, j);
    }

    // no-slip at top and bottom
    for (int i = 0; i < grid->imax + 3; i++) {
        grid->v(i, 0) = 0.0;
        grid->v(i, grid->jmax + 1) = 0.0;
    }
}

void setBoundaryConditionsVelocityGeometry(StaggeredGrid *grid) {
    // Geometry boundaries
    for (int i = 1; i < grid->imax + 1; i++) {
        for (int j = 1; j < grid->jmax + 1; j++) {
            // check if is obstacle
            if ((grid->flag_field(i, j) & FlagFieldMask::MASK_CELL_TYPE) == (FlagFieldMask::CELL_OBSTACLE & FlagFieldMask::MASK_CELL_TYPE)) {
                // obstacle cell
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_NORTH) {
                    grid->u(i, j) = -grid->u(i, j+1);
                    grid->u(i-1, j) = -grid->u(i-1, j+1);
                    grid->v(i, j) = 0.0;
                }
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_SOUTH) {
                    grid->u(i, j) = -grid->u(i, j-1);
                    grid->u(i-1, j) = -grid->u(i-1, j-1);
                    grid->v(i, j) = 0.0;
                }
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_WEST) {
                    grid->u(i-1, j) = 0.0;
                    grid->v(i, j) = -grid->v(i-1, j);
                    grid->v(i, j-1) = -grid->v(i-1, j-1);
                }
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_EAST) {
                    grid->u(i, j) = 0.0;
                    grid->v(i, j) = -grid->v(i+1, j);
                    grid->v(i, j-1) = -grid->v(i+1, j-1);
                }
                if (grid->flag_field(i, j) == FlagFieldMask::CELL_OBSTACLE) {
                    // interior obstacle cell, so no-slip
                    grid->u(i,j) = 0.0;
                    grid->v(i,j) = 0.0;
                }
            }
        }
    }
}

void setBoundaryConditionsP(StaggeredGrid *grid) {
    for (int i = 0; i < grid->imax + 2; i++) {
        grid->p(i, 0) = grid->p(i, 1);
        grid->p(i, grid->jmax + 1) = grid->p(i, grid->jmax);
    }
    for (int j = 0; j < grid->jmax + 2; j++) {
        grid->p(0, j) = grid->p(1, j);
        grid->p(grid->imax + 1, j) = grid->p(grid->imax, j);
    }

    // Geometry boundaries
    int tmp_p = 0;
    int counter = 0;
    for (int i = 1; i < grid->imax + 1; i++) {
        for (int j = 1; j < grid->jmax + 1; j++) {
            tmp_p = 0;
            counter = 0;
            // check if is obstacle
            if ((grid->flag_field(i, j) & FlagFieldMask::MASK_CELL_TYPE) == (FlagFieldMask::CELL_OBSTACLE & FlagFieldMask::MASK_CELL_TYPE)) {
                // obstacle cell
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_NORTH) {
                    tmp_p += grid->p(i, j+1);
                    counter++;
                }
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_SOUTH) {
                    tmp_p += grid->p(i, j-1);
                    counter++;
                }
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_WEST) {
                    tmp_p += grid->p(i+1, j);
                    counter++;
                }
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_EAST) {
                    tmp_p += grid->p(i-1, j);
                    counter++;
                }
                if (counter > 0) {
                    grid->p(i, j) = tmp_p / counter;
                }
                if (grid->flag_field(i, j) == FlagFieldMask::CELL_OBSTACLE) {
                    // interior obstacle cell, so no-slip
                    grid->p(i,j) = 0.0;
                }
            }
        }
    }
}

void setBoundaryConditionsInterpolatedVelocityGeometry(StaggeredGrid *grid) {
    // Geometry boundaries
    int tmp_u = 0;
    int counter_u = 0;
    int tmp_v = 0;
    int counter_v = 0;
    for (int i = 1; i < grid->imax + 1; i++) {
        for (int j = 1; j < grid->jmax + 1; j++) {
            tmp_u = 0;
            counter_u = 0;
            tmp_v = 0;
            counter_v = 0;
            // check if is obstacle
            if ((grid->flag_field(i, j) & FlagFieldMask::MASK_CELL_TYPE) == (FlagFieldMask::CELL_OBSTACLE & FlagFieldMask::MASK_CELL_TYPE)) {
                // obstacle cell
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_NORTH) {
                    tmp_u += grid->u_interpolated(i, j+1);
                    counter_u++;
                    tmp_v += grid->v_interpolated(i, j+1);
                    counter_v++;
                }
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_SOUTH) {
                    tmp_u += grid->u_interpolated(i, j-1);
                    counter_u++;
                    tmp_v += grid->v_interpolated(i, j-1);
                    counter_v++;
                }
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_WEST) {
                    tmp_u += grid->u_interpolated(i+1, j);
                    counter_u++;
                    tmp_v += grid->v_interpolated(i+1, j);
                    counter_v++;
                }
                if (grid->flag_field(i, j) & FlagFieldMask::FLUID_EAST) {
                    tmp_u += grid->u_interpolated(i-1, j);
                    counter_u++;
                    tmp_v += grid->v_interpolated(i-1, j);
                    counter_v++;
                }
                if (counter_u > 0) {
                    grid->u_interpolated(i, j) = tmp_u / counter_u;
                }
                if (counter_v > 0) {
                    grid->v_interpolated(i, j) = tmp_v / counter_v;
                }
                if (grid->flag_field(i, j) == FlagFieldMask::CELL_OBSTACLE) {
                    // interior obstacle cell, so no-slip
                    grid->u_interpolated(i,j) = 0.0;
                    grid->v_interpolated(i,j) = 0.0;
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    // CONSTANTS
    const int imax = 100;
    const int jmax = 20;
    const double xlength = 10.0;
    const double ylength = 2.0;
    const double t_end = 20;
    const double tau = 0.5;
    const double eps = 1e-3;
    const double omg = 1.7;
    const int itermax = 2000;
    const double alpha = 0.9;
    const double Re = 10000.0;

    // VARIABLES
    double t = 0;
    double dt = 0.05;
    double res = 99999;
    
    Simulation sim(imax, jmax, xlength, ylength, t_end, tau, eps, omg, itermax, alpha, Re, t, dt, res);
    StaggeredGrid grid(imax, jmax, xlength, ylength);

    // Set initially u to 1.0
    for (int i = 0; i < grid.imax + 2; i++) {
        for (int j = 0; j < grid.jmax + 3; j++) {
            grid.u(i, j) = 1.0;
        }
    }

    // Manage Flag Field with Bitmasks
    for (int i = 1; i < grid.imax + 1; i++) {
        for (int j = 1; j < grid.jmax + 1; j++) {
            double centerX = (grid.jmax + 2) / 2.0;
            double centerY = (grid.jmax + 2) / 2.0;

            // Calculate the distance to the center of the circle
            double distanceToCenter = std::sqrt((i - centerX) * (i - centerX) +
                                                (j - centerY) * (j - centerY));

            if (distanceToCenter <= std::min((grid.imax+2) / 6.0, (grid.jmax+2) / 6.0)) {
                // Inside the circular obstacle
                grid.flag_field(i, j) = FlagFieldMask::CELL_OBSTACLE;

                // Determine quadrant
                bool inFirstQuadrant = (i > (grid.imax + 2) / 2.0 && j > (grid.jmax + 2) / 2.0);
                bool inSecondQuadrant = (i < (grid.imax + 2) / 2.0 && j > (grid.jmax + 2) / 2.0);
                bool inThirdQuadrant = (i < (grid.imax + 2) / 2.0 && j < (grid.jmax + 2) / 2.0);
                bool inFourthQuadrant = (i > (grid.imax + 2) / 2.0 && j < (grid.jmax + 2) / 2.0);

                // Set fluid neighbors based on quadrant

                if (inFirstQuadrant) {
                    // West neighbor
                    grid.flag_field(i, j) |= FlagFieldMask::FLUID_EAST;
                    // South neighbor
                    grid.flag_field(i, j) |= FlagFieldMask::FLUID_SOUTH;
                }

                if (inSecondQuadrant) {
                    // East neighbor
                    grid.flag_field(i, j) |= FlagFieldMask::FLUID_EAST;
                    // South neighbor
                    grid.flag_field(i, j) |= FlagFieldMask::FLUID_NORTH;
                }

                if (inThirdQuadrant) {
                    // East neighbor
                    grid.flag_field(i, j) |= FlagFieldMask::FLUID_WEST;
                    // North neighbor
                    grid.flag_field(i, j) |= FlagFieldMask::FLUID_NORTH;
                }

                if (inFourthQuadrant) {
                    // West neighbor
                    grid.flag_field(i, j) |= FlagFieldMask::FLUID_WEST;
                    // North neighbor
                    grid.flag_field(i, j) |= FlagFieldMask::FLUID_SOUTH;
                }
            } else {
                // Fluid cells have fluid neighbors
                grid.flag_field(i, j) = FlagFieldMask::CELL_FLUID;
            }
        }
    }

    // loop through flag field and print with std::bitset
    /*for (int i = 0; i < grid.imax + 2; i++) {
        for (int j = 0; j < grid.jmax + 2; j++) {
            std::cout << std::bitset<5>(grid.flag_field(i, j)) << " ";
        }
        std::cout << std::endl;
    }
    std::exit(0);*/

    int n = 0;

    while(sim.t < sim.t_end) {
        n = 0;
        selectDtAccordingToStabilityCondition(&grid, &sim);
        // print dt and residual
        std::cout << "t: " << sim.t << " dt: " << sim.dt << " res: " << sim.res << std::endl;
        setBoundaryConditionsU(&grid);
        setBoundaryConditionsV(&grid);
        setBoundaryConditionsVelocityGeometry(&grid);
        computeF(&grid, &sim);
        computeG(&grid, &sim);
        computeRHS(&grid, &sim);
        while ((sim.res > eps || sim.res == 0) && n < itermax) {
            setBoundaryConditionsP(&grid);
            updateStepLGLS(&grid, &sim);
            computeResidual(&grid, &sim);
            n++;
        }
        computeU(&grid, &sim);
        computeV(&grid, &sim);
        grid.po = grid.p;
        sim.t = sim.t + sim.dt;
    }

    setBoundaryConditionsU(&grid);
    setBoundaryConditionsV(&grid);
    setBoundaryConditionsVelocityGeometry(&grid);

    grid.interpolateVelocity();

    setBoundaryConditionsInterpolatedVelocityGeometry(&grid);
    setBoundaryConditionsP(&grid);

    saveMatrix("u.dat", &grid.u_interpolated);
    saveMatrix("v.dat", &grid.v_interpolated);
    saveMatrix("p.dat", &grid.p);

    return 0;
}
