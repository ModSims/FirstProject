#include <bitset>
#include <iostream>
#include "cfd.h"
#include "kernel.h"

using namespace CFD;
using namespace Kernel; 

void KarmanVortexStreet2D::setBoundaryConditionsU() {
    // Inflow at left boundary
    for (int j = 0; j < this->grid.jmax + 3; j++) {
        this->grid.u(0, j) = 1.0;
    }

    // Outflow at right boundary
    for (int j = 0; j < this->grid.jmax + 3; j++) {
        this->grid.u(this->grid.imax, j) = this->grid.u(this->grid.imax-1, j);
    }

    // no-slip at top and bottom
    for (int i = 0; i < this->grid.imax + 2; i++) {
        this->grid.u(i, 0) = -this->grid.u(i, 1);
        this->grid.u(i, this->grid.jmax + 1) = -this->grid.u(i, this->grid.jmax);
    }
}

void KarmanVortexStreet2D::setBoundaryConditionsV() {
    // Inflow at left boundary
    for (int j = 0; j < this->grid.jmax + 2; j++) {
        this->grid.v(0, j) = 0.0;
    }

    // Outflow at right boundary
    for (int j = 0; j < this->grid.jmax + 2; j++) {
        this->grid.v(this->grid.imax + 1, j) = this->grid.v(this->grid.imax, j);
    }

    // no-slip at top and bottom
    for (int i = 0; i < this->grid.imax + 3; i++) {
        this->grid.v(i, 0) = 0.0;
        this->grid.v(i, this->grid.jmax + 1) = 0.0;
    }
}

void KarmanVortexStreet2D::setBoundaryConditionsVelocityGeometry() {
    // Geometry boundaries
    for (int i = 1; i < this->grid.imax + 1; i++) {
        for (int j = 1; j < this->grid.jmax + 1; j++) {
            // check if is obstacle
            if ((this->grid.flag_field(i, j) & FlagFieldMask::MASK_CELL_TYPE) == (FlagFieldMask::CELL_OBSTACLE & FlagFieldMask::MASK_CELL_TYPE)) {
                // obstacle cell
                if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_NORTH) {
                    this->grid.u(i, j) = -this->grid.u(i, j+1);
                    this->grid.u(i-1, j) = -this->grid.u(i-1, j+1);
                    this->grid.v(i, j) = 0.0;
                }
                if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_SOUTH) {
                    this->grid.u(i, j) = -this->grid.u(i, j-1);
                    this->grid.u(i-1, j) = -this->grid.u(i-1, j-1);
                    this->grid.v(i, j) = 0.0;
                }
                if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_WEST) {
                    this->grid.u(i-1, j) = 0.0;
                    this->grid.v(i, j) = -this->grid.v(i-1, j);
                    this->grid.v(i, j-1) = -this->grid.v(i-1, j-1);
                }
                if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_EAST) {
                    this->grid.u(i, j) = 0.0;
                    this->grid.v(i, j) = -this->grid.v(i+1, j);
                    this->grid.v(i, j-1) = -this->grid.v(i+1, j-1);
                }
                if (this->grid.flag_field(i, j) == FlagFieldMask::CELL_OBSTACLE) {
                    // interior obstacle cell, so no-slip
                    this->grid.u(i,j) = 0.0;
                    this->grid.v(i,j) = 0.0;
                }
            }
        }
    }
}

void KarmanVortexStreet2D::setBoundaryConditionsP() {
    for (int i = 0; i < this->grid.imax + 2; i++) {
        this->grid.p(i, 0) = this->grid.p(i, 1);
        this->grid.p(i, this->grid.jmax + 1) = this->grid.p(i, this->grid.jmax);
    }
    for (int j = 0; j < this->grid.jmax + 2; j++) {
        this->grid.p(0, j) = this->grid.p(1, j);
        this->grid.p(this->grid.imax + 1, j) = this->grid.p(this->grid.imax, j);
    }

    // Geometry boundaries
    int tmp_p = 0;
    int counter = 0;
    for (int i = 1; i < this->grid.imax + 1; i++) {
        for (int j = 1; j < this->grid.jmax + 1; j++) {
            tmp_p = 0;
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
                if (this->grid.flag_field(i, j) == FlagFieldMask::CELL_OBSTACLE) {
                    // interior obstacle cell, so no-slip
                    this->grid.p(i,j) = 0.0;
                }
            }
        }
    }
}

void KarmanVortexStreet2D::setBoundaryConditionsInterpolatedVelocityGeometry() {
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

void KarmanVortexStreet2D::run() {
    // Set initially u to 1.0
    for (int i = 0; i < this->grid.imax + 2; i++) {
        for (int j = 0; j < this->grid.jmax + 3; j++) {
            this->grid.u(i, j) = 1.0;
        }
    }

    // Manage Flag Field with Bitmasks
    // Square in the middle with fifth of the size of the domain
    int width = floor(this->grid.jmax / 5.0);
    int distanceTop = floor((this->grid.jmax - width) / 2.0);
    int distanceBottom = distanceTop + width;
    int distanceLeft = distanceTop;
    int distanceRight = distanceBottom;
    for (int i = 1; i < this->grid.imax + 1; i++) {
        for (int j = 1; j < this->grid.jmax + 1; j++) {

            // Check if cell is inside the square
            if (i >= distanceLeft && i <= distanceRight && j >= distanceTop && j <= distanceBottom) {
                // Inside the square
                this->grid.flag_field(i, j) = FlagFieldMask::CELL_OBSTACLE;

                // Top Layer of the square
                if (j == distanceTop) {
                    // South neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_SOUTH;
                }
                if (j == distanceBottom) {
                    // North neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_NORTH;
                }
                if (i == distanceLeft) {
                    // West neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_WEST;
                }
                if (i == distanceRight) {
                    // East neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_EAST;
                } 
            } else {
                // Fluid cells have fluid neighbors
                this->grid.flag_field(i, j) = FlagFieldMask::CELL_FLUID;
            }
        }
    }

    /*// loop through flag field and print with std::bitset
    for (int i = 0; i < this->grid.imax + 2; i++) {
        for (int j = 0; j < this->grid.jmax + 2; j++) {
            std::cout << std::bitset<5>(this->grid.flag_field(i, j)) << " ";
        }
        std::cout << std::endl;
    }
    std::exit(0);*/

    int n = 0;

    while(this->t < this->t_end) {
        n = 0;
        selectDtAccordingToStabilityCondition();
        // print dt and residual
        std::cout << "t: " << this->t << " dt: " << this->dt << " res: " << this->res << std::endl;
        setBoundaryConditionsU();
        setBoundaryConditionsV();
        setBoundaryConditionsVelocityGeometry();
        computeF();
        computeG();
        computeRHS();
        while ((this->res > this->eps || this->res == 0) && n < this->itermax) {
            setBoundaryConditionsP();
            updateStepLGLS();
            computeResidual();
            n++;
        }
        computeU();
        computeV();
        this->grid.po = this->grid.p;
        this->t = this->t + this->dt;
        this->setBoundaryConditionsU();
        this->setBoundaryConditionsV();
        if (std::abs(t - std::round(t)) < 0.1) {
            this->grid.interpolateVelocity();
            setBoundaryConditionsInterpolatedVelocityGeometry();
            setBoundaryConditionsP();
            saveVTK(this);
        }
    }

    setBoundaryConditionsU();
    setBoundaryConditionsV();
    setBoundaryConditionsVelocityGeometry();

    this->grid.interpolateVelocity();

    setBoundaryConditionsInterpolatedVelocityGeometry();
    setBoundaryConditionsP();
    

    return;
}
