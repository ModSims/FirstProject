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
                    tmp_p += this->grid.p(i+1, j);
                    counter++;
                }
                if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_EAST) {
                    tmp_p += this->grid.p(i-1, j);
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
    int tmp_u = 0;
    int counter_u = 0;
    int tmp_v = 0;
    int counter_v = 0;
    for (int i = 1; i < this->grid.imax + 1; i++) {
        for (int j = 1; j < this->grid.jmax + 1; j++) {
            tmp_u = 0;
            counter_u = 0;
            tmp_v = 0;
            counter_v = 0;
            // check if is obstacle
            if ((this->grid.flag_field(i, j) & FlagFieldMask::MASK_CELL_TYPE) == (FlagFieldMask::CELL_OBSTACLE & FlagFieldMask::MASK_CELL_TYPE)) {
                // obstacle cell
                if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_NORTH) {
                    tmp_u += this->grid.u_interpolated(i, j+1);
                    counter_u++;
                    tmp_v += this->grid.v_interpolated(i, j+1);
                    counter_v++;
                }
                if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_SOUTH) {
                    tmp_u += this->grid.u_interpolated(i, j-1);
                    counter_u++;
                    tmp_v += this->grid.v_interpolated(i, j-1);
                    counter_v++;
                }
                if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_WEST) {
                    tmp_u += this->grid.u_interpolated(i+1, j);
                    counter_u++;
                    tmp_v += this->grid.v_interpolated(i+1, j);
                    counter_v++;
                }
                if (this->grid.flag_field(i, j) & FlagFieldMask::FLUID_EAST) {
                    tmp_u += this->grid.u_interpolated(i-1, j);
                    counter_u++;
                    tmp_v += this->grid.v_interpolated(i-1, j);
                    counter_v++;
                }
                if (counter_u > 0) {
                    this->grid.u_interpolated(i, j) = tmp_u / counter_u;
                }
                if (counter_v > 0) {
                    this->grid.v_interpolated(i, j) = tmp_v / counter_v;
                }
                if (this->grid.flag_field(i, j) == FlagFieldMask::CELL_OBSTACLE) {
                    // interior obstacle cell, so no-slip
                    this->grid.u_interpolated(i,j) = 0.0;
                    this->grid.v_interpolated(i,j) = 0.0;
                }
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
    for (int i = 1; i < this->grid.imax + 1; i++) {
        for (int j = 1; j < this->grid.jmax + 1; j++) {
            double centerX = (this->grid.jmax + 2) / 2.0;
            double centerY = (this->grid.jmax + 2) / 2.0;

            // Calculate the distance to the center of the circle
            double distanceToCenter = std::sqrt((i - centerX) * (i - centerX) +
                                                (j - centerY) * (j - centerY));

            if (distanceToCenter <= std::min((this->grid.imax+2) / 6.0, (this->grid.jmax+2) / 6.0)) {
                // Inside the circular obstacle
                this->grid.flag_field(i, j) = FlagFieldMask::CELL_OBSTACLE;

                // Determine quadrant
                bool inFirstQuadrant = (i > (this->grid.imax + 2) / 2.0 && j > (this->grid.jmax + 2) / 2.0);
                bool inSecondQuadrant = (i < (this->grid.imax + 2) / 2.0 && j > (this->grid.jmax + 2) / 2.0);
                bool inThirdQuadrant = (i < (this->grid.imax + 2) / 2.0 && j < (this->grid.jmax + 2) / 2.0);
                bool inFourthQuadrant = (i > (this->grid.imax + 2) / 2.0 && j < (this->grid.jmax + 2) / 2.0);

                // Set fluid neighbors based on quadrant

                if (inFirstQuadrant) {
                    // West neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_EAST;
                    // South neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_SOUTH;
                }

                if (inSecondQuadrant) {
                    // East neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_EAST;
                    // South neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_NORTH;
                }

                if (inThirdQuadrant) {
                    // East neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_WEST;
                    // North neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_NORTH;
                }

                if (inFourthQuadrant) {
                    // West neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_WEST;
                    // North neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_SOUTH;
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
