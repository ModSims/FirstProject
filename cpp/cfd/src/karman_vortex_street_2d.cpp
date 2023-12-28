#include <bitset>
#include <iostream>
#include "cfd.h"

using namespace CFD;

void KarmanVortexStreet2D::setBoundaryConditionsU() {
    // Inflow and outflow at left and right boundary
    for (int j = 0; j < this->grid.jmax + 3; j++) {
        // Inflow at left boundary (Left wall)
        this->grid.u(0, j) = 1.0;
        // Outflow at right boundary (Right wall)
        this->grid.u(this->grid.imax + 1, j) = this->grid.u(this->grid.imax, j);
    }

    // no-slip at top and bottom
    for (int i = 0; i < this->grid.imax + 2; i++) {
        // Bottom wall
        this->grid.u(i, 0) = -this->grid.u(i, 1);
        // Top wall
        this->grid.u(i, this->grid.jmax + 2) = -this->grid.u(i, this->grid.jmax + 1);
    }
}

void KarmanVortexStreet2D::setBoundaryConditionsV() {
    // Inflow and outflow at left and right boundary
    for (int j = 0; j < this->grid.jmax + 2; j++) {
        // Inflow at left boundary (Left wall)
        this->grid.v(0, j) = -this->grid.v(1, j);
        // Outflow at right boundary (Right wall)
        this->grid.v(this->grid.imax + 1, j) = this->grid.v(this->grid.imax, j);
    }

    // no-slip at top and bottom
    for (int i = 0; i < this->grid.imax + 3; i++) {
        // Bottom wall
        this->grid.v(i, 0) = 0.0;
        // Top wall
        this->grid.v(i, this->grid.jmax + 1) = 0.0;
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
}

void KarmanVortexStreet2D::run() {
    // Manage Flag Field with Bitmasks
    // Square in the middle with fifth of the size of the domain
    int width = floor(this->grid.jmax / 4.0);
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
        this->selectDtAccordingToStabilityCondition();
        // print dt and residual
        std::cout << "t: " << this->t << " dt: " << this->dt << " res: " << this->res_norm << std::endl;
        this->setBoundaryConditionsU();
        this->setBoundaryConditionsV();
        this->setBoundaryConditionsVelocityGeometry();
        this->computeF();
        this->computeG();
        this->setBoundaryConditionsVelocityGeometry();
        this->computeRHS();
        while ((this->res_norm > this->eps || this->res_norm == 0) && n < this->itermax) {
            this->setBoundaryConditionsP();
            this->setBoundaryConditionsPGeometry();
            this->solveWithJacobi();
            this->computeResidual();
            n++;
        }
        this->computeU();
        this->computeV();
        this->grid.po = this->grid.p;
        this->t = this->t + this->dt;
        this->setBoundaryConditionsU();
        this->setBoundaryConditionsV();
        this->setBoundaryConditionsVelocityGeometry();
        this->setBoundaryConditionsP();
        this->setBoundaryConditionsPGeometry();
        if (std::abs(t - std::round(t)) < 0.1) {
            this->grid.interpolateVelocity();
            this->setBoundaryConditionsInterpolatedVelocityGeometry();
            saveVTK(this);
        }
    }

    this->setBoundaryConditionsU();
    this->setBoundaryConditionsV();
    this->setBoundaryConditionsVelocityGeometry();

    this->grid.interpolateVelocity();

    this->setBoundaryConditionsInterpolatedVelocityGeometry();
    this->setBoundaryConditionsP();
    this->setBoundaryConditionsPGeometry();

    return;
}
