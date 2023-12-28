#include <bitset>
#include <iostream>
#include "cfd.h"

using namespace CFD;

void FlowOverStep2D::setBoundaryConditionsU() {
    // Inflow at left boundary (Only half)
    for (int j = this->grid.jmax/2; j < this->grid.jmax + 3; j++) {
        this->grid.u(0, j) = 1.0;
    }

    for (int j = 0; j < this->grid.jmax/2; j++) {
        this->grid.u(0, j) = 0.0;
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

void FlowOverStep2D::setBoundaryConditionsV() {
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

void FlowOverStep2D::setBoundaryConditionsP() {
    for (int i = 0; i < this->grid.imax + 2; i++) {
        this->grid.p(i, 0) = this->grid.p(i, 1);
        this->grid.p(i, this->grid.jmax + 1) = this->grid.p(i, this->grid.jmax);
    }
    for (int j = 0; j < this->grid.jmax + 2; j++) {
        this->grid.p(0, j) = this->grid.p(1, j);
        this->grid.p(this->grid.imax + 1, j) = this->grid.p(this->grid.imax, j);
    }
}

void FlowOverStep2D::run() {
    // Set initially u to 1.0
    for (int i = 0; i < this->grid.imax + 2; i++) {
        for (int j = 0; j < this->grid.jmax + 3; j++) {
            this->grid.u(i, j) = 1.0;
        }
    }

    // Manage Flag Field with Bitmasks
    for (int i = 1; i < this->grid.imax + 1; i++) {
        for (int j = 1; j < this->grid.jmax + 1; j++) {
            // Fluid cells have fluid neighbors
            this->grid.flag_field(i, j) = FlagFieldMask::CELL_FLUID;

            if (i < (this->grid.imax+2) / 8 && j < (this->grid.jmax+2) / 2) {
                // Reset flag field
                this->grid.flag_field(i, j) = FlagFieldMask::CELL_OBSTACLE;

                // Obstacle cells have fluid neighbors
                if (i == ((this->grid.imax+2)/8 - 1)) {
                    // East neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_EAST;
                }

                if (j == ((this->grid.jmax+2)/2 - 1)) {
                    // North neighbor
                    this->grid.flag_field(i, j) |= FlagFieldMask::FLUID_NORTH;
                }
            }
            else if (i == (this->grid.imax+2) / 8 && j < (this->grid.jmax+2) / 2) {
                // Fluid cells have obstacle neighbors at West
                this->grid.flag_field(i, j) &= FlagFieldMask::OBSTACLE_WEST;


            }
            else if (i < (this->grid.imax+2) / 8 && j == (this->grid.jmax+2) / 2) {
                // Fluid cells have obstacle neighbors at South
                this->grid.flag_field(i, j) &= FlagFieldMask::OBSTACLE_SOUTH;
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
