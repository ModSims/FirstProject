#include <bitset>
#include <iostream>
#include "simulations.h"

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

    FluidSimulation::run();

    return;
}
