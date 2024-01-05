#include <bitset>
#include <iostream>
#include "simulations.h"

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
    int height = floor(this->grid.jmax / 8.0);
    int distanceTop = floor((this->grid.jmax - height) / 2.0);
    int distanceBottom = distanceTop + height;
    int distanceLeft = floor((this->grid.jmax - width) / 2.0);
    int distanceRight = distanceLeft + width;
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

    FluidSimulation::run();

    return;
}
