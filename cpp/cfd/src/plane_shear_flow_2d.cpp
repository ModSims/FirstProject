#include <iostream>
#include "simulations.h"

using namespace CFD;

void PlaneShearFlow2D::setBoundaryConditionsU() {
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

void PlaneShearFlow2D::setBoundaryConditionsV() {
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

void PlaneShearFlow2D::setBoundaryConditionsP() {
    for (int i = 0; i < this->grid.imax + 2; i++) {
        this->grid.p(i, 0) = this->grid.p(i, 1);
        this->grid.p(i, this->grid.jmax + 1) = this->grid.p(i, this->grid.jmax);
    }
    for (int j = 0; j < this->grid.jmax + 2; j++) {
        this->grid.p(0, j) = this->grid.p(1, j);
        this->grid.p(this->grid.imax + 1, j) = this->grid.p(this->grid.imax, j);
    }
}
