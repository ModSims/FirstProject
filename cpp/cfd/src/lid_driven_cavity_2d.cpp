#include <iostream>
#include "cfd.h"


#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace CFD;

void LidDrivenCavity2D::setBoundaryConditionsU() {
    // Everything no-slip in u-rows
    for (int j = 0; j < this->grid.jmax + 3; j++) {
        // Left wall
        this->grid.u(0, j) = 0.0;
        // Right wall
        this->grid.u(this->grid.imax, j) = 0.0;
    }

    // Interpolate with inner wall
    for (int i = 0; i < this->grid.imax + 2; i++) {
        // Bottom wall
        this->grid.u(i, 0) = -this->grid.u(i, 1);
    }

    // Moving wall u_i,jmax+1 = 2.0 - ui,jmax
    for (int i = 0; i < this->grid.imax + 2; i++) {
        // Top wall
        this->grid.u(i,this->grid.jmax+1) = 2.0 - this->grid.u(i,this->grid.jmax);
    }
}

void LidDrivenCavity2D::setBoundaryConditionsV() {
    // Everything no-slip in v-rows
    for (int j = 0; j < this->grid.jmax + 2; j++) {
        // Left wall
        this->grid.v(0, j) = -this->grid.v(1, j);
        // Right wall
        this->grid.v(this->grid.imax + 1, j) = -this->grid.v(this->grid.imax, j);
    }
    
    // Interpolate with inner wall
    for (int i = 0; i < this->grid.imax + 3; i++) {
        // Bottom wall
        this->grid.v(i, 0) = 0.0;
        // Top wall
        this->grid.v(i, this->grid.jmax) = 0.0;
    }
}

void LidDrivenCavity2D::setBoundaryConditionsP() {
    for (int i = 0; i < this->grid.imax + 2; i++) {
        this->grid.p(i, 0) = this->grid.p(i, 1);
        this->grid.p(i, this->grid.jmax + 1) = this->grid.p(i, this->grid.jmax);
    }
    for (int j = 0; j < this->grid.jmax + 2; j++) {
        this->grid.p(0, j) = this->grid.p(1, j);
        this->grid.p(this->grid.imax + 1, j) = this->grid.p(this->grid.imax, j);
    }
}

