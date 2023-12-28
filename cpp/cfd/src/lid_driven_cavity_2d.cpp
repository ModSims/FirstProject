#include <iostream>
#include "cfd.h"

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

void LidDrivenCavity2D::run() {
    int n = 0;
    while(this->t < this->t_end) {
        n = 0;
        this->selectDtAccordingToStabilityCondition();
        // print dt and residual
        std::cout << "t: " << this->t << " dt: " << this->dt << " res: " << this->res_norm << std::endl;
        this->setBoundaryConditionsU();
        this->setBoundaryConditionsV();
        this->computeF();
        this->computeG();
        this->computeRHS();
        while ((this->res_norm > this->eps || this->res_norm == 0) && n < this->itermax) {
            this->setBoundaryConditionsP();
            this->solveWithJacobi();
            this->computeResidual();
            n++;
        }
        this->computeU();
        this->computeV();
        this->grid.po = this->grid.p;
        this->t = this->t + this->dt;
        if (std::abs(t - std::round(t)) < 0.1) {
            this->grid.interpolateVelocity();
            saveVTK(this);
        }

    }

    this->setBoundaryConditionsU();
    this->setBoundaryConditionsV();

    this->grid.interpolateVelocity();

    return;
}
