#include <iostream>
#include "cfd.h"
#include "kernel.h"

using namespace CFD;
using namespace Kernel; 

void setBoundaryConditionsU(StaggeredGrid *grid) {
    // Everything no-slip in u-rows
    for (int j = 0; j < grid->jmax + 3; j++) {
        grid->u(0, j) = 0.0;
        grid->u(grid->imax, j) = 0.0;
    }

    // Interpolate with inner wall
    for (int i = 0; i < grid->imax + 2; i++) {
        grid->u(i, 0) = -grid->u(i, 1);
        grid->u(i, grid->jmax + 1) = -grid->u(i, grid->jmax);
    }

    // Moving wall u_i,jmax+1 = 2.0 - ui,jmax
    for (int i = 0; i < grid->jmax + 2; i++) {
        grid->u(i,grid->jmax+1) = 2.0 - grid->u(i,grid->jmax);
    }
}

void setBoundaryConditionsV(StaggeredGrid *grid) {
    // Everything no-slip in v-rows
    for (int j = 0; j < grid->jmax + 2; j++) {
        grid->v(0, j) = -grid->v(1, j);
        grid->v(grid->imax + 1, j) = -grid->v(grid->imax, j);
    }
    
    // Interpolate with inner wall
    for (int i = 0; i < grid->imax + 3; i++) {
        grid->v(i, 0) = 0.0;
        grid->v(i, grid->jmax) = 0.0;
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
}

int main(int argc, char** argv) {
    // CONSTANTS
    const int imax = 50;
    const int jmax = 50;
    const double xlength = 1.0;
    const double ylength = 1.0;
    const double t_end = 50.0;
    const double tau = 0.5;
    const double eps = 1e-3;
    const double omg = 1.7;
    const int itermax = 100;
    const double alpha = 0.5;
    const double Re = 100.0; 

    // VARIABLES
    double t = 0;
    double dt = 0.05;
    double res = 99999;
    
    Simulation sim(imax, jmax, xlength, ylength, t_end, tau, eps, omg, itermax, alpha, Re, t, dt, res);
    StaggeredGrid grid(imax, jmax, xlength, ylength);

    int n = 0;

    while(sim.t < sim.t_end) {
        n = 0;
        selectDtAccordingToStabilityCondition(&grid, &sim);
        // print dt and residual
        std::cout << "t: " << sim.t << " dt: " << sim.dt << " res: " << sim.res << std::endl;
        setBoundaryConditionsU(&grid);
        setBoundaryConditionsV(&grid);
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

    saveMatrix("u.dat", &grid.u);
    saveMatrix("v.dat", &grid.v);
    saveMatrix("p.dat", &grid.p);

    return 0;
}
