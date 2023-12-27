#pragma once
#include "cfd.h"
#include "kernel.h"

using namespace CFD;
using namespace Kernel;

class LidDrivenCavity2D {
    public:
    void setBoundaryConditionsU();
    void setBoundaryConditionsV();
    void setBoundaryConditionsP();
    void run(const int imax, const int jmax, const double Re, const double xlength, const double ylength,
             const double t_end, const double tau, const double eps, const double omg, const int itermax,
             const double alpha);
    Simulation sim;
    StaggeredGrid grid;
};