#include <iostream>
#include "simulations.h"

using namespace CFD;

int main(int argc, char* argv[]) {
    FluidParams flow_params = FluidParams("lid_driven_cavity_2d", argc, argv);
    LidDrivenCavity2D sim = LidDrivenCavity2D(flow_params);
    sim.run();
    sim.saveData();
}
