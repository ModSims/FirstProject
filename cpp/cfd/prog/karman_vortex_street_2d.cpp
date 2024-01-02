#include <iostream>
#include "simulations.h"

using namespace CFD;

int main(int argc, char* argv[]) {
    FluidParams flow_params = FluidParams("karman_vortex_street_2d", argc, argv);
    KarmanVortexStreet2D sim = KarmanVortexStreet2D(flow_params);
    sim.run();
    sim.saveData();
}
