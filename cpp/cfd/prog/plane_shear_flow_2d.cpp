#include <iostream>
#include "simulations.h"

using namespace CFD;

int main(int argc, char* argv[]) {
    FluidParams flow_params = FluidParams("plane_shear_flow_2d", argc, argv);
    PlaneShearFlow2D sim = PlaneShearFlow2D(flow_params);
    sim.run();
    sim.saveData();
}
