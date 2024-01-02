#include <iostream>
#include "simulations.h"

using namespace CFD;

int main(int argc, char* argv[]) {
    FluidParams flow_params = FluidParams("flow_over_step_2d", argc, argv);
    FlowOverStep2D sim = FlowOverStep2D(flow_params);
    sim.run();
    sim.saveData();
}
