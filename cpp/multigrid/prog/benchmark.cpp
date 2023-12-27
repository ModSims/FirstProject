#include <iostream>
#include "benchmark.h"

int main() {
    std::cout << "Multigrid Benchmark with Lid Driven Cavity Simulation" << std::endl;

    LidDrivenCavity2D obj;

    std::cout << "Running simulation..." << std::endl;
    obj.run(
        // CONSTANTS
        50, // imax
        50, // jmax
        100.0, // Re
        1.0, // xlength
        1.0, // ylength
        5.0, // t_end
        0.5, // tau
        1e-3, // eps
        1.7, // omg
        100, // itermax
        0.5 // alpha
    );
    std::cout << "Simulation residual: " << obj.sim.res << std::endl;
    std::cout << "Simulation finished." << std::endl;

    return 0;
}
