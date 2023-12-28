#include <iostream>
#include "cfd.h"
#include "kernel.h"
#include "argparse.h"

using namespace CFD;

int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("plane_shear_flow_2d");
    program.add_argument("-i", "--imax").help("imax").default_value(100).action([](const std::string& value) { return std::stoi(value); });
    program.add_argument("-j", "--jmax").help("jmax").default_value(20).action([](const std::string& value) { return std::stoi(value); });
    program.add_argument("-x", "--xlength").help("xlength").default_value(10.0).action([](const std::string& value) { return std::stod(value); });
    program.add_argument("-y", "--ylength").help("ylength").default_value(2.0).action([](const std::string& value) { return std::stod(value); });
    program.add_argument("-z", "--t_end").help("t_end").default_value(10.0).action([](const std::string& value) { return std::stod(value); });
    program.add_argument("-u", "--tau").help("tau").default_value(0.5).action([](const std::string& value) { return std::stod(value); });
    program.add_argument("-e", "--eps").help("eps").default_value(1e-3).action([](const std::string& value) { return std::stod(value); });
    program.add_argument("-o", "--omg").help("omg").default_value(1.7).action([](const std::string& value) { return std::stod(value); });
    program.add_argument("-m", "--itermax").help("itermax").default_value(100).action([](const std::string& value) { return std::stoi(value); });
    program.add_argument("-a", "--alpha").help("alpha").default_value(0.9).action([](const std::string& value) { return std::stod(value); });
    program.add_argument("-r", "--Re").help("Re").default_value(100.0).action([](const std::string& value) { return std::stod(value); });
    program.add_argument("-t", "--t").help("t").default_value(0.0).action([](const std::string& value) { return std::stod(value); });
    program.add_argument("-d", "--dt").help("dt").default_value(0.05).action([](const std::string& value) { return std::stod(value); });

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    PlaneShearFlow2D sim = PlaneShearFlow2D(
        program.get<int>("--imax"),
        program.get<int>("--jmax"),
        program.get<double>("--xlength"),
        program.get<double>("--ylength"),
        program.get<double>("--t_end"),
        program.get<double>("--tau"),
        program.get<double>("--eps"),
        program.get<double>("--omg"),
        program.get<int>("--itermax"),
        program.get<double>("--alpha"),
        program.get<double>("--Re"),
        program.get<double>("--t"),
        program.get<double>("--dt")
    );
    sim.run();
    sim.saveMatrices();
}
