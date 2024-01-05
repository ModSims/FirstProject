#include "cfd.h"

using namespace CFD;

SolverType CFD::convertSolverType(const std::string& solver) {
    if (solver == "jacobi") {
        return SolverType::JACOBI;
    }
    else if (solver == "multigrid_jacobi") {
        return SolverType::MULTIGRID_JACOBI;
    }
    else if (solver == "conjugated_gradient") {
        return SolverType::CONJUGATED_GRADIENT;
    }
    else if (solver == "multigrid_pcg") {
        return SolverType::MULTIGRID_PCG;
    }
    else {
        throw std::invalid_argument("Invalid solver type");
    }
}

FluidParams::FluidParams(std::string name, int argc, char* argv[]) 
    : argument_parser(name)
{
    this->argument_parser.add_argument("-i", "--imax").help("imax").default_value(this->imax).action([](const std::string& value) { return std::stoi(value); });
    this->argument_parser.add_argument("-j", "--jmax").help("jmax").default_value(this->jmax).action([](const std::string& value) { return std::stoi(value); });
    this->argument_parser.add_argument("-x", "--xlength").help("xlength").default_value(this->xlength).action([](const std::string& value) { return std::stod(value); });
    this->argument_parser.add_argument("-y", "--ylength").help("ylength").default_value(this->ylength).action([](const std::string& value) { return std::stod(value); });
    this->argument_parser.add_argument("-z", "--t_end").help("t_end").default_value(this->t_end).action([](const std::string& value) { return std::stod(value); });
    this->argument_parser.add_argument("-u", "--tau").help("tau").default_value(this->tau).action([](const std::string& value) { return std::stod(value); });
    this->argument_parser.add_argument("-e", "--eps").help("eps").default_value(this->eps).action([](const std::string& value) { return std::stod(value); });
    this->argument_parser.add_argument("-o", "--omg").help("omg").default_value(this->omg).action([](const std::string& value) { return std::stod(value); });
    this->argument_parser.add_argument("-a", "--alpha").help("alpha").default_value(this->alpha).action([](const std::string& value) { return std::stod(value); });
    this->argument_parser.add_argument("-r", "--Re").help("Re").default_value(this->Re).action([](const std::string& value) { return std::stod(value); });
    this->argument_parser.add_argument("-t", "--t").help("t").default_value(this->t).action([](const std::string& value) { return std::stod(value); });
    this->argument_parser.add_argument("-d", "--dt").help("dt").default_value(this->dt).action([](const std::string& value) { return std::stod(value); });
    this->argument_parser.add_argument("-l", "--save_interval").help("VTK save interval").default_value(this->save_interval).action([](const std::string& value) { return std::stod(value); });
    this->argument_parser.add_argument("-s", "--solver").help("solver").default_value("jacobi").action([](const std::string& value) { return value; });


    try {
        this->argument_parser.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << this->argument_parser;
        exit(1);
    }

    this->imax = this->argument_parser.get<int>("--imax"),
    this->jmax = this->argument_parser.get<int>("--jmax"),
    this->xlength = this->argument_parser.get<double>("--xlength"),
    this->ylength = this->argument_parser.get<double>("--ylength"),
    this->t_end = this->argument_parser.get<double>("--t_end"),
    this->tau = this->argument_parser.get<double>("--tau"),
    this->eps = this->argument_parser.get<double>("--eps"),
    this->omg = this->argument_parser.get<double>("--omg"),
    this->alpha = this->argument_parser.get<double>("--alpha"),
    this->Re = this->argument_parser.get<double>("--Re"),
    this->t = this->argument_parser.get<double>("--t"),
    this->dt = this->argument_parser.get<double>("--dt"),
    this->save_interval = this->argument_parser.get<double>("--save_interval"),
    this->solver_type = convertSolverType(this->argument_parser.get<std::string>("--solver"));
}