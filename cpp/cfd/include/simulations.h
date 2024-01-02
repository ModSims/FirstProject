#pragma once
#include "cfd.h"

namespace CFD {
    using namespace Eigen;

    class LidDrivenCavity2D : public FluidSimulation {
        public:
            LidDrivenCavity2D(const FluidParams& params) : FluidSimulation(params) {
                grid = StaggeredGrid(imax, jmax, xlength, ylength);
            }
            void setBoundaryConditionsU() override;
            void setBoundaryConditionsV() override;
            void setBoundaryConditionsP() override;
            void setBoundaryConditionsVelocityGeometry() override {};
            void setBoundaryConditionsInterpolatedVelocityGeometry() override {};
            void setBoundaryConditionsPGeometry() override {};
    };

    class FlowOverStep2D : public FluidSimulation {
        public:
            FlowOverStep2D(const FluidParams& params) : FluidSimulation(params) {
                grid = StaggeredGrid(imax, jmax, xlength, ylength);
            }
            void setBoundaryConditionsU() override;
            void setBoundaryConditionsV() override;
            void setBoundaryConditionsP() override;
            void run();
    };

    class KarmanVortexStreet2D : public FluidSimulation {
        public:
            KarmanVortexStreet2D(const FluidParams& params) : FluidSimulation(params) {
                grid = StaggeredGrid(imax, jmax, xlength, ylength);
            }
            void setBoundaryConditionsU() override;
            void setBoundaryConditionsV() override;
            void setBoundaryConditionsP() override;
            void run();
    };

    class PlaneShearFlow2D : public FluidSimulation {
        public:
            PlaneShearFlow2D(const FluidParams& params) : FluidSimulation(params) {
                grid = StaggeredGrid(imax, jmax, xlength, ylength);
            }
            void setBoundaryConditionsU() override;
            void setBoundaryConditionsV() override;
            void setBoundaryConditionsP() override;
            void setBoundaryConditionsVelocityGeometry() override {};
            void setBoundaryConditionsInterpolatedVelocityGeometry() override {};
            void setBoundaryConditionsPGeometry() override {};
    };
}
