#pragma once

#include <vector>
#include <iostream>
#include <memory>
#include "staggered_grid.h"

namespace CFD {
    using namespace Eigen;

    class MultigridHierarchy {
    public:
        MultigridHierarchy(int numLevels, const StaggeredGrid* finestGrid) {
            grids.reserve(numLevels + 1);

            for (int i = numLevels; i >= 1; i--) {
                int imax = finestGrid->imax / (1 << i);
                int jmax = finestGrid->jmax / (1 << i);
                std::cout << "Creating grid with imax = " << imax << " and jmax = " << jmax << std::endl;
                grids.push_back(std::make_unique<StaggeredGrid>(imax, jmax, finestGrid->xlength, finestGrid->ylength));
            }

            grids.push_back(std::unique_ptr<StaggeredGrid>(const_cast<StaggeredGrid*>(finestGrid)));
        }

        int numLevels() const {
            return grids.size();
        }

        std::vector<std::unique_ptr<StaggeredGrid>> grids;
    };

    namespace Multigrid {
        void vcycle(MultigridHierarchy *hierarchy, int currentLevel, double omg);
        void restrict_operator(const StaggeredGrid* fine, StaggeredGrid* coarse);
        void prolongate_operator(const StaggeredGrid* coarse, StaggeredGrid* fine);
        void relax(StaggeredGrid* grid, int numSweeps, double omg);
    }
}
