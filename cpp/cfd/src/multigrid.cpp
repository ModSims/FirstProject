#include "multigrid.h"

using namespace CFD;

void Multigrid::vcycle(MultigridHierarchy *hierarchy, int currentLevel, double omg) {
    // Briggs, Multigrid Tutorial, p. 41
    // reusing p and RHS for error and residual for simplicity and memory efficiency
    // (Ae = res) and (Ap = RHS)

    if (currentLevel == 0) {
        // Relax on the coarset grid
        relax(hierarchy->grids[0].get(), 31, omg);
    } else {
        // Relax on the current grid
        relax(hierarchy->grids[currentLevel].get(), 1, omg);

        // Restrict the residual to the coarser grid
        restrict_operator(hierarchy->grids[currentLevel].get(), hierarchy->grids[currentLevel-1].get());

        // Set the error on the coarser grid to zero
        hierarchy->grids[currentLevel-1].get()->p.setZero();

        // Recursively call vcycle
        vcycle(hierarchy, currentLevel-1, omg);

        // Prolongate the error to the finer grid
        prolongate_operator(hierarchy->grids[currentLevel-1].get(), hierarchy->grids[currentLevel].get());

        // correct the current approximation with the error on the current grid
        hierarchy->grids[currentLevel].get()->p += hierarchy->grids[currentLevel].get()->res;

        // Post-smooth on the current grid
        relax(hierarchy->grids[currentLevel].get(), 1, omg);
    }
}

void Multigrid::restrict_operator(const StaggeredGrid *fine, StaggeredGrid *coarse) {
    // Restrict with full weighting
    // Briggs, Multigrid Tutorial, p. 36

    // Restrict res^h to RHS^{2h} but saving on RHS^{2h}
    for (int i = 1; i <= coarse->imax; i++) {
        for (int j = 1; j <= coarse->jmax; j++) {
            coarse->RHS(i,j) = 0.0625*(
                fine->res(2*i-1,2*j-1) + fine->res(2*i-1,2*j+1) + fine->res(2*i+1,2*j-1) + fine->res(2*i+1,2*j+1)
                + 2 * (fine->res(2*i,2*j-1) + fine->res(2*i,2*j+1) + fine->res(2*i-1,2*j) + fine->res(2*i+1,2*j))
                + 4 * fine->res(2*i,2*j)
            );
        }
    }
}

void Multigrid::prolongate_operator(const StaggeredGrid *coarse, StaggeredGrid *fine) {
    // Prolongate with linear interpolation
    // Briggs, Multigrid Tutorial, p. 35

    // Prolongate p^{2h} to res^h
    fine->res.setZero();
    for (int i = 0; i <= coarse->imax; i++) {
        for (int j = 0; j <= coarse->jmax; j++) {
            fine->res(2*i,2*j) = coarse->p(i,j);
            fine->res(2*i+1,2*j) = 0.5 * (coarse->p(i,j) + coarse->p(i+1,j));
            fine->res(2*i,2*j+1) = 0.5 * (coarse->p(i,j) + coarse->p(i,j+1));
            fine->res(2*i+1,2*j+1) = 0.25 * (coarse->p(i,j) + coarse->p(i+1,j) + coarse->p(i,j+1) + coarse->p(i+1,j+1));
        }
    }
}

void Multigrid::relax(StaggeredGrid *grid, int numSweeps, double omg) {
    // Relaxation with Jacobi
    // Its in abstract form, so it can be used for (Ae = res) or (Ap = RHS)
    
    for (int sweep = 0; sweep < numSweeps; sweep++) {
        // Jacobi smoother with relaxation factor (omega)
        for (int i = 1; i <= grid->imax; i++) {
            for (int j = 1; j <= grid->jmax; j++) {
                grid->p(i, j) = (1 - omg) * grid->p(i, j) +
                                    omg * 0.25 * (
                                        grid->p(i - 1, j) + grid->p(i + 1, j) +
                                        grid->p(i, j - 1) + grid->p(i, j + 1)
                                        - grid->dxdy * (grid->RHS(i, j))
                                    ) / (1 + 2 * (grid->dx2 + grid->dy2));
            }
        }
    }
}