#include "cfd.h"
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>

void CFD::saveVTK(FluidSimulation* sim) {
    vtkSmartPointer<vtkStructuredGrid> vtk_grid = vtkSmartPointer<vtkStructuredGrid>::New();
    vtk_grid->SetDimensions(sim->imax, sim->jmax, 1);

    // Create vtk points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int j = 0; j < sim->jmax; ++j) {
        for (int i = 0; i < sim->imax; ++i) {
            points->InsertNextPoint(i * sim->grid.dx, j * sim->grid.dy, 0);
        }
    }

    vtk_grid->SetPoints(points);

    // Add pressure data to the grid
    vtkSmartPointer<vtkDoubleArray> pressure_array = vtkSmartPointer<vtkDoubleArray>::New();
    pressure_array->SetNumberOfComponents(1);
    pressure_array->SetName("Pressure");

    // Add velocity data to the grid as a vector field
    vtkSmartPointer<vtkDoubleArray> velocity_array = vtkSmartPointer<vtkDoubleArray>::New();
    velocity_array->SetNumberOfComponents(3);  // 3 components for the vector (VelocityX, VelocityY, 0)
    velocity_array->SetComponentName(0, "VelocityX");
    velocity_array->SetComponentName(1, "VelocityY");
    velocity_array->SetComponentName(2, "VelocityZ");  // Z component is set to 0
    velocity_array->SetName("Velocity");

    for (int j = 0; j < sim->jmax; ++j) {
        for (int i = 0; i < sim->imax; ++i) {
            pressure_array->InsertNextValue(sim->grid.p.coeffRef(i, j));
            velocity_array->InsertNextTuple3(sim->grid.u_interpolated.coeffRef(i, j),
                                             sim->grid.v_interpolated.coeffRef(i, j),
                                             0.0);  // Z component is set to 0
        }
    }

    vtk_grid->GetPointData()->AddArray(pressure_array);
    vtk_grid->GetPointData()->AddArray(velocity_array);

    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    char filename[100];
    sprintf(filename, "output_%d.vts", sim->it);
    writer->SetFileName(filename);
    writer->SetInputData(vtk_grid);
    writer->Write();
}

void CFD::saveVTKGeometry(FluidSimulation* sim) {
    vtkSmartPointer<vtkStructuredGrid> vtk_grid = vtkSmartPointer<vtkStructuredGrid>::New();
    vtk_grid->SetDimensions(sim->imax, sim->jmax, 1);

    // Create vtk points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int j = 0; j < sim->jmax; ++j) {
        for (int i = 0; i < sim->imax; ++i) {
            points->InsertNextPoint(i * sim->grid.dx, j * sim->grid.dy, 0);
        }
    }

    vtk_grid->SetPoints(points);

    // Add geometry
    vtkSmartPointer<vtkDoubleArray> geometry_array = vtkSmartPointer<vtkDoubleArray>::New();
    geometry_array->SetNumberOfComponents(1);
    geometry_array->SetName("Geometry");

    bool has_geometry = false;

    for (int j = 0; j < sim->jmax; ++j) {
        for (int i = 0; i < sim->imax; ++i) {
            if ((sim->grid.flag_field(i, j) & FlagFieldMask::MASK_CELL_TYPE) == (FlagFieldMask::CELL_OBSTACLE & FlagFieldMask::MASK_CELL_TYPE)) {
                geometry_array->InsertNextValue(1.0);
                has_geometry = true;
            } else {
                geometry_array->InsertNextValue(0.0);
            
            }
        }
    }

    if (has_geometry) {
        vtk_grid->GetPointData()->AddArray(geometry_array);

        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        char filename[100];
        sprintf(filename, "geometry.vts");
        writer->SetFileName(filename);
        writer->SetInputData(vtk_grid);
        writer->Write();
    }
}