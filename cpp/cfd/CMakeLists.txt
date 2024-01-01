project(cfd)

find_package(Eigen3 REQUIRED COMPONENTS CXX)
find_package(VTK REQUIRED COMPONENTS
    CommonDataModel
    CommonCore
    CommonExecutionModel
    CommonMath
    CommonMisc
    CommonSystem
    CommonTransforms
    IOCore
    IOLegacy
    IOParallelXML
    IOXML
    IOXMLParser
    ParallelCore
    FiltersCore
    vtksys
)

set(CFD_SOURCES
    src/cfd.cpp
    src/vtk.cpp
    src/lid_driven_cavity_2d.cpp
    src/flow_over_step_2d.cpp
    src/karman_vortex_street_2d.cpp
    src/plane_shear_flow_2d.cpp
)

# Compile the CFD library
add_library(cfd STATIC ${CFD_SOURCES})
target_compile_options(cfd PRIVATE)
target_include_directories(cfd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cfd Eigen3::Eigen kernel solvers ${VTK_LIBRARIES})

# Compile CFD Programs 
## Lid Driven Cavity 2D
add_executable(lid_driven_cavity_2d prog/lid_driven_cavity_2d.cpp)
target_compile_options(lid_driven_cavity_2d PRIVATE -fPIC)
target_include_directories(lid_driven_cavity_2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(lid_driven_cavity_2d cfd Eigen3::Eigen kernel)

# Flow over step 2D
add_executable(flow_over_step_2d prog/flow_over_step_2d.cpp)
target_compile_options(flow_over_step_2d PRIVATE -fPIC)
target_include_directories(flow_over_step_2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(flow_over_step_2d cfd Eigen3::Eigen kernel)

# Karman Vortex Street 2D
add_executable(karman_vortex_street_2d prog/karman_vortex_street_2d.cpp)
target_compile_options(karman_vortex_street_2d PRIVATE -fPIC)
target_include_directories(karman_vortex_street_2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(karman_vortex_street_2d cfd Eigen3::Eigen kernel)

# Plane Shear Flow 2D
add_executable(plane_shear_flow_2d prog/plane_shear_flow_2d.cpp)
target_compile_options(plane_shear_flow_2d PRIVATE -fPIC)
target_include_directories(plane_shear_flow_2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(plane_shear_flow_2d cfd Eigen3::Eigen kernel)
