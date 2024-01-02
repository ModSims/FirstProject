cmake_minimum_required(VERSION 3.25)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Debug)
set(
  CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -Wno-register -Wno-comment -fPIC \
  -march=native -Wfatal-errors -fvisibility=hidden"
)
set(CMAKE_CTEST_ARGUMENTS "--output-on-failure")

project(ModSims VERSION 1.0.0 LANGUAGES CXX)

# Set the binary and runtime output directories
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}/cpp/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../out/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../out/lib")

enable_testing()

# Add extern projects
set(VTK_GROUP_ENABLE_StandAlone DONT_WANT)
set(VTK_GROUP_ENABLE_Rendering DONT_WANT)
set(VTK_MODULE_ENABLE_VTK_CommonDataModel YES)
set(VTK_MODULE_ENABLE_VTK_CommonCore YES)
set(VTK_MODULE_ENABLE_VTK_CommonExecutionModel YES)
set(VTK_MODULE_ENABLE_VTK_CommonMath YES)
set(VTK_MODULE_ENABLE_VTK_CommonMisc YES)
set(VTK_MODULE_ENABLE_VTK_CommonSystem YES)
set(VTK_MODULE_ENABLE_VTK_CommonTransforms YES)
set(VTK_MODULE_ENABLE_VTK_IOCore YES)
set(VTK_MODULE_ENABLE_VTK_IOLegacy YES)
set(VTK_MODULE_ENABLE_VTK_IOParallelXML YES)
set(VTK_MODULE_ENABLE_VTK_IOXML YES)
set(VTK_MODULE_ENABLE_VTK_IOXMLParser YES)
set(VTK_MODULE_ENABLE_VTK_ParallelCore YES)
set(VTK_MODULE_ENABLE_VTK_FiltersCore YES)
set(VTK_MODULE_ENABLE_VTK_vtksysYES)
add_subdirectory(extern/vtk)

# Include libraries
add_subdirectory(kernel)
add_subdirectory(solvers)
add_subdirectory(maths)
add_subdirectory(cfd)
add_subdirectory(homework)
# add_subdirectory(pybind)