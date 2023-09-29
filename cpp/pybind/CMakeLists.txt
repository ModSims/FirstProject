cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(
  CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -Wno-register -Wno-comment -fPIC \
  -march=native -Wfatal-errors -fvisibility=hidden"
)

project(ModSims VERSION 1.0.0 LANGUAGES CXX)

# Set the binary and runtime output directories
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")

# Include external libraries
find_package(Eigen3 REQUIRED COMPONENTS CXX)

# Include pybind11
find_package(pybind11 REQUIRED)
pybind11_add_module(${PROJECT_NAME} src/main.cpp)

# Link Eigen3 to the python library
target_link_libraries(${PROJECT_NAME} 
    PRIVATE Eigen3::Eigen
    PRIVATE solvers
)