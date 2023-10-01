cmake_minimum_required(VERSION 3.25)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
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

# Include libraries
add_subdirectory(solvers)
add_subdirectory(pybind)
