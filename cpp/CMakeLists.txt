cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(FirstProject VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_subdirectory(geometry)
add_subdirectory(solvers)
add_subdirectory(kernel)
