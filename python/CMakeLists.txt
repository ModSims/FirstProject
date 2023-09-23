cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(ModSims VERSION 1.0.0 LANGUAGES CXX)

# Set the binary and runtime output directories
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}/python/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")

# Include pybind11
find_package(pybind11 REQUIRED)
pybind11_add_module(modsims src/modsims.cpp)