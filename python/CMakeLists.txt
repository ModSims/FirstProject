cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(ModSims VERSION 1.0.0 LANGUAGES CXX)

# Set the binary and runtime output directories
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")

# Include external libraries
find_package(Eigen3 REQUIRED COMPONENTS CXX)

# Include pybind11
find_package(pybind11 REQUIRED)
pybind11_add_module(${PROJECT_NAME} src/main.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../cpp/solvers/include)

# Link Eigen3 to the python library
target_link_libraries(${PROJECT_NAME} 
    PRIVATE Eigen3::Eigen
    PRIVATE "${CMAKE_SOURCE_DIR}/../cpp/build/solvers/libsolvers.a"
)