# CMake minimum version
cmake_minimum_required(VERSION 3.10)

# Project name
project(solvers)

# Define the source files for the library
set(SOLVERS_SOURCES
    src/linear.cpp
    src/trivial.cpp
    src/jacobi.cpp
    src/gauss_seidel.cpp
)

# External packages
find_package(Eigen3 REQUIRED COMPONENTS CXX)

# Create a library from the source files
add_library(solvers STATIC ${SOLVERS_SOURCES})

target_compile_options(solvers PRIVATE -fPIC) # Add this line to compile with -fPIC

# Set include directories for the library
target_include_directories(solvers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link Eigen3 to the solvers library
target_link_libraries(solvers PRIVATE Eigen3::Eigen)

# List of executable targets and their respective source files
set(EXECUTABLES
    trivial
    jacobi
    gauss_seidel
)

foreach(EXECUTABLE ${EXECUTABLES})
    add_executable(test_solver_${EXECUTABLE} tests/${EXECUTABLE}.cpp)
    target_link_libraries(test_solver_${EXECUTABLE} PRIVATE solvers PRIVATE Eigen3::Eigen)
endforeach()
