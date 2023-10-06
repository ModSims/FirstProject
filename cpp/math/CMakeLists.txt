# CMake minimum version
cmake_minimum_required(VERSION 3.25)

# Project name
project(math)

# Define the source files for the library
set(MATH_SOURCES
    src/radians.cpp
)

# External packages
find_package(Eigen3 REQUIRED COMPONENTS CXX)
find_package(Catch2 3 REQUIRED)

# Create a library from the source files
add_library(math STATIC ${MATH_SOURCES})

target_compile_options(math PRIVATE -fPIC) # Add this line to compile with -fPIC

# Set include directories for the library
target_include_directories(math PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link Eigen3 to the solvers library
target_link_libraries(math PRIVATE Eigen3::Eigen)

# List of executable targets and their respective source files
set(EXECUTABLES
    radians
)

include(CTest)
include(Catch)
foreach(EXECUTABLE ${EXECUTABLES})
    add_executable(test_math_${EXECUTABLE} tests/${EXECUTABLE}.cpp)
    target_link_libraries(test_math_${EXECUTABLE}
                          PRIVATE math
                          PRIVATE Eigen3::Eigen
                          PRIVATE Catch2::Catch2WithMain
    )
    catch_discover_tests(test_math_${EXECUTABLE})
endforeach()