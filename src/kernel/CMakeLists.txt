# Includes
target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)


# External packages
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(HDF5 REQUIRED COMPONENTS CXX)


# Build lib
add_library(kernel STATIC kernel.cpp)