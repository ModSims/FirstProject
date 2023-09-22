cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Capricorn-CFD)

add_subdirectory(src/geometry)
add_subdirectory(src/kernel)

add_executable(Capricorn-CFD main.cpp)
target_link_libraries(Capricorn-CFD PRIVATE geometry PRIVATE kernel PRIVATE MPI::MPI_CXX PRIVATE ${HDF5_LIBRARIES})
set_target_properties(Capricorn-CFD PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PROJECT_SOURCE_DIR}/bin/debug
                      RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROJECT_SOURCE_DIR}/bin/release
)