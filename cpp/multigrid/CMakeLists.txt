project(multigrid)

find_package(Eigen3 REQUIRED COMPONENTS CXX)

set(MULTIGRID_SOURCES
    src/multigrid.cpp
)

# Compile the Multigrid library
add_library(multigrid STATIC ${MULTIGRID_SOURCES})
target_compile_options(multigrid PRIVATE -fPIC)
target_include_directories(multigrid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(multigrid Eigen3::Eigen kernel solvers cfd)

# Compile Multigrid Programs 
add_executable(multigrid_benchmark prog/benchmark.cpp)
target_compile_options(multigrid_benchmark PRIVATE -fPIC)
target_include_directories(multigrid_benchmark PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(multigrid_benchmark Eigen3::Eigen kernel maths solvers cfd multigrid)