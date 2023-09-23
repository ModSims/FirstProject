project(kernel)

# External packages
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Build lib
set(KERNEL_SOURCES
    src/kernel.cpp
)
add_executable(kernel ${KERNEL_SOURCES})
target_link_libraries(kernel
                        PRIVATE geometry 
                        PRIVATE solvers 
                        PRIVATE MPI::MPI_CXX 
                        PRIVATE ${HDF5_LIBRARIES}
)