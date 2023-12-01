project(cfd)

find_package(Eigen3 REQUIRED COMPONENTS CXX)

set(CFD_SOURCES
    src/cfd.cpp
)

# Compile the CFD library
add_library(cfd STATIC ${CFD_SOURCES})
target_compile_options(cfd PRIVATE -fPIC)
target_include_directories(cfd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cfd Eigen3::Eigen)

# Compile CFD Programs 
add_executable(lid_driven_cavity_2d prog/lid_driven_cavity_2d.cpp)
target_compile_options(lid_driven_cavity_2d PRIVATE -fPIC)
target_include_directories(lid_driven_cavity_2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(lid_driven_cavity_2d cfd Eigen3::Eigen kernel)