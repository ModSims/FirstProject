project(kernel)

set(KERNEL_SOURCES
    src/timer.cpp
    src/filehandler.cpp
)

# External packages
find_package(Eigen3 REQUIRED COMPONENTS CXX)

add_library(kernel STATIC ${KERNEL_SOURCES})

target_compile_options(kernel PRIVATE -fPIC)

target_include_directories(kernel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(kernel PRIVATE Eigen3::Eigen)