project(kernel)

set(KERNEL_SOURCES
    src/timer.cpp
)

add_library(kernel STATIC ${KERNEL_SOURCES})

target_compile_options(kernel PRIVATE -fPIC)

target_include_directories(kernel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)