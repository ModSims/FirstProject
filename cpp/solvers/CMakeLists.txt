project(solvers)

set(SOLVERS_SOURCES
    src/trivial.cpp
    src/jacobi.cpp
    src/gauss_seidel.cpp
)

add_library(solvers STATIC ${SOLVERS_SOURCES})

target_include_directories(solvers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
