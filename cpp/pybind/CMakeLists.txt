project(PyModSims)

# Include external libraries
find_package(Eigen3 REQUIRED COMPONENTS CXX)

file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Include pybind11
find_package(pybind11 REQUIRED)
pybind11_add_module(${PROJECT_NAME} ${SOURCES})

target_compile_options(${PROJECT_NAME} PRIVATE -fPIC) # Add this line to compile with -fPIC

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link Eigen3 to the python library
target_link_libraries(${PROJECT_NAME} 
    PRIVATE Eigen3::Eigen
    PRIVATE kernel
    PRIVATE solvers
    PRIVATE math
)

# Add pytest test cases
add_test(NAME python_tests COMMAND pytest ${CMAKE_SOURCE_DIR}/pybind/tests)
