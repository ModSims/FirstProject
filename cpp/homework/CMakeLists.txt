project(homework)

find_package(Eigen3 REQUIRED COMPONENTS CXX)

# Compile CFD Programs 
add_executable(hausaufgabe_2 prog/hausaufgabe_2.cpp)
target_compile_options(hausaufgabe_2 PRIVATE -fPIC)
target_include_directories(hausaufgabe_2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(hausaufgabe_2 Eigen3::Eigen kernel maths solvers)