# Includes
target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Build lib
add_library(geometry STATIC sphere.cpp)