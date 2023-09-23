project(geometry)

set(GEOMETRY_SOURCES
    src/sphere.cpp
)

add_library(geometry STATIC ${GEOMETRY_SOURCES})

target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
