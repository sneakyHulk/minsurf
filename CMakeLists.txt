set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

PROJECT(minsurf)

cmake_minimum_required(VERSION 3.8)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(minsurf MACOSX_BUNDLE src/main.cpp)

SET(GCC_COVERAGE_COMPILE_FLAGS "-std=c++11 -O3 -fopenmp")
SET(GCC_COVERAGE_LINK_FLAGS    "-fopenmp")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )


target_link_libraries(minsurf ${VTK_LIBRARIES})


