# CMakeList.txt: CMake-Projekt für "minsurface". Schließen Sie die Quelle ein, und definieren Sie
# projektspezifische Logik hier.
#
cmake_minimum_required (VERSION 3.8)

# Fügen Sie der ausführbaren Datei für dieses Projekt eine Quelle hinzu.


# TODO: Fügen Sie bei Bedarf Tests und Installationsziele hinzu.
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories("/usr/include/eigen3")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable (minsurface main.cpp)

SET(GCC_COVERAGE_COMPILE_FLAGS "-march=native -O3 -fopenmp -frename-registers -funroll-loops")
SET(GCC_COVERAGE_LINK_FLAGS    "-fopenmp")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )


target_link_libraries(minsurface ${VTK_LIBRARIES})