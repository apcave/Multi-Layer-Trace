cmake_minimum_required(VERSION 3.14)
project(MultiLayerTrace LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add source files
set(SRC_FILES
    Wave.cpp
    Medium.cpp
    Material.cpp
    Interface.cpp
    Composite.cpp
    PythonInterface.cpp
)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/src)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
execute_process(
    COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
    OUTPUT_VARIABLE pybind11_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
list(APPEND CMAKE_PREFIX_PATH "${pybind11_dir}")
find_package(pybind11 REQUIRED)

# Python module target
pybind11_add_module(composite_wave_response PythonInterface.cpp Wave.cpp Medium.cpp Material.cpp Interface.cpp Composite.cpp)
