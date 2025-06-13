cmake_minimum_required(VERSION 3.14)
project(MultiLayerTrace LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add source files
set(SRC_FILES
    src/Wave.cpp
    src/Medium.cpp
    src/Material.cpp
    src/Interface.cpp
    src/Composite.cpp
    src/PythonInterface.cpp
)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/src)

# pybind11
find_package(pybind11 REQUIRED)

# Python module target
pybind11_add_module(composite_wave_response src/python_interface.cpp src/Wave.cpp src/Medium.cpp src/Material.cpp src/Interface.cpp src/Composite.cpp)
