cmake_minimum_required(VERSION 3.8)
project(Assignment_01_Histogram)

include_directories(/usr/include/openmpi-x86_64)
set(CMAKE_C_COMPILER mpiCC)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_COMPILER mpic++)


set(SOURCE_FILES main.cpp)
add_executable(Assignment_01_Histogram ${SOURCE_FILES})