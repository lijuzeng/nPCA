cmake_minimum_required(VERSION 3.22)

project(npca VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 11)

add_compile_options(-O3 -Wall -march=native)

add_executable(npca main.cpp npca.hpp)
