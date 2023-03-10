cmake_minimum_required(VERSION 3.21)
project(compare)

set(CMAKE_CXX_STANDARD 11)

add_compile_options(-O3 -Wall -march=native)

add_executable(compare main.cpp net.hpp wyhash.h)
