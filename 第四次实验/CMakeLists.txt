cmake_minimum_required(VERSION 3.21)
project(FG)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

add_executable(FG main.cpp matrix.h iterative.h utils.h)
