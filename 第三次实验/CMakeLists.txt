cmake_minimum_required(VERSION 3.21)
project(NumberAnalysisThird)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

add_executable(NumberAnalysisThird main.cpp gauss.h matrix.h utils.h decomposition.h)
