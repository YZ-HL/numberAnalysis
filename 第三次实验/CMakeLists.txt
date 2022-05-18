cmake_minimum_required(VERSION 3.21)
project(NumberAnalysisThird)

set(CMAKE_CXX_STANDARD 14)

add_executable(NumberAnalysisThird main.cpp gauss.h matrix.h utils.h decomposition.h)
