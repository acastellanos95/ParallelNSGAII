cmake_minimum_required(VERSION 3.22)
project(NSGAII)

set(CMAKE_CXX_STANDARD 17)

add_executable(NSGAII main.cpp Individual.h Rand.h Utils.h Problem.h Sort.h NSGA.cpp NSGA.h lib/WFG/WFG.cpp lib/WFG/WFG.h lib/hypervolume/HypervolumeUtils.h)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(NSGAII PUBLIC OpenMP::OpenMP_CXX)
endif()
