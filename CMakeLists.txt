cmake_minimum_required(VERSION 3.9)

project(percfort)
enable_language(Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(generalflags "-Wall -Wextra")
set(speedflags "-Ofast -fopenmp")
set(debugflags "-O0 -g -fcheck=all")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${generalflags}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE} ${speedflags}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG} ${debugflags}")

add_library(utilities ${CMAKE_CURRENT_SOURCE_DIR}/src/utilities.f90)
add_library(percolation ${CMAKE_CURRENT_SOURCE_DIR}/src/percolation.f90)
add_library(hk ${CMAKE_CURRENT_SOURCE_DIR}/src/hk.f90)
add_library(randomwalk ${CMAKE_CURRENT_SOURCE_DIR}/src/randomwalk.f90)

target_link_libraries(percolation hk utilities)
target_link_libraries(randomwalk utilities percolation)

set(mylibraries utilities percolation lapack hk randomwalk)

