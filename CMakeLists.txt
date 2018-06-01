cmake_minimum_required(VERSION 3.10)

project(percfort)
enable_language(Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(generalflags "-Wall -Wextra")
set(speedflags "-Ofast -fopenmp")
set(debugflags "-O0 -g  -fcheck=all")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${generalflags}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE} ${speedflags}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG} ${debugflags}")

add_library(utilities ${CMAKE_SOURCE_DIR}/src/utilities.f90)
add_library(percolation ${CMAKE_SOURCE_DIR}/src/percolation.f90)
add_library(hk ${CMAKE_SOURCE_DIR}/src/hk.f90)

add_dependencies(percolation utilities hk)

set(mylibraries utilities percolation lapack hk)

