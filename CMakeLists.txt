cmake_minimum_required(VERSION 3.9)

project(fortperc Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP)

find_package(LAPACK REQUIRED)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all -Wall")
endif()

file(GLOB fortperc_src ${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90)

add_library(fortperc ${fortperc_src})

if(OpenMP_Fortran_FOUND)
    target_link_libraries(fortperc OpenMP::OpenMP_Fortran)
else()
    message(WARNING "No OpenMP support found!")
endif()

target_link_libraries(fortperc ${LAPACK_LIBRARIES})
