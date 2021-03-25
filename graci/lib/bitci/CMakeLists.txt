#
# Minimum version of CMake required
#
cmake_minimum_required(VERSION 3.2)

#
# Project name and language
#
project(bitci)
enable_language(Fortran)

#
# Conditional compilation flags
#
option(CBINDING "Optional C binding for public routines" OFF)
if (CBINDING)
  add_definitions(-DCBINDING)
endif()

#
# Compiler-specific flags
#
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-O0 -cpp -ffree-form -fimplicit-none")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-g -O2 -cpp -stand f08 -free -implicitnone -traceback")
    set(bounds "-check bounds")
endif()

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

#
# Compile
#
file(GLOB_RECURSE sources  src/*.f90 src/*.F90)
add_library(bitci SHARED ${sources})

#
# Install
#
install(TARGETS bitci DESTINATION lib)
