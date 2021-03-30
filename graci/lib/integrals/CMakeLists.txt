#
# Minimum version of CMake required
#
cmake_minimum_required(VERSION 3.2)

#
# Project name and language
#
project(int_pyscf)
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
#    set(dialect "-O2 -ffree-form -fimplicit-none -mavx512f")
    set(dialect "-O2 -ffree-form -cpp -fimplicit-none -mavx2")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-g -O2 -cpp -stand f08 -free -implicitnone -traceback \
        -march=corei7-avx")
    set(bounds "-check bounds")
endif()

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

#
# HDF5 include and lib directories
#
if (HDF5_INC_DIR)
  include_directories(${HDF5_INC_DIR})
endif()
if (HDF5_LIB_DIR)
  link_directories(${HDF5_LIB_DIR})
endif()

#
# Set the modules directory
#
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/modules)

#
# Add the library target
#
file(GLOB_RECURSE sources  src/*.f90 src/*.F90)
add_library(int_pyscf SHARED ${sources})

# Link the target to the HDF5 Fortran library
target_link_libraries(int_pyscf -lhdf5_fortran)

#
# Install
#
install(TARGETS int_pyscf DESTINATION lib)
