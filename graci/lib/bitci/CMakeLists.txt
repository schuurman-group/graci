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
# HDF5 include and lib directories
#
if (HDF5_INC_DIR)
  include_directories(${HDF5_INC_DIR})
endif()
if (HDF5_LIB_DIR)
  link_directories(${HDF5_LIB_DIR})
endif()

#
# Compiler-specific flags
#
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-g -O2 -ffast-math -fbacktrace -cpp -ffree-form -fimplicit-none \
        -mavx -lblas -llapack")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-g -O3 -cpp -free -implicitnone -traceback \
        -xCORE-AVX2 -mkl -diag-disable 7423 -qopenmp")
  set(bounds "-check bounds")
endif()

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

#
# int_pyscf include and lib directories
#
set(INT_PYSCF_INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../integrals/modules)
include_directories(${INT_PYSCF_INC_DIR})
set(INT_PYSCF_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../integrals/lib)
link_directories(${INT_PYSCF_LIB_DIR})

#
# Add the library target 
#
file(GLOB_RECURSE sources src/*)
add_library(bitci SHARED ${sources})

# Link to int_pyscf
target_link_libraries(bitci -lint_pyscf)

#
# Install
#
install(TARGETS bitci DESTINATION lib)
