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
# Profile-guided optimisation (ifort only)
#
option(PGOGEN "Profile-guided optimisation generation" OFF)
option(PGOUSE "Profile-guided optimisation use" OFF)
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  if (PGOGEN)
    add_definitions(-prof-gen -prof-dir=../)
  elseif (PGOUSE)
    add_definitions(-prof-use -prof-dir=../)
  endif()
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
  set(dialect "-g -O0 -fcheck=all -fbacktrace -cpp -ffree-form -fimplicit-none \
        -mavx -lblas -llapack")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")  
  set(dialect "-g -O3 -auto -ipo -cpp -free -implicitnone -traceback \
        -march=native -mkl -w -qopenmp")
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
# Add bitCI library target
#
file(GLOB_RECURSE CIsources
  src/ci/* src/include/* src/io/* src/utils/* src/spin_coupling/*
  src/slater_condon/*)
add_library(bitci SHARED ${CIsources})

#
# Add bitSI library target
#
file(GLOB_RECURSE SIsources
  src/si/* src/include/* src/io/* src/utils/* src/spin_coupling/*
  src/slater_condon/*)
add_library(bitsi SHARED ${SIsources})


# Link to int_pyscf
target_link_libraries(bitci -lint_pyscf)

#
# Install
#
install(TARGETS bitci bitsi DESTINATION lib)
