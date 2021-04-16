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
option(PGO_GEN "Profile-guided optimisation generation" OFF)
option(PGO_USE "Profile-guided optimisation use" OFF)
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  if (PGO_GEN)
    add_definitions(-prof-gen -prof-dir=../)
  elseif (PGO_USE)
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
  set(dialect "-g -O2 -ffast-math -fbacktrace -cpp -ffree-form -fimplicit-none \
        -mavx -lblas -llapack")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")  
  set(dialect "-w -g -O3 -cpp -free -implicitnone -traceback -mkl -qopenmp \
    -march=native -auto -ipo -qopt-prefetch")
  set(bounds "-check bounds")
endif()

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

#
# Add bitCI library target
#
file(GLOB_RECURSE CIsources
  src/ci/* src/include/* src/io/* src/utils/* src/spin_coupling/*
  src/slater_condon/* ../integrals/src/*)
add_library(bitci SHARED ${CIsources})

#
# Add bitSI library target
#
file(GLOB_RECURSE SIsources
  src/si/* src/include/* src/io/* src/utils/* src/spin_coupling/*
  src/slater_condon/* ../integrals/src/*)
add_library(bitsi SHARED ${SIsources})

# Link the target to the HDF5 Fortran library
target_link_libraries(bitci -lhdf5_fortran)

#
# Install
#
install(TARGETS bitci bitsi DESTINATION lib)
