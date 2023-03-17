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
option(CBINDING "Optional C binding for public routines" ON)
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
# overlap library directory
#
link_directories(../lib)

#
# Compiler-specific flags
#
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-g -O3 -fbacktrace -cpp -ffree-form -fimplicit-none \
               -mavx -lblas -llapack")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-w -g -00 -check all -cpp -free -implicitnone \
               -traceback -mkl -qopenmp -xhost -auto -ipo \
               -qopt-prefetch -check noarg_temp_created")
  set(bounds "-check bounds")
endif()

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

#
# Add bitCI library target
#
file(GLOB_RECURSE CIsources src/ci/* src/common/* src/integrals/*)
add_library(bitci SHARED ${CIsources})

#
# Add bitSI library target
#
file(GLOB_RECURSE SIsources src/si/* src/common/* src/integrals/*)
add_library(bitsi SHARED ${SIsources})

#
# Add bitWF library target
#
file(GLOB_RECURSE SIsources src/wf/* src/common/* src/integrals/*)
add_library(bitwf SHARED ${SIsources})

# Link the targets to the HDF5 Fortran library
target_link_libraries(bitci -lhdf5_fortran)
target_link_libraries(bitsi -lhdf5_fortran)
target_link_libraries(bitwf -lhdf5_fortran)

# Link the bitci and bitwf target to the overlap library
target_link_libraries(bitci -loverlap)
target_link_libraries(bitwf -loverlap)

#
# Install
#
install(TARGETS bitci bitsi bitwf DESTINATION lib)
