#
# Minimum version of CMake required
#
cmake_minimum_required(VERSION 3.2)

#
# Project name and language
#
project(overlap)
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
# Compiler-specific flags
#
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-g -O3 -fbacktrace -cpp -ffree-form -fimplicit-none \
               -mavx -lblas -llapack")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-w -g -O3 -cpp -free -implicitnone -traceback \
               -qmkl -qopenmp -xhost -auto -ipo -qopt-prefetch \
               -check noarg_temp_created")
  set(bounds "-check bounds")
endif()

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

#
# Add overlap library target
#
file(GLOB_RECURSE OVERLAPsources src/*)
add_library(overlap SHARED ${OVERLAPsources})

#
# Install
#
install(TARGETS overlap DESTINATION lib)
