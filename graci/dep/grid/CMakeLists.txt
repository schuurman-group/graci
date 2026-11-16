#
# Minimum version of CMake required
#
cmake_minimum_required(VERSION 4.0)

#
# Project name and language
#
project(grid)
enable_language(Fortran)

#
# Conditional compilation flags
#
option(CBINDING "Optional C binding for public routines" ON)
if (CBINDING)
  add_definitions(-DCBINDING)
endif()

#
# Compiler-specific flags
#
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-g -O3 -fbacktrace -cpp -ffree-form -fimplicit-none \
               -mavx2 -mfma -fopenmp")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-w -g -O3 -cpp -free -implicitnone -traceback \
               -qopenmp -xhost -auto -ipo -qopt-prefetch -fp-model precise \
               -check noarg_temp_created")
  set(bounds "-check bounds")
endif()

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

#
# Add grid library target (molecular quadrature grids: atomic grids,
# Delaunay neighbour lists, partition weights and their nuclear
# derivatives; no external dependencies)
#
file(GLOB_RECURSE GRIDsources src/*)
add_library(grid SHARED ${GRIDsources})

#
# Install
#
install(TARGETS grid DESTINATION lib)
