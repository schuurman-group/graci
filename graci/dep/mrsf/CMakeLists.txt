#
# Minimum version of CMake required
#
cmake_minimum_required(VERSION 4.0)

#
# Project name and language
#
project(mrsf)
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
               -mavx2 -mfma -fopenmp -lblas -llapack")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-w -g -O3 -cpp -free -implicitnone -traceback \
               -qmkl -qopenmp -xhost -auto -ipo -qopt-prefetch \
               -check noarg_temp_created")
  set(bounds "-check bounds")
  add_definitions(-DUSE_MKL)
endif()

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

#
# Add mrsf library target
#
file(GLOB_RECURSE MRSFsources src/*)
add_library(mrsf SHARED ${MRSFsources})

#
# Install
#
install(TARGETS mrsf DESTINATION lib)
