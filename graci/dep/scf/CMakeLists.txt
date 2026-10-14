#
# Minimum version of CMake required
#
cmake_minimum_required(VERSION 4.0)

#
# Project name and language
#
project(scf)
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
# libxc (the build that the installed PySCF ships; install_scf sets the
# directory from the pyscf package)
#
if(NOT LIBXC_DIR)
  message(FATAL_ERROR "LIBXC_DIR is not set (install_scf sets it from the installed pyscf)")
endif()

#
# Add scf library target
#
file(GLOB_RECURSE SCFsources src/*)
add_library(scf SHARED ${SCFsources})
target_link_libraries(scf ${LIBXC_DIR}/libxc.so)
set_target_properties(scf PROPERTIES INSTALL_RPATH "${LIBXC_DIR}" BUILD_WITH_INSTALL_RPATH TRUE)

#
# Install
#
install(TARGETS scf DESTINATION lib)
