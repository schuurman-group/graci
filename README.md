# graci
General Reference Configuration Interaction package

# Build and use
In the following, $TOPDIR will refer to the path to the top graci directory

## Dependencies
CMake v3.2 or higher

PySCF

HDF5 build **including the Fortran libraries**

## Recomendations
Compile using Intel ifort and MKL for optimal performance

Note that both are now freely available through the Intel OneAPI suite

## Build
(1) cd $TOPDIR/graci/lib

(2) export FC=**fname** (**fname** \in {ifort, gfortran})

(3) ./install_all -DHDF5_INC_DIR=**include_path** -DHDF5_LIB_DIR=**lib_path**

where **include_path** and **lib_path** are the paths to the HDF5 include and lib directories
 
## Environment variables
A small number of environment variables need to be
set/appended before graci can be executed.

In bash, this would take the form:

export GRACI=$TOPDIR

export PATH=$PATH:$GRACI/bin

export PYTHONPATH=$GRACI

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GRACI/graci/lib/lib

## Running graci
After setting the above environment variables, simply use the command

graci file.inp

to run a graci calculation, where file.inp is a graci input file
