# graci
General Reference Configuration Interaction package

# Python dependencies
GRaCI has a fair few Python dependencies. These may most easily be
handled by using the provided graci.yml Anaconda environment
file. Running

conda env create -f graci.yml

will create an Anaconda environment named 'graci' in which GRaCI may be run

Note, however, that this environment does not include the PySCF
dependency, which must be installed separately

# Build and use
In the following, $TOPDIR will refer to the path to the top graci directory

## Other dependencies
CMake v3.2 or higher

PySCF

## Recomendations
Compile using Intel ifort and MKL for optimal performance

Note that both are now freely available through the Intel OneAPI suite

## Build
(1) cd $TOPDIR/graci/dep

(2) export FC=**fname** (**fname** \in {ifort, gfortran})

(3) ./install_all

## Environment variables
A small number of environment variables need to be
set/appended before graci can be executed.

In bash, this would take the form:

export GRACI=$TOPDIR

export PATH=$PATH:$GRACI/bin

export PYTHONPATH=$GRACI

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GRACI/graci/dep/lib

## Running graci
After setting the above environment variables, simply use the command

graci file.inp

to run a graci calculation, where file.inp is a graci input file
