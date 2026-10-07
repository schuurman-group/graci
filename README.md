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

## MRSF-TDDFT
Mixed-reference spin-flip TDDFT (MRSF-TDDFT, Lee, Filatov, Lee and Choi,
J. Chem. Phys. 149, 104101 (2018)) is available through the `$mrsftddft`
section. It requires an ROKS/ROHF triplet reference (`mult = 3` in the
`$scf` section) and a global hybrid functional (e.g. `xc = bhlyp`); the
density-fitted MO integrals of the `$scf` section are reused. One section
computes one spin manifold:

    $mrsftddft section
     label   = singlets
     mult    = 1               # 1: singlets, 3: triplets
     nstates = [3 1 1 2]       # roots per irrep
    $end

The irreps are those of the MRSF states themselves (the irrep of the
triplet reference, i.e. the product of the two SOMO irreps, times the irrep
of the spin-flip excitation), so the closed-shell ground state S0 belongs
to the totally symmetric irrep.

Optional keywords: `nextra` (extra Davidson roots tracked, 3), `diag_tol`
(residual norm, 1e-5), `diag_iter` (100), `diag_maxvec` (subspace size per
root, 20), `hfx` (override the fraction of HF exchange), `spc` (spin-pair
coupling scale factors, default = HF exchange fraction), `conf_thresh`
(amplitude print threshold, 0.05), `precision` (`double`/`single` integral
storage), `mem_budget` (GB for the sigma-vector work arrays, 1.0),
`vv_storage` (`paired`/`full` storage of the virtual-virtual DF block),
`print_orbitals`, `ref_state`, `scf_label`, `label`. Total energies are
E(ROKS triplet) + omega; the lowest singlet root is S0. State densities,
natural orbitals and moments are produced as for DFT/MRCI, and the
`$transition` section works between two `$mrsftddft` objects of the same
multiplicity (see examples/h2o_mrsf). Range-separated functionals and
spin-orbit coupling are not yet supported.

Analytic nuclear gradients of MRSF-TDDFT states are requested with a
`$mrsfgradient` section referring to an `$mrsftddft` section:

    $mrsfgradient section
     label      = grad
     mrsf_label = singlets     # the $mrsftddft section (optional if unique)
     states     = [1 2]        # adiabatic state indices (default: all)
    $end

The gradients (Hartree/Bohr, in the PySCF frame of the input geometry)
are printed, stored in the checkpoint file (`Mrsfgradient.<label>`:
`grad`, `grad_ref`, `energies`, `state_syms`) and can be written to files
with `gextract <chkpt> -grad <label>` (see examples/ch2o_mrsf_gradient). The implementation solves one
Z-vector equation per state (block preconditioned conjugate gradients on
the ROKS orbital Hessian; `zvec_tol` (1e-8), `zvec_iter` (100),
`zvec_solver = pcg|dense`), evaluates all MO-space two-electron
contractions with the density-fitted integrals of the MRSF library, and
uses PySCF for the AO derivative integrals, the XC terms and the
reference ROKS gradient. The XC quadrature grid is treated as fixed for
the response term (the ROKS reference part can include the grid
response with `grid_response = True`); the resulting error is of the
order of 1e-5 Hartree/Bohr at the default grid level and vanishes with
finer grids. Requirements: `precision = double`, the full MO space
(`mo_cutoff` not truncating), global hybrid or HF functionals.
