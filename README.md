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
`zvec_solver = pcg|dense`) and evaluates everything but the integral
and functional-derivative evaluation inside the MRSF library: the
MO-space two-electron contractions with the density-fitted integrals, the
Z-vector operator, the XC kernel and the XC gradient terms on the
quadrature grid (hole-width factor pairs, no AO matrices) and the
contraction of the derivative integrals with the three-index densities.
PySCF supplies the AO derivative integrals, the grid and the libxc
derivatives. `mem_budget` (GB, 2.0) bounds the cached AO values on the
grid (re-evaluated block by block when they do not fit) and the
auxiliary blocks of the derivative-integral contraction. The XC
quadrature grid is treated as fixed for the response term (the ROKS
reference part can include the grid response with
`grid_response = True`); the resulting error is of the order of 1e-5
Hartree/Bohr at the default grid level and vanishes with finer grids.
Requirements: `precision = double`, the full MO space (`mo_cutoff` not
truncating), global hybrid or HF functionals.

Overlaps between MRSF-TDDFT states of two geometries (for nonadiabatic
dynamics) are computed by the existing `$overlap` section when its
`bra_label` and `ket_label` refer to two `$mrsftddft` objects (two
`$molecule`/`$scf`/`$mrsftddft` sets in one input, same basis and
reference occupation pattern; see examples/h2o_mrsf_overlap):

    $overlap section
     bra_label   = singlets1
     ket_label   = singlets2
     bra_states  = [1:4]
     ket_states  = [1:4]
     mrsf_method = exact         # exact | tlf0 | tlf1 | tlf2
    $end

The MRSF states are expanded in the determinants of the two reference
determinants and their overlaps reduce to two-index determinants of the
MO overlap matrix (determinant factorisation, Lee, Kim, Lee, Choi, JCTC
15, 882 (2019)). `mrsf_method = exact` (default) evaluates these
determinants exactly from one inverse of the hole block and one of the
core block of the MO overlap matrix; `tlf0`, `tlf1`, `tlf2` apply the
truncated Leibniz formula of that paper to the same determinants (orders
0, 1, 2 in the off-diagonal MO overlaps; `tlf0` drops the terms coupling
the two reference determinants). The TLF variants assume that the MOs of
the two geometries are nearly orthonormal in the same order (as for
consecutive MD steps with a consistent MO phase convention); MO
reorderings or mixings between the two calculations, e.g. a
symmetry-blocked and a C1 calculation, break them, whereas `exact` is
invariant. Everything but the AO overlap matrix is evaluated in the MRSF
library; the cost is negligible compared to the two MRSF-TDDFT
calculations (pyrazine/aug-cc-pVTZ, 8 x 8 states: 0.01 s). Overlaps between singlet and triplet
objects vanish identically and are returned as zero; the sign of each
state (and therefore of each row/column) is that of the MRSF
eigenvectors.
