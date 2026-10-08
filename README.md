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

Extended MRSF-TDDFT (EMRSF-TDDFT, Oh, Kim, Jung, Choi and Lee, ChemRxiv
2026, doi 10.26434/chemrxiv.15000818) adds the core-to-virtual (CV)
configurations of the closed-shell configuration G = |C O1^2| (the
excitations present in linear-response TDDFT but absent from MRSF-TDDFT),
which stabilise charge-transfer states. It is selected with
`extended = True` in the `$mrsftddft` section; `ccp` scales the coupling
block between the MRSF and CV configurations (default: the fraction of
HF exchange, as in the paper). The CV block uses the closed-shell KS
matrix of G evaluated with the triplet orbitals, corrected on its
diagonal by (1 - c_HF)[(pp|O2O2) - (pp|O1O1)], the shift A_G of the paper
and, for singlets, the XC kernel of the functional at the G density; the
coupling block contains the exact Hamiltonian matrix elements between the
MRSF and CV configuration state functions (the paper's eq. 8), scaled by
`ccp`. Energies, amplitudes, the CV weight gamma_CV, state densities
(natural orbitals, moments) and transition densities (oscillator
strengths through `$transition` between two extended objects of the same
multiplicity) are available. The densities are the expectation values
<Psi_I|E_pq|Psi_J> of the expansion in MRSF and CV configuration state
functions, as for MRSF-TDDFT. Overlaps and gradients are not yet
available for extended objects. The kernel needs the AO
values of the DFT grid cached in memory: `mem_budget` must cover
8*ncomp*ngrid*(nao+nocc) bytes (about 1 GB for pyrazine/aug-cc-pVTZ at
grid level 2; the run stops with the required value otherwise). The cost
per sigma vector is about twice that of MRSF-TDDFT for HF or triplets and
about four to five times with the kernel (see examples/h2o_mrsf_extended).

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
     mrsf_align  = True          # align the ket MOs to the bra MOs
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
the two geometries are nearly orthonormal in the same order; with
`mrsf_align = True` (default) the ket MOs are first aligned to the bra
MOs within each orbital class (core, SOMO pair, virtuals) by orthogonal
Procrustes rotations obtained from the MO overlap blocks, which removes
reorderings, sign flips and mixings of (near-)degenerate orbitals
between the two calculations, e.g. between a symmetry-blocked and a C1
calculation. The MRSF states are invariant under these rotations, so
the ket amplitudes are transformed exactly and the `exact` overlaps do
not change; the alignment only makes the TLF expansions valid.
Everything but the AO overlap matrix is evaluated in the MRSF library;
the cost is negligible compared to the two MRSF-TDDFT calculations
(pyrazine/aug-cc-pVTZ, 8 x 8 states: 0.01 s). Overlaps between singlet and triplet
objects vanish identically and are returned as zero; the sign of each
state (and therefore of each row/column) is that of the MRSF
eigenvectors.
