!**********************************************************************
! Parameters of the various DFT/MRCI Hamiltonians
!**********************************************************************
module hparam

  use constants
  
  implicit none

  save

  ! Number of Hamiltonians implemented
  integer(is), parameter :: nham=22

  ! Hamiltonian labels
  character(len=20), parameter, dimension(nham) :: hlbl= &
       ['abinitio            ', &
        'grimme              ', &
        'grimme_short        ', &
        'r2016               ', &
        'r2016_short         ', &
        'r2017               ', &
        'r2017_short         ', &
        'r2018               ', &
        'r2018_short         ', &
        'r2022               ', &
        'qe8                 ', &
        'qe8_short           ', &
        'cvs-qe8             ', &
        'cvs-test            ', &
        'r2026               ', &
        'cvs-r2026           ', &
        'qe8_sym             ', &
        'rc_dftmrci          ', &
        'rc_dftmrci_5p       ', &
        'rc_dftmrci_7p       ', &
        'rc_dftmrci_8p       ', &
        'rc_dftmrci_kj       ']

  ! Hamiltonian integer label
  integer(is)           :: ihamiltonian

  ! Hamiltonian parameter array
  real(dp), allocatable :: hpar(:)

  ! Number of Hamiltonian parameters
  integer(is)           :: nhpar

  ! Configuration selection energy cutoff
  real(dp)              :: desel
  
!----------------------------------------------------------------------
! Grimme's original DFT/MRCI Hamiltonian
! J. Chem. Phys., 111, 5645 (1999)
!----------------------------------------------------------------------
  ! Singlet parameters, delta E_sel = 1.0
  real(dp), parameter, dimension(5) :: grimme1= &
       [0.6195d0, & ! p1
       3.2719d0, &  ! p2
       0.5102d0, &  ! pJ
       0.5945d0, &  ! p[0]
       0.1058d0]    ! alpha
  
  ! Singlet parameters, delta E_sel = 0.8
  real(dp), parameter, dimension(5) :: grimme1_short= &
       [0.6290d0, & ! p1
        8.0000d0, & ! p2
        0.5030d0, & ! pJ
        0.6110d0, & ! p[0]
        0.1190d0]   ! alpha
  
  ! Triplet parameters, delta E_sel = 1.0
  real(dp), parameter, dimension(5) :: grimme3= &
       [0.6195d0, & ! p1
       3.2719d0, &  ! p2
       0.4930d0, &  ! pJ
       0.0000d0, &  ! p[0]
       0.0563d0]    ! alpha

  ! Triplet parameters, delta E_sel = 0.8
  real(dp), parameter, dimension(5) :: grimme3_short= &
       [0.6290d0, & ! p1
       8.0000d0, &  ! p2
       0.4860d0, &  ! pJ
       0.0000d0, &  ! p[0]
       0.0630d0]    ! alpha

!----------------------------------------------------------------------
! Lyskov's 2016 redesigned DFT/MRCI Hamiltonian
! J. Chem. Phys., 144, 034104 (2016)
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(4) :: r2016= &
       [0.507894d0, & ! pJ
       0.355895d0, &  ! pF
       0.568168d0, &  ! p1
       1.788d0]       ! p2
       
  ! delta E_sel = 0.8
  real(dp), parameter, dimension(4) :: r2016_short= &
       [0.503506d0, & ! pJ
       0.368122d0, &  ! pF
       0.579809d0, &  ! p1
       2.187d0]       ! p2

!----------------------------------------------------------------------
! Heil's 2017 DFT/MRCI Hamiltonian for odd and even electron numbers
! J. Chem. Phys., 147, 194104 (2017)
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(4) :: r2017= &
       [0.503001d0, & ! pJ
       0.358727d0, &  ! pF
       0.563893d0, &  ! p1
       1.8571d0]      ! p2

  ! delta E_sel = 0.8
  real(dp), parameter, dimension(4) :: r2017_short= &
       [0.500779d0, & ! pJ
       0.356986d0, &  ! pF
       0.573523d0, &  ! p1
       1.9266d0]      ! p2

!----------------------------------------------------------------------
! Heil's DFT/MRCI Hamiltonian for transition metal complexes
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(4) :: r2018= &
       [0.508918d0, & ! pJ
       0.362362d0, &  ! pF
       0.558411d0, &  ! p1
       4.47165d0]     ! p2

  ! delta E_sel = 0.8
  real(dp), parameter, dimension(4) :: r2018_short= &
       [0.505808d0, & ! pJ
       0.359626d0, &  ! pF
       0.577732d0, &  ! p1
       11.499113d0]   ! p2

!----------------------------------------------------------------------
! R2022 DFT/MRCI Hamiltonian
! J. Phys. Chem A, 127, 2011 (2023)
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(5) :: r2022= &
       [3.4673d0, & ! p2
       0.5085d0,  & ! pJ^he
       0.4649d0,  & ! pJ^hhee
       0.3426d0,  & ! px^he
       0.5416d0]    ! px^hhee
  
!----------------------------------------------------------------------
! QE8 Hamiltonians
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(5) :: qe8= &
       [0.425623d0, & ! pJ
       0.252259d0, &  ! pF
       0.692173d0, &  ! p1
       4.611269d0, &  ! p2
       8.0d0]         ! n

  ! delta E_sel = 0.8
  real(dp), parameter, dimension(5) :: qe8_short= &
       [0.419d0, & ! pJ
       0.258d0,  &  ! pF
       0.712d0,  &  ! p1
       4.69d0,   &  ! p2
       8.0d0]         ! n
  
  real(dp), parameter, dimension(7) :: cvs_qe8= &
       [0.425623d0, & ! pJ_vv
       0.252259d0, &  ! pF_vv
       0.499646d0, &  ! p1
       0.214962d0, &  ! p2
       8.0d0, &       ! n
       0.560644d0, &  ! pJ_cv
       0.252259d0]    ! pF_cv
  
  real(dp), parameter, dimension(6) :: cvs_test= &
       [0.425623d0, & ! pJ_vv
       0.252259d0, &  ! pF
       0.499646d0, &  ! p1
       0.214962d0, &  ! p2
       8.0d0, &       ! n
       0.560644d0]    ! pJ_cv

!----------------------------------------------------------------------
! 2026 DFT/MRCI Hamiltonian
! Heil diagonal form with separate he/hhee exchange scaling;
! off-diagonal prefactor p1 = 1 - 2*pJ + pF_hhee (derived).
! *** Preliminary parameters: these need to be re-fitted ***
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(5) :: r2026= &
       [0.425623d0, & ! pJ
       0.252259d0, &  ! pF_he
       0.252259d0, &  ! pF_hhee
       4.611269d0, &  ! p2
       8.0d0]         ! n

!----------------------------------------------------------------------
! CVS-2026 DFT/MRCI Hamiltonian
! As r2026 but with separate core-valence Coulomb (pJ_cv) and
! exchange (pF_cv) scaling.
! *** Preliminary parameters: these need to be re-fitted ***
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(7) :: cvs_r2026= &
       [0.425623d0, & ! pJ_vv
       0.252259d0, &  ! pF_he_vv
       0.252259d0, &  ! pF_hhee_vv
       4.611269d0, &  ! p2
       8.0d0, &       ! n
       0.425623d0, &  ! pJ_cv
       0.252259d0]    ! pF_cv

!----------------------------------------------------------------------
! RC DFT/MRCI Hamiltonian
! QE8_SYM with range-separated exchange: pF_SR scales the SR part
! K_SR = K_full - K_LR, pF_LR scales the LR part K_LR.
! ω is not a parameter — it comes from the XC functional.
! Degenerates to QE8_SYM when ω=0 (global hybrid, K_LR=0).
! *** Preliminary parameters: initialized to QE8 values ***
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
! Off-diagonal elements carry two independent damping functions:
!
!   H_IJ = damp * H_SR + damp_LR * H_LR
!   damp    = p1    * exp(-p2    * |dE|**n)
!   damp_LR = p1_LR * exp(-p2_LR * |dE|**n_LR)
!
! p1_LR is held at 1 so that the LR contribution is the undamped ab initio
! matrix element as dE -> 0, with p2_LR setting how fast it falls away.
  real(dp), parameter, dimension(10) :: rc_dftmrci_p= &
       [0.425623d0, & ! pJ_SR
       0.425623d0, &  ! pJ_LR
       0.252259d0, &  ! pF_SR
       0.252259d0, &  ! pF_LR
       1.0d0,      &  ! p1_LR (LR off-diagonal pre-factor, held at 1)
       4.611269d0, &  ! p2_LR (LR off-diagonal exponent coefficient)
       8.0d0,      &  ! n_LR  (LR off-diagonal exponent power)
       0.692173d0, &  ! p1    (SR off-diagonal pre-factor)
       4.611269d0, &  ! p2    (SR off-diagonal exponent coefficient)
       8.0d0]         ! n     (SR off-diagonal exponent power)

!----------------------------------------------------------------------
! rc_dftmrci_5p: the reduced, physically-constrained form.
!
! Two changes relative to rc_dftmrci (ihamiltonian=18):
!
!  1. ONE Coulomb scaling, not two. The Hartree term is not range
!     separated in any hybrid -- J[rho] uses the full 1/r12 -- so there
!     is nothing in the functional that would justify pJ_SR /= pJ_LR.
!     Exchange IS range separated, so pF_SR and pF_LR are kept apart.
!
!  2. The long-range off-diagonal contribution carries no damping
!     parameters at all: H_IJ = p1*exp(-p2*dE^n)*H_SR + H_LR. A
!     semilocal correlation functional is local in r12 and has nothing
!     to double count between separated electrons, so the long-range
!     coupling is taken at its ab initio value. The p1_LR/p2_LR/n_LR of
!     ihamiltonian=18 were degenerate with p1 and numerically inert.
!
! Seeded from the CAM-QTP00 fit of the reduced form (31 QUEST
! transitions, RMSD 0.129 eV).
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(6) :: rc_dftmrci_5p_p= &
       [0.420240d0, & ! pJ    (single, full-range Coulomb scaling)
       0.259689d0, &  ! pF_SR
       0.302801d0, &  ! pF_LR
       0.521930d0, &  ! p1    (SR off-diagonal pre-factor)
       2.786858d0, &  ! p2    (SR off-diagonal exponent coefficient)
       8.0d0]         ! n     (SR off-diagonal exponent power)

!----------------------------------------------------------------------
! rc_dftmrci_7p: the physically-motivated form.
!
! Every parameter is answerable to something in the functional:
!
!  pJ        ONE Coulomb scaling. Neither the Hartree operator nor the
!            correlation functional is range separated, so there is no
!            basis for scaling short- and long-range Coulomb differently.
!
!  pK_SR,    Exchange IS range separated, and the long-range channel
!  pK_LR     carries less DFT exchange, so less of it needs correcting.
!
!  p1_SR,    Two prefactors answering two independent questions about
!  p1_LR     the SAME functional: what fraction of each channel's
!            correlation does E_c already hold? H_SR and H_LR partition
!            the operator exactly and disjointly, so this is not double
!            counting. p1_LR is expected near 1 -- but as a prediction to
!            be tested, not a constraint imposed.
!
!  p2, n     SHARED between the channels. dE discriminates static from
!            dynamic correlation, a property of the configuration pair;
!            range discriminates where in r12 the coupling lives. The two
!            axes are orthogonal -- measured directly: mean dE is flat
!            across the per-element LR fraction -- so the dE dependence
!            cannot legitimately differ by channel.
!
! Note the limits: p1_LR -> p1_SR recovers uniform (canonical) damping and
! makes the off-diagonal split vacuous; p1_LR -> 1 is the ab initio
! long-range limit. The fit is free to choose between them.
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(7) :: rc_dftmrci_7p_p= &
       [0.350000d0, & ! pJ     (single, full-range Coulomb scaling)
       0.260000d0, &  ! pK_SR
       0.350000d0, &  ! pK_LR
       0.520000d0, &  ! p1_SR  (SR off-diagonal pre-factor)
       0.850000d0, &  ! p1_LR  (LR off-diagonal pre-factor)
       2.790000d0, &  ! p2     (shared exponent coefficient)
       8.0d0]         ! n      (shared exponent power)

!----------------------------------------------------------------------
! rc_dftmrci_8p: split Coulomb, split exchange, shared damping shape.
!
! Differs from rc_dftmrci_7p by splitting the Coulomb scaling again. That
! split is NOT justified by range separation of the functional -- neither
! J[rho] nor E_c[rho] knows about omega -- but by what the split actually
! does. The LR fraction of the Coulomb integral (hh|aa),
!
!     f_J = (hh|aa)_LR / (hh|aa)
!
! measures how diffuse the target orbital is: 0.67 for a compact valence
! pi* in formaldehyde, 0.93 for a Rydberg orbital. So
!
!     pJ_eff = pJ_SR (1 - f_J) + pJ_LR f_J
!
! is a state-dependent Coulomb correction -- 0.23 for valence, 0.16 for
! Rydberg with the fitted values -- applied with no state classification
! anywhere in the code. The DFT correlation error genuinely does differ
! between compact and diffuse excited states, so this is physical; it is
! simply not the physics the SR/LR label advertises. Constraining
! pJ_SR = pJ_LR removes the state dependence entirely and costs a factor
! of 1.9 in RMSD.
!
! Exchange keeps its split too, though the same analysis shows it does
! very little: K is 97% short range for valence and still 83% for
! Rydberg, so pK_LR multiplies almost nothing and is weakly determined.
! It is retained because that may not hold for other state types.
!
! The damping SHAPE (p2, n) is shared between the channels: dE separates
! static from dynamic correlation, a property of the configuration pair,
! while range separates where in r12 the coupling sits. The two are
! measurably independent -- mean dE is flat across the per-element LR
! fraction -- so the shape cannot legitimately differ by channel. Only
! the prefactors do.
!
! p1_LR went to 1.000 in all seven functionals tested (omega 0.29-0.42,
! SR exact exchange 0.16-0.54, with and without VV10 nonlocal
! correlation), so it may reasonably be frozen at 1.
!
! Seeded from the CAM-QTP00 six-parameter fit (RMSD 0.129 eV).
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(8) :: rc_dftmrci_8p_p= &
       [0.420240d0, & ! pJ_SR
       0.134828d0, &  ! pJ_LR
       0.259689d0, &  ! pK_SR
       0.302801d0, &  ! pK_LR
       0.521930d0, &  ! p1_SR  (SR off-diagonal pre-factor)
       1.000000d0, &  ! p1_LR  (LR off-diagonal pre-factor)
       2.786858d0, &  ! p2     (shared exponent coefficient)
       8.0d0]         ! n      (shared exponent power)

!----------------------------------------------------------------------
! rc_dftmrci_kj: rc_dftmrci_8p with a K/J-dependent LR Coulomb scaling.
!
! Measured 2026-09-09 on the NTOs of the fit-set states (aug-cc-pVTZ,
! CAM-QTP00, omega=0.29): f_J does NOT separate a 90-degree twisted CT
! state from a Rydberg state --
!
!   twisted DMABN CT   f_J = 0.835, 0.896
!   Rydberg (CO,H2CO)  f_J = 0.836, 0.841, 0.874
!
! -- the CT values straddle the Rydberg range, so no threshold in f_J can
! tell them apart.  The pJ_LR scans nevertheless show the two want
! opposite pJ_LR (DMABN 0.000, Rydberg 0.050, nitro sCT 0.100).
!
! The ratio K/J = symvx(i,j)/Vc(i,j) DOES separate them, with a clean gap
! and in the order the scans require:
!
!   nitro sCT   K/J = 0.134, 0.150
!   Rydberg     K/J = 0.018, 0.023, 0.080
!   DMABN CT    K/J = 0.0031, 0.0094
!
! because exchange requires hole/particle overlap and Coulomb does not: a
! displaced CT pair kills K while J survives, whereas a concentric diffuse
! Rydberg pair keeps K small but finite.  J itself does no work here --
! it is 0.185-0.232 for every state measured -- so K/J is really K made
! dimensionless.
!
! The LR Coulomb scaling therefore becomes pair-dependent:
!
!   pJ_LR_eff(i,j) = pJ_LR + pJ_K * K(i,j)/sqrt(J(i,i) J(j,j))
!
! The normalisation is Cauchy-Schwarz, not K/J: the ratio is then bounded
! in [0,1] and reads as the fraction of charge density lying where hole and
! particle coincide.  See kj_ratio in dftmrci.f90.
!
! One extra parameter.  pJ_K = 0 reduces this EXACTLY to rc_dftmrci_8p,
! so the 8p results are recoverable and the seed is the 8p vector with a
! zero appended.
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(9) :: rc_dftmrci_kj_p= &
       [0.601449d0, & ! pJ_SR
       0.034748d0, &  ! pJ_LR
       0.347687d0, &  ! pK_SR
       0.000293d0, &  ! pK_LR
       0.562784d0, &  ! p1_SR
       0.827828d0, &  ! p1_LR
       5.966700d0, &  ! p2
       8.0d0, &       ! n
       0.0d0]         ! pJ_K   (0 => identical to rc_dftmrci_8p)





contains

!######################################################################
! load_hpar: loads the DFT/MRCI Hamiltonian parameters
!######################################################################
  subroutine load_hpar(iham)

    use constants
    use bitglobal
    use iomod
    
    implicit none

    ! Hamiltonian integer label
    integer(is), intent(in) :: iham

!----------------------------------------------------------------------
! Set the Hamiltonian integer label
!----------------------------------------------------------------------
    ihamiltonian=iham

!----------------------------------------------------------------------
! Load the Hamiltonian parameters
!----------------------------------------------------------------------
    select case(iham)

    case(1)
       ! Canonical: do nothing
       ldftmrci=.false.
       desel=999.9d0
       return
       
    case(2)
       ! Grimme, standard
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       if (imult == 1) then
          hpar=grimme1
          desel=1.0d0
       else if (imult == 3) then
          hpar=grimme3
          desel=1.0d0
       else
          errmsg='Only singlet and triplet states are supported using' &
               //' Grimme''s Hamiltonian'
          call error_control
       endif
       
    case(3)
       ! Grimme, short
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       if (imult == 1) then
          hpar=grimme1_short
          desel=0.8d0
       else if (imult == 3) then
          hpar=grimme3_short
          desel=0.8d0
       else
          errmsg='Only singlet and triplet states are supported using' &
               //' Grimme''s Hamiltonian'
          call error_control
       endif
       
    case(4)
       ! Lyskov, standard
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2016
       desel=1.0d0
       
    case(5)
       ! Lyskov, short
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2016_short
       desel=0.8d0
       
    case(6)
       ! Heil17, standard
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2017
       desel=1.0d0
       
    case(7)
       ! Heil17, short
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2017_short
       desel=0.8d0
       
    case(8)
       ! Heil18, standard
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2018
       desel=1.0d0
              
    case(9)
       ! Heil18, short
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2018_short
       desel=0.8d0

    case(10)
       ! R2022
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       hpar=r2022
       desel=1.0d0
       
    case(11)
       ! QE8
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       hpar=qe8
       desel=1.0d0

    case(12)
       ! QE8, short
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       hpar=qe8_short
       desel=0.8d0
       
    case(13)
       ! CVS-QE8
       ldftmrci=.true.
       nhpar=7
       allocate(hpar(nhpar))
       hpar=cvs_qe8
       desel=1.0d0

    case(14)
       ! CVS-TEST
       ldftmrci=.true.
       nhpar=6
       allocate(hpar(nhpar))
       hpar=cvs_test
       desel=1.0d0

    case(15)
       ! R2026
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       hpar=r2026
       desel=1.0d0

    case(16)
       ! CVS-R2026
       ldftmrci=.true.
       nhpar=7
       allocate(hpar(nhpar))
       hpar=cvs_r2026
       desel=1.0d0

    case(17)
       ! QE8 with symmetry-averaged exchange (qe8_sym)
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       hpar=qe8
       desel=1.0d0

    case(18)
       ! RC DFT/MRCI: range-corrected exchange and Coulomb with separate SR/LR scaling
       ldftmrci=.true.
       nhpar=10
       allocate(hpar(nhpar))
       hpar=rc_dftmrci_p
       desel=1.0d0

    case(19)
       ! RC DFT/MRCI, reduced form: one Coulomb scaling, split exchange,
       ! undamped long-range off-diagonal contribution
       ldftmrci=.true.
       nhpar=6
       allocate(hpar(nhpar))
       hpar=rc_dftmrci_5p_p
       desel=1.0d0

    case(20)
       ! RC DFT/MRCI, physically-motivated form: one Coulomb scaling,
       ! split exchange, two off-diagonal prefactors sharing one damping
       ! function
       ldftmrci=.true.
       nhpar=7
       allocate(hpar(nhpar))
       hpar=rc_dftmrci_7p_p
       desel=1.0d0

    case(21)
       ! RC DFT/MRCI: split Coulomb and exchange, two off-diagonal
       ! prefactors sharing one damping function
       ldftmrci=.true.
       nhpar=8
       allocate(hpar(nhpar))
       hpar=rc_dftmrci_8p_p
       desel=1.0d0

    case(22)
       ! RC DFT/MRCI with K/J-dependent LR Coulomb scaling
       ldftmrci=.true.
       nhpar=9
       allocate(hpar(nhpar))
       hpar=rc_dftmrci_kj_p
       desel=1.0d0

    case default
       ! Unrecognised Hamiltonian
       write(errmsg,'(a,x,i0)') &
            'Error in load_hpar: unrecognised Hamiltonian number',&
            iham
       call error_control
       
    end select
       
    return
    
  end subroutine load_hpar

!######################################################################
  
end module hparam
