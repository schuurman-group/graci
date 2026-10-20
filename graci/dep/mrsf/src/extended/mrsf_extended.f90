!**********************************************************************
! mrsf_extended: extended MRSF-TDDFT (EMRSF-TDDFT; Oh, Kim, Jung, Choi,
! Lee, ChemRxiv 2026) and its O1-O2 rotation-invariant form EMRSF-Pi
! (plan Part XVIII). The response vector gains the core-to-virtual (CV)
! amplitudes y(a,i), a in P (SOMO rows masked), i in C, of a closed-shell-
! like parent of the SOMO pair: G = |C Cbar O1 O1bar| (published model)
! or the invariant lowest pair state Pi = wG G + wD D + wS S, with
! D = |C Cbar O2 O2bar| and S the open-shell singlet pair (the MRSF OO
! slot); stored as the extra columns nocca+1..ncol of the vector [x | y].
!
! Sigma contributions handled here (sigma_core does the one-body terms
! and the exchange sweep, which also delivers the exchange-type coupling
! terms through Scv = sum_Q T_cv B^Q_{HH}, T_cv = sum_b B^Q_ab y_bj):
!   CV block:  2 s_J (jb|ld) y (singlet; s_J = 1, or kappa for the
!              kernel-free model), + shift y, + K^xc y (singlet, DFT);
!              -c_H (jl|bd) y and the one-body terms come from sigma_core
!   couplings: C y into the MRSF slots and C^T x into the CV columns,
!              scaled by c_cp. Closed forms per parent X in {G, D, S},
!              pinned against the CI in ~/calculations/mrsf_dev/inv/
!              coef_inv.py (the G forms are those of emrsf_ref.py), and
!              combined with the parent weights w_X.  MRSF slot classes:
!              CO1 (row 1 of the core columns), CO2 (row 2), O1V/O2V
!              (columns hO1/hO2), CV, OO (1,hO1), G slot (1,hO2), D slot
!              (2,hO1). Term types: Coulomb-type sum_Q B^Q_{p O_x} j_Y,
!              exchange-type elements of Scv, Fock-type elements of the
!              parent's one-body operator F^XX, SOMO-pair integrals
!              (W-type vectors, the Mp/Gp blocks).  Example (singlet CV
!              CSFs of G, slot (j in C -> b in V)):
!     CO1 (p):  -2 (jb|O2p) + (jp|bO2) - d_pj F_{O2 b}
!     CO2 (p):  -d_pj (O2O1|O2b)
!     O1V (q):  +d_qb (jO1|O2O1)
!     O2V (q):  +2 (jb|qO1) - (jO1|qb) - d_qb F_{j O1}
!     CV (p,q): -d_pj (qO1|O2b) + d_qb (jO1|O2p)
!     OO:       +(jO1|O2b) - 2 (jb|O2O1)
!     G:        +sqrt2 F_{jb};  D: 0
! The published model is the special case w = (1,0,0) with F^G (eq. 17)
! as the one-body operator and F' (eq. 16, diagonal or covariant) in the
! CV block; EMRSF-Pi takes the invariant operators built by the Python
! driver (graci/interfaces/mrsf/mrsf_pi.py).
!**********************************************************************
module mrsf_extended

  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_space
  use mrsf_xcgrid, only: xc_ready, ao_cached, ref_ready, cs_ready, xc_closed_shell_init, xc_kernel_cv

  implicit none

  real(dp), allocatable :: jbjb(:,:)        ! (nvirb,nC) (jb|jb)
  real(dp), allocatable :: NpT(:,:,:,:)     ! (nV,nC,2,2) NpT(b,j,x,a) = (j O_x|O_a b)
  real(dp), allocatable :: Fvec(:,:,:)      ! (nV,2,3)  F^XX(O_x, b), parent X
  real(dp), allocatable :: Fhol(:,:,:)      ! (nC,2,3)  F^XX(j, O_x)
  real(dp), allocatable :: Fblk3(:,:,:)     ! (nV,nC,3) F^XX(b, j)
  real(dp), allocatable :: Wv(:,:,:)        ! (nV,3,2)  (O_x O_y|O_z b), xy = 1: 11, 2: 22, 3: 12
  real(dp), allocatable :: Wc(:,:,:)        ! (nC,3,2)  (O_x O_y|O_z j)
  ! combined (w-weighted) coupling data of the current multiplicity
  real(dp)              :: cJco(2,2) = 0.0_dp, cKco(2,2) = 0.0_dp   ! (slot row s, SOMO x)
  real(dp)              :: cJov(2,2) = 0.0_dp, cKov(2,2) = 0.0_dp   ! (slot column y, SOMO x)
  real(dp)              :: cJoo(3) = 0.0_dp, cKoo(3,4) = 0.0_dp     ! slots OO, G, D; Scv elements
  real(dp), allocatable :: fco(:,:)         ! (nV,2) Fock + W vectors of the CO slots
  real(dp), allocatable :: fov(:,:)         ! (nC,2) Fock + W vectors of the OV slots
  real(dp), allocatable :: Mcomb(:,:)       ! (nV,nV) CV A-type block
  real(dp), allocatable :: Gcomb(:,:)       ! (nC,nC) CV B-type block
  real(dp), allocatable :: Fblk(:,:,:)      ! (nV,nC,3) full-block Fock terms of the slots OO, G, D
  integer(is)           :: coef_mult = 0
  ! kernel-free operator A0 for the inner correction equations of the
  ! Davidson solver (set around the inner sigma batches)
  logical               :: kernel_off = .false.
  ! F' correction of the published model: 1 = diagonal, 2 = covariant VV and
  ! CC blocks of (1 - c_H)(J[D_O2] - J[D_O1]) (Part XII); 0 = none
  integer(is)           :: fprime_mode = 1
  ! parent: 1 = G (published model), 2 = Pi (invariant lowest pair state)
  integer(is)           :: ext_mode = 1
  real(dp)              :: wpi(3) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
  real(dp)              :: coul_scale = 1.0_dp    ! scale of the singlet CV Coulomb term

contains

!######################################################################
! ext_initialise: one-body operators, coupling data, shift, diagonal
!   fone(nmo,nmo): one-body operator of the CV block (published: F^DFT of
!                  G, corrected to F' here; Pi: the invariant operator)
!   fcpl(nmo,nmo,3): one-body operators of the Fock-type coupling terms
!                  per parent G, D, S (published: F^DFT in column 1)
!   wpi1(3): parent weights; shift1: CV shift (Pi; the published A_G is
!   computed here); coul1: singlet Coulomb scale; mode1: 1 published, 2 Pi
!######################################################################
  subroutine ext_initialise(fone, fcpl, wpi1, shift1, ccp1, use_kernel1, coul1, mode1, fprime1)

    real(dp), intent(in)    :: fone(nmo,nmo), fcpl(nmo,nmo,3), wpi1(3), shift1, ccp1, coul1
    logical, intent(in)     :: use_kernel1
    integer(is), intent(in) :: mode1, fprime1

    real(dp), allocatable :: corr(:), dO1(:), dO2(:), Dpart(:,:), dQ(:), Jvv(:,:), Jcc(:,:), Bq(:,:)
    integer(is)           :: p, a, j, b, blk, Ql, Q, x, y, z, xy, iOx(2)
    real(dp)              :: o2o2o1o1

    if (.not. extended) call mrsf_error('ext_initialise: extended flag not set')
    if (.not. ints_loaded) call mrsf_error('ext_initialise: integrals not loaded')
    if (.not. allocated(Bcv)) call mrsf_error('ext_initialise: particle-core block missing')
    if (mode1 /= 1 .and. mode1 /= 2) call mrsf_error('ext_initialise: mode must be 1 or 2')
    call ext_free()

    ccp = ccp1
    use_kernel = use_kernel1
    ext_mode = mode1
    wpi = wpi1
    coul_scale = coul1
    if (use_kernel) then
       if (.not. xc_ready) call mrsf_error('ext_initialise: the grid kernel is not initialised')
       if (.not. ao_cached) call mrsf_error('ext_initialise: the extended method needs the '// &
            'AO values cached on the grid (raise mem_budget)')
       ! the closed-shell kernel of the CV block: installed by the driver
       ! (xc_cs_from_g) or built here from a reference-type kernel
       if (.not. cs_ready) then
          if (.not. ref_ready) call mrsf_error('ext_initialise: closed-shell kernel not installed')
          call xc_closed_shell_init()
       endif
    endif

    allocate(corr(nmo), dO1(naux), dO2(naux))
    dO1 = Dall(:,iO1)
    dO2 = Dall(:,iO2)
    o2o2o1o1 = dot_product(dO2, dO1)
    corr = 0.0_dp
    fprime_mode = 0
    if (ext_mode == 1) then
       if (fprime1 /= 1 .and. fprime1 /= 2) call mrsf_error('ext_initialise: fprime must be 1 or 2')
       fprime_mode = fprime1
       ! diagonal correction (1 - c_H)[(pp|O2O2) - (pp|O1O1)]
       call dgemv('T', naux, nmo, 1.0_dp, Dall, naux, dO2, 1, 0.0_dp, corr, 1)
       call dgemv('T', naux, nmo, -1.0_dp, Dall, naux, dO1, 1, 1.0_dp, corr, 1)
       corr = (1.0_dp - chf) * corr
       if (fprime_mode == 2) then
          if (store_sp) call mrsf_error('ext_initialise: the covariant F'' needs double-precision integrals')
          allocate(dQ(naux), Jvv(nvirb,nvirb), Jcc(max(nC,1_is),max(nC,1_is)))
          dQ = (1.0_dp - chf) * (dO2 - dO1)
          call ext_coulomb_blocks(dQ, Jvv, Jcc)
       endif
    endif

    ! one-body operator of the CV block (V x V embedded in the particle
    ! layout, C x C), with the F' correction of the published model
    allocate(Fp_emb(nvirb,nvirb), source=0.0_dp)
    allocate(Fp_cc(max(nC,1_is),max(nC,1_is)), source=0.0_dp)
    do b = 3, nvirb
       do a = 3, nvirb
          Fp_emb(a,b) = fone(Pmap(a),Pmap(b))
          if (fprime_mode == 2) Fp_emb(a,b) = Fp_emb(a,b) + Jvv(a,b)
       enddo
       if (fprime_mode == 1) Fp_emb(b,b) = Fp_emb(b,b) + corr(Pmap(b))
    enddo
    do j = 1, nC
       do p = 1, nC
          Fp_cc(p,j) = fone(Hmap(p),Hmap(j))
          if (fprime_mode == 2) Fp_cc(p,j) = Fp_cc(p,j) + Jcc(p,j)
       enddo
       if (fprime_mode == 1) Fp_cc(j,j) = Fp_cc(j,j) + corr(Hmap(j))
    enddo
    if (allocated(dQ)) deallocate(dQ, Jvv, Jcc)

    ! Fock-type elements of the parents' one-body operators
    iOx = (/ iO1, iO2 /)
    allocate(Fvec(max(nV,1_is),2,3), Fhol(max(nC,1_is),2,3), Fblk3(max(nV,1_is),max(nC,1_is),3), source=0.0_dp)
    do z = 1, 3
       do x = 1, 2
          do b = 1, nV
             Fvec(b,x,z) = fcpl(iOx(x),Pmap(2+b),z)
          enddo
          do j = 1, nC
             Fhol(j,x,z) = fcpl(Hmap(j),iOx(x),z)
          enddo
       enddo
       do j = 1, nC
          do b = 1, nV
             Fblk3(b,j,z) = fcpl(Pmap(2+b),Hmap(j),z)
          enddo
       enddo
    enddo

    ! B^Q_{O1O2} and the SOMO-pair integral vectors (O_x O_y|O_z b), (O_x O_y|O_z j)
    allocate(Bo12(naux), Bq(naux,3))
    do Q = 1, naux
       blk = (Q-1)/nQ + 1
       Ql  = Q - (blk-1)*nQ
       Bo12(Q) = Boo(Ql,nC+1,nC+2,blk)
    enddo
    Bq(:,1) = dO1
    Bq(:,2) = dO2
    Bq(:,3) = Bo12
    allocate(Wv(max(nV,1_is),3,2), Wc(max(nC,1_is),3,2), source=0.0_dp)
    do z = 1, 2
       do xy = 1, 3
          if (nV > 0) call dgemv('T', naux, nV, 1.0_dp, Bvo(1,1,z), naux, Bq(1,xy), 1, 0.0_dp, Wv(1,xy,z), 1)
          if (nC > 0) call dgemv('T', naux, nC, 1.0_dp, Bco(1,1,z), naux, Bq(1,xy), 1, 0.0_dp, Wc(1,xy,z), 1)
       enddo
    enddo

    ! transposed SOMO-pair blocks (contiguous columns for the rank-1 updates)
    allocate(NpT(max(nV,1_is),max(nC,1_is),2,2))
    do y = 1, 2
       do x = 1, 2
          do j = 1, nC
             do b = 1, nV
                NpT(b,j,x,y) = Np(j,b,x,y)
             enddo
          enddo
       enddo
    enddo

    ! shift of the CV block: A_G (eq. 13) for the published model, the
    ! lowest pair-block eigenvalue for Pi
    if (ext_mode == 1) then
       A_G = FbPP(1,1) - FaHH(nC+2,nC+2) - chf * o2o2o1o1
    else
       A_G = shift1
    endif

    ! CV diagonal without the shift: F_bb - F_jj - c_H (jj|bb) [+ 2 s_J (jb|jb) singlet]
    allocate(diag_cv(nvirb,max(nC,1_is)), source=1.0e20_dp)
    allocate(jbjb(nvirb,max(nC,1_is)), source=0.0_dp)
    allocate(Dpart(naux,nvirb))
    do a = 1, nvirb
       Dpart(:,a) = Dall(:,Pmap(a))
    enddo
    do j = 1, nC
       do a = 3, nvirb
          diag_cv(a,j) = Fp_emb(a,a) - Fp_cc(j,j) - chf * dot_product(Dpart(:,a), Dall(:,Hmap(j)))
          jbjb(a,j) = sum(Bcv(a,j,:)**2)
       enddo
    enddo
    deallocate(corr, dO1, dO2, Dpart, Bq)

    allocate(fco(max(nV,1_is),2), fov(max(nC,1_is),2), Mcomb(max(nV,1_is),max(nV,1_is)), &
         Gcomb(max(nC,1_is),max(nC,1_is)), Fblk(max(nV,1_is),max(nC,1_is),3))
    coef_mult = 0
    ext_ready = .true.
    if (verbose) then
       write(6,'(/,2x,a)') 'Extended MRSF-TDDFT initialised'
       write(6,'(2x,a,i0,a,i0,a,f8.4,a,f12.6,a,l1)') 'CV columns = ', nC, ', xdim_tot = ', &
            xdim_tot, ', c_cp = ', ccp, ', shift = ', A_G, ' Ha, kernel = ', use_kernel
       if (ext_mode == 2) then
          write(6,'(2x,a,3f9.5,a,f6.3)') 'parent: invariant lowest pair state, w(G, D, S) = ', wpi, &
               ', Coulomb scale ', coul_scale
       else if (fprime_mode == 2) then
          write(6,'(2x,a)') "parent G; F' correction: covariant (full VV and CC blocks)"
       else
          write(6,'(2x,a)') "parent G; F' correction: diagonal (published form)"
       endif
    endif

  end subroutine ext_initialise

!######################################################################
! ext_coulomb_blocks: Jvv(a,b) = sum_Q B^Q_ab jq(Q) (a, b in P) and
! Jcc(p,j) = sum_Q B^Q_pj jq(Q) (p, j in C) from the resident blocks
! (paired or full planes, Q-blocked Boo); double precision only
!######################################################################
  subroutine ext_coulomb_blocks(jq, Jvv, Jcc)

    real(dp), intent(in)  :: jq(naux)
    real(dp), intent(out) :: Jvv(nvirb,nvirb)
    real(dp), intent(out) :: Jcc(:,:)
    real(dp), allocatable :: Jlow(:,:), Jup(:,:)
    integer(is)           :: k, Q, a, b, blk, Ql, p, j

    Jvv = 0.0_dp
    if (vv_full) then
       do Q = 1, naux
          Jvv = Jvv + jq(Q) * Bvv(:,:,Q)
       enddo
    else
       ! plane k: B^{2k-1} in the lower triangle (incl. the diagonal),
       ! B^{2k} in the strict upper triangle
       allocate(Jlow(nvirb,nvirb), Jup(nvirb,nvirb), source=0.0_dp)
       do k = 1, nplane
          Jlow = Jlow + jq(2*k-1) * Bvv(:,:,k)
          if (2*k <= naux) Jup = Jup + jq(2*k) * Bvv(:,:,k)
       enddo
       do b = 1, nvirb
          do a = b+1, nvirb
             Jvv(a,b) = Jlow(a,b) + Jup(b,a)
             Jvv(b,a) = Jvv(a,b)
          enddo
       enddo
       deallocate(Jlow, Jup)
    endif
    ! diagonal from the stored diagonals (both parities)
    do a = 1, nvirb
       Jvv(a,a) = dot_product(Dall(:,Pmap(a)), jq)
    enddo
    Jcc = 0.0_dp
    do Q = 1, naux
       blk = (Q-1)/nQ + 1
       Ql  = Q - (blk-1)*nQ
       do j = 1, nC
          do p = 1, nC
             Jcc(p,j) = Jcc(p,j) + jq(Q) * Boo(Ql,p,j,blk)
          enddo
       enddo
    enddo

  end subroutine ext_coulomb_blocks

!######################################################################
! parent_tables: pinned closed-form coefficients of the couplings
! between the MRSF slots and the CV CSFs of parent X (1 G, 2 D, 3 S) for
! the multiplicity mult (inv/coef_inv.py). Index conventions:
!   tJco(s,x): Coulomb-type  sum_Q B^Q_{p O_x} j_Y  into the CO slot row s
!   tKco(s,x): exchange-type Scv(x,p)               into the CO slot row s
!   tFco(s,x): Fock-type d_pj F^XX(O_x b)           into the CO slot row s
!   tWco(s,xy,z): d_pj (O_x O_y|O_z b)              into the CO slot row s
!   tJov/tKov/tFov(y,x), tWov(y,xy,z): the same for the OV slot column hO_y
!   tA(x,y): d_pj (q O_x|b O_y) Mp(q,b,x,y);  tB(x,y): d_qb (j O_x|p O_y) Gp(j,p,x,y)
!   tJoo(slot): Bo12-type Coulomb into the slots OO, G, D
!   tKoo(slot,k): Scv(2,hO1), Scv(1,hO2), Scv(1,hO1), Scv(2,hO2) into the slots
!   tFblk(slot): sum_jb F^XX(b,j) y_bj into the slots
!######################################################################
  subroutine parent_tables(mult, X, tJco, tKco, tFco, tWco, tJov, tKov, tFov, tWov, tA, tB, tJoo, tKoo, tFblk)

    integer(is), intent(in) :: mult, X
    real(dp), intent(out)   :: tJco(2,2), tKco(2,2), tFco(2,2), tWco(2,3,2)
    real(dp), intent(out)   :: tJov(2,2), tKov(2,2), tFov(2,2), tWov(2,3,2)
    real(dp), intent(out)   :: tA(2,2), tB(2,2), tJoo(3), tKoo(3,4), tFblk(3)

    real(dp), parameter :: r2 = 1.4142135623730951_dp, h2 = 0.7071067811865476_dp
    real(dp), parameter :: q2 = 0.35355339059327373_dp, t2 = 1.0606601717798212_dp
    logical :: sgl

    sgl = (mult == 1)
    tJco = 0.0_dp; tKco = 0.0_dp; tFco = 0.0_dp; tWco = 0.0_dp
    tJov = 0.0_dp; tKov = 0.0_dp; tFov = 0.0_dp; tWov = 0.0_dp
    tA = 0.0_dp; tB = 0.0_dp; tJoo = 0.0_dp; tKoo = 0.0_dp; tFblk = 0.0_dp
    select case (X)
    case (1)   ! parent G
       if (sgl) tJco(1,2) = -2.0_dp
       tKco(1,2) = 1.0_dp
       tFco(1,2) = -1.0_dp
       tWco(2,3,2) = -1.0_dp
       tWov(1,3,1) = merge(1.0_dp, -1.0_dp, sgl)
       if (sgl) tJov(2,1) = 2.0_dp
       tKov(2,1) = merge(-1.0_dp, 1.0_dp, sgl)
       tFov(2,1) = merge(-1.0_dp, 1.0_dp, sgl)
       tA(1,2) = -1.0_dp
       tB(1,2) = merge(1.0_dp, -1.0_dp, sgl)
       if (sgl) tJoo(1) = -2.0_dp
       tKoo(1,1) = 1.0_dp
       if (sgl) tFblk(2) = r2
    case (2)   ! parent D (the O1 <-> O2 mirror of G)
       tWco(1,3,1) = -1.0_dp
       if (sgl) tJco(2,1) = -2.0_dp
       tKco(2,1) = 1.0_dp
       tFco(2,1) = -1.0_dp
       if (sgl) tJov(1,2) = 2.0_dp
       tKov(1,2) = merge(-1.0_dp, 1.0_dp, sgl)
       tFov(1,2) = merge(-1.0_dp, 1.0_dp, sgl)
       tWov(2,3,2) = merge(1.0_dp, -1.0_dp, sgl)
       tA(2,1) = -1.0_dp
       tB(2,1) = merge(1.0_dp, -1.0_dp, sgl)
       if (sgl) tJoo(1) = 2.0_dp
       tKoo(1,2) = merge(-1.0_dp, 1.0_dp, sgl)
       if (sgl) tFblk(3) = r2
    case (3)   ! parent S (open-shell singlet pair)
       if (sgl) tJco(1,1) = -r2
       tKco(1,1) = h2
       tWco(1,1,1) = -q2
       tWco(1,3,2) = t2
       tFco(1,1) = -h2
       if (sgl) tJco(2,2) = r2
       tKco(2,2) = -h2
       tWco(2,3,1) = -t2
       tWco(2,2,2) = q2
       tFco(2,2) = h2
       if (sgl) tJov(1,1) = r2
       tKov(1,1) = merge(-h2, h2, sgl)
       tWov(1,1,1) = merge(q2, -q2, sgl)
       tWov(1,3,2) = merge(-t2, t2, sgl)
       tFov(1,1) = merge(-h2, h2, sgl)
       if (sgl) tJov(2,2) = -r2
       tKov(2,2) = merge(h2, -h2, sgl)
       tWov(2,3,1) = merge(t2, -t2, sgl)
       tWov(2,2,2) = merge(-q2, q2, sgl)
       tFov(2,2) = merge(h2, -h2, sgl)
       tA(1,1) = -h2
       tA(2,2) = h2
       tB(1,1) = merge(h2, -h2, sgl)
       tB(2,2) = merge(-h2, h2, sgl)
       if (sgl) then
          tFblk(1) = r2
          tJoo(2) = -2.0_dp
          tKoo(2,2) = 1.0_dp
          tJoo(3) = 2.0_dp
          tKoo(3,1) = -1.0_dp
       else
          tKoo(1,3) = h2
          tKoo(1,4) = -h2
       endif
    end select

  end subroutine parent_tables

!######################################################################
! set_coefficients: w-weighted coupling data of the multiplicity
!######################################################################
  subroutine set_coefficients(mult)

    integer(is), intent(in) :: mult

    real(dp)    :: tJco(2,2), tKco(2,2), tFco(2,2), tWco(2,3,2)
    real(dp)    :: tJov(2,2), tKov(2,2), tFov(2,2), tWov(2,3,2)
    real(dp)    :: tA(2,2), tB(2,2), tJoo(3), tKoo(3,4), tFblk(3)
    integer(is) :: ip, ir, ix, y, xy, z, slot
    real(dp)    :: w

    if (mult == coef_mult) return
    coef_mult = mult
    cJco = 0.0_dp; cKco = 0.0_dp; cJov = 0.0_dp; cKov = 0.0_dp; cJoo = 0.0_dp; cKoo = 0.0_dp
    fco = 0.0_dp; fov = 0.0_dp; Mcomb = 0.0_dp; Gcomb = 0.0_dp; Fblk = 0.0_dp
    do ip = 1, 3
       w = wpi(ip)
       if (w == 0.0_dp) cycle
       call parent_tables(mult, ip, tJco, tKco, tFco, tWco, tJov, tKov, tFov, tWov, tA, tB, tJoo, tKoo, tFblk)
       cJco = cJco + w*tJco
       cKco = cKco + w*tKco
       cJov = cJov + w*tJov
       cKov = cKov + w*tKov
       cJoo = cJoo + w*tJoo
       cKoo = cKoo + w*tKoo
       do ir = 1, 2
          do ix = 1, 2
             if (tFco(ir,ix) /= 0.0_dp) fco(1:nV,ir) = fco(1:nV,ir) + w*tFco(ir,ix)*Fvec(1:nV,ix,ip)
             if (tFov(ir,ix) /= 0.0_dp) fov(1:nC,ir) = fov(1:nC,ir) + w*tFov(ir,ix)*Fhol(1:nC,ix,ip)
             if (tA(ir,ix) /= 0.0_dp) Mcomb(1:nV,1:nV) = Mcomb(1:nV,1:nV) + w*tA(ir,ix)*Mp(1:nV,1:nV,ir,ix)
             if (tB(ir,ix) /= 0.0_dp) Gcomb(1:nC,1:nC) = Gcomb(1:nC,1:nC) + w*tB(ir,ix)*Gp(1:nC,1:nC,ir,ix)
          enddo
          do z = 1, 2
             do xy = 1, 3
                if (tWco(ir,xy,z) /= 0.0_dp) fco(1:nV,ir) = fco(1:nV,ir) + w*tWco(ir,xy,z)*Wv(1:nV,xy,z)
                if (tWov(ir,xy,z) /= 0.0_dp) fov(1:nC,ir) = fov(1:nC,ir) + w*tWov(ir,xy,z)*Wc(1:nC,xy,z)
             enddo
          enddo
       enddo
       do slot = 1, 3
          if (tFblk(slot) /= 0.0_dp) Fblk(1:nV,1:nC,slot) = Fblk(1:nV,1:nC,slot) + w*tFblk(slot)*Fblk3(1:nV,1:nC,ip)
       enddo
    enddo

  end subroutine set_coefficients

!######################################################################
! ext_block_terms: exchange-type C^T x terms of one Q block of the DF
! sweep (Twork holds T(a,Q,col) of the block): for the OV slot columns
! S(:,nocca+j) += ccp cKov(y,x) sum_Q T(:,Q,hO_y) B^Q_{j O_x} (SOMO-row
! corrections in ext_terms) and for the CO slot rows
! S(3:,nocca+j) += ccp cKco(s,x) sum_p (jp|V O_x) X(s,p)
!######################################################################
  subroutine ext_block_terms(blk, nvec, X, S)

    integer(is), intent(in) :: blk, nvec
    real(dp), intent(in)    :: X(nvirb,ncol,nvec)
    real(dp), intent(inout) :: S(nvirb,ncol,nvec)

    real(dp), allocatable :: xs(:,:), t(:,:,:)
    integer(is)           :: Q0, nQl, v, p, ir, ix, iy
    real(dp)              :: t0

    if (nC == 0 .or. nV == 0) return
    t0 = wall_time()
    Q0  = (blk-1)*nQ
    nQl = min(nQ, naux - Q0)

    do iy = 1, 2
       do ix = 1, 2
          if (cKov(iy,ix) == 0.0_dp) cycle
          do v = 1, nvec
             call dgemm('N','N', nvirb, nC, nQl, ccp*cKov(iy,ix), Twork(1,1,nC+iy,v), nvirb, &
                  Boo(1,1,nC+ix,blk), nQ, 1.0_dp, S(1,nocca+1,v), nvirb)
          enddo
       enddo
    enddo

    allocate(xs(nC,nvec), t(nQ,nC,nvec))
    do ir = 1, 2
       if (all(cKco(ir,:) == 0.0_dp)) cycle
       do v = 1, nvec
          do p = 1, nC
             xs(p,v) = X(ir,p,v)
          enddo
       enddo
       call dgemm('N','N', nQ*nC, nvec, nC, 1.0_dp, Boo(1,1,1,blk), nQ*nocca, xs, nC, 0.0_dp, t, nQ*nC)
       do ix = 1, 2
          if (cKco(ir,ix) == 0.0_dp) cycle
          do v = 1, nvec
             call dgemm('T','N', nV, nC, nQl, ccp*cKco(ir,ix), Bvo(Q0+1,1,ix), naux, t(1,1,v), nQ, &
                  1.0_dp, S(3,nocca+1,v), nvirb)
          enddo
       enddo
    enddo
    deallocate(xs, t)
    time_ext = time_ext + wall_time() - t0

  end subroutine ext_block_terms

!######################################################################
! ext_terms: extended-method sigma contributions of a batch
!   X:   compressed vectors (nvirb,ncol,nvec)
!   Xt:  expanded vectors
!   Scv: raw exchange-sweep output of the CV columns, (nvirb,nocca,nvec):
!        Scv(a,i) = sum_{Q,j in C} T_cv(a,Q,j) B^Q_{ji}
!######################################################################
  subroutine ext_terms(nvec, mult, X, Xt, S, Scv)

    integer(is), intent(in) :: nvec, mult
    real(dp), intent(in)    :: X(nvirb,ncol,nvec), Xt(nvirb,ncol,nvec)
    real(dp), intent(inout) :: S(nvirb,ncol,nvec)
    real(dp), intent(in)    :: Scv(nvirb,nocca,nvec)

    real(dp), allocatable :: Ypk(:,:), jY(:,:), jx(:,:), jtot(:,:), Yout(:,:), tmp(:,:)
    real(dp), allocatable :: xs(:,:), xo(:,:), Vcv(:,:,:)
    integer(is)           :: v, p, q, j, ir, ix, iy, k, slot, Q0, nQc, nQl, cv0, b
    integer(is)           :: srow(3), scol(3), krow(4), kcol(4)
    real(dp)              :: t0, fac, sgl
    real(dp), external    :: ddot

    if (.not. ext_ready) call mrsf_error('ext_terms: ext_initialise not called')
    t0 = wall_time()
    call set_coefficients(mult)
    cv0 = nocca + 1
    sgl = merge(1.0_dp, 0.0_dp, mult == 1)
    if (nC == 0 .or. nV == 0) return
    ! slots OO, G, D and the Scv elements of the exchange-type slot terms
    srow = (/ 1_is, 1_is, 2_is /); scol = (/ nC+1, nC+2, nC+1 /)
    krow = (/ 2_is, 1_is, 1_is, 2_is /); kcol = (/ nC+1, nC+2, nC+1, nC+2 /)

    ! ---------------------------------------------------------------
    ! Coulomb pass over Bcv (one streaming pass, Q-blocked so that the
    ! block read by the first dgemm is still cached for the second):
    !   j_Y(Q) = sum_{ai} B^Q_ai y_ai   (from the CV columns)
    !   j_x(Q) = Coulomb-type C^T ix terms of all classes
    !   Yout   = sum_Q B^Q_ai [2 sgl s_J j_Y(Q) + ccp j_x(Q)]
    ! ---------------------------------------------------------------
    allocate(Ypk(nvirb*nC,nvec), jY(naux,nvec), jx(naux,nvec), jtot(naux,nvec))
    allocate(Yout(nvirb*nC,nvec), xs(nC,nvec), xo(nV,nvec), tmp(max(nV,nC),nvec))
    do v = 1, nvec
       call dcopy(nvirb*nC, Xt(1,cv0,v), 1, Ypk(1,v), 1)
    enddo
    jx = 0.0_dp
    do ir = 1, 2
       if (all(cJco(ir,:) == 0.0_dp)) cycle
       do v = 1, nvec
          do p = 1, nC
             xs(p,v) = X(ir,p,v)
          enddo
       enddo
       do ix = 1, 2
          if (cJco(ir,ix) /= 0.0_dp) call dgemm('N','N', naux, nvec, nC, cJco(ir,ix), Bco(1,1,ix), naux, &
               xs, nC, 1.0_dp, jx, naux)
       enddo
    enddo
    do iy = 1, 2
       if (all(cJov(iy,:) == 0.0_dp)) cycle
       do v = 1, nvec
          do q = 1, nV
             xo(q,v) = X(2+q,nC+iy,v)
          enddo
       enddo
       do ix = 1, 2
          if (cJov(iy,ix) /= 0.0_dp) call dgemm('N','N', naux, nvec, nV, cJov(iy,ix), Bvo(1,1,ix), naux, &
               xo, nV, 1.0_dp, jx, naux)
       enddo
    enddo
    do slot = 1, 3
       if (cJoo(slot) == 0.0_dp) cycle
       do v = 1, nvec
          call daxpy(naux, cJoo(slot)*X(srow(slot),scol(slot),v), Bo12, 1, jx(1,v), 1)
       enddo
    enddo
    nQc = max(1_is, min(naux, int(4.0e6_dp / (8.0_dp*real(nvirb,dp)*real(nC,dp)), is)))
    Yout = 0.0_dp
    do Q0 = 0, naux-1, nQc
       nQl = min(nQc, naux - Q0)
       call dgemm('T','N', nQl, nvec, nvirb*nC, 1.0_dp, Bcv(1,1,Q0+1), nvirb*nC, Ypk, nvirb*nC, &
            0.0_dp, jY(Q0+1,1), naux)
       jtot(Q0+1:Q0+nQl,:) = 2.0_dp*sgl*coul_scale*jY(Q0+1:Q0+nQl,:) + ccp*jx(Q0+1:Q0+nQl,:)
       call dgemm('N','N', nvirb*nC, nvec, nQl, 1.0_dp, Bcv(1,1,Q0+1), nvirb*nC, jtot(Q0+1,1), naux, &
            1.0_dp, Yout, nvirb*nC)
    enddo
    do v = 1, nvec
       call daxpy(nvirb*nC, 1.0_dp, Yout(1,v), 1, S(1,cv0,v), 1)
    enddo
    ! Coulomb-type C iy terms
    do ix = 1, 2
       if (any(cJco(:,ix) /= 0.0_dp)) then
          call dgemm('T','N', nC, nvec, naux, ccp, Bco(1,1,ix), naux, jY, naux, 0.0_dp, tmp, max(nV,nC))
          do ir = 1, 2
             if (cJco(ir,ix) == 0.0_dp) cycle
             do v = 1, nvec
                do p = 1, nC
                   S(ir,p,v) = S(ir,p,v) + cJco(ir,ix)*tmp(p,v)
                enddo
             enddo
          enddo
       endif
       if (any(cJov(:,ix) /= 0.0_dp)) then
          call dgemm('T','N', nV, nvec, naux, ccp, Bvo(1,1,ix), naux, jY, naux, 0.0_dp, tmp, max(nV,nC))
          do iy = 1, 2
             if (cJov(iy,ix) == 0.0_dp) cycle
             do v = 1, nvec
                do q = 1, nV
                   S(2+q,nC+iy,v) = S(2+q,nC+iy,v) + cJov(iy,ix)*tmp(q,v)
                enddo
             enddo
          enddo
       endif
    enddo
    do slot = 1, 3
       if (cJoo(slot) == 0.0_dp) cycle
       do v = 1, nvec
          S(srow(slot),scol(slot),v) = S(srow(slot),scol(slot),v) + ccp*cJoo(slot)*dot_product(Bo12, jY(:,v))
       enddo
    enddo
    deallocate(Ypk, jY, jx, jtot, Yout, xs, xo, tmp)

    do v = 1, nvec
       ! ------------------------------------------------------------
       ! exchange-type C iy terms from the sweep output and the CV
       ! exchange term itself
       ! ------------------------------------------------------------
       do ir = 1, 2
          do ix = 1, 2
             if (cKco(ir,ix) == 0.0_dp) cycle
             do p = 1, nC
                S(ir,p,v) = S(ir,p,v) + ccp*cKco(ir,ix)*Scv(ix,p,v)
             enddo
          enddo
       enddo
       do iy = 1, 2
          do ix = 1, 2
             if (cKov(iy,ix) == 0.0_dp) cycle
             do q = 1, nV
                S(2+q,nC+iy,v) = S(2+q,nC+iy,v) + ccp*cKov(iy,ix)*Scv(2+q,nC+ix,v)
             enddo
          enddo
       enddo
       do slot = 1, 3
          do k = 1, 4
             if (cKoo(slot,k) /= 0.0_dp) S(srow(slot),scol(slot),v) = S(srow(slot),scol(slot),v) &
                  + ccp*cKoo(slot,k)*Scv(krow(k),kcol(k),v)
          enddo
       enddo
       do j = 1, nC
          call daxpy(nV, -chf, Scv(3,j,v), 1, S(3,nocca+j,v), 1)
       enddo
       ! ------------------------------------------------------------
       ! Fock-type and SOMO-pair C iy terms
       ! ------------------------------------------------------------
       do ir = 1, 2
          call dgemv('T', nV, nC, ccp, Xt(3,cv0,v), nvirb, fco(1,ir), 1, 1.0_dp, S(ir,1,v), nvirb)
       enddo
       do iy = 1, 2
          call dgemv('N', nV, nC, ccp, Xt(3,cv0,v), nvirb, fov(1,iy), 1, 1.0_dp, S(3,nC+iy,v), 1)
       enddo
       call dgemm('N','N', nV, nC, nV, ccp, Mcomb, nV, Xt(3,cv0,v), nvirb, 1.0_dp, S(3,1,v), nvirb)
       call dgemm('N','N', nV, nC, nC, ccp, Xt(3,cv0,v), nvirb, Gcomb, nC, 1.0_dp, S(3,1,v), nvirb)
       do slot = 1, 3
          if (all(Fblk(:,:,slot) == 0.0_dp)) cycle
          fac = 0.0_dp
          do j = 1, nC
             fac = fac + ddot(nV, Fblk(1,j,slot), 1, Xt(3,nocca+j,v), 1)
          enddo
          S(srow(slot),scol(slot),v) = S(srow(slot),scol(slot),v) + ccp*fac
       enddo
       ! ------------------------------------------------------------
       ! C^T ix terms into the CV columns (the exchange-type parts were
       ! accumulated in the sweep, the Coulomb-type ones through j_x)
       ! ------------------------------------------------------------
       do ir = 1, 2
          call dger(nV, nC, ccp, fco(1,ir), 1, X(ir,1,v), nvirb, S(3,cv0,v), nvirb)
       enddo
       do iy = 1, 2
          call dger(nV, nC, ccp, X(3,nC+iy,v), 1, fov(1,iy), 1, S(3,cv0,v), nvirb)
       enddo
       call dgemm('T','N', nV, nC, nV, ccp, Mcomb, nV, X(3,1,v), nvirb, 1.0_dp, S(3,cv0,v), nvirb)
       call dgemm('N','T', nV, nC, nC, ccp, X(3,1,v), nvirb, Gcomb, nC, 1.0_dp, S(3,cv0,v), nvirb)
       do slot = 1, 3
          if (all(Fblk(:,:,slot) == 0.0_dp)) cycle
          do j = 1, nC
             call daxpy(nV, ccp*X(srow(slot),scol(slot),v), Fblk(1,j,slot), 1, S(3,nocca+j,v), 1)
          enddo
       enddo
       ! transposed exchange-type terms of the slots OO, G, D:
       ! S(3:,cv0+j) += ccp cKoo(slot,k) X(slot) (j O_c|O_r b) for the Scv
       ! element k = (row r, column hO_c)
       do slot = 1, 3
          do k = 1, 4
             if (cKoo(slot,k) == 0.0_dp) cycle
             fac = ccp*cKoo(slot,k)*X(srow(slot),scol(slot),v)
             if (fac == 0.0_dp) cycle
             do j = 1, nC
                call daxpy(nV, fac, NpT(1,j,kcol(k)-nC,krow(k)), 1, S(3,nocca+j,v), 1)
             enddo
          enddo
       enddo
       ! SOMO-row corrections of the transposed OV exchange terms (the
       ! sweep contracted the full hole columns of X~, including the
       ! OO-slot entries X~(1:2,hO_y))
       do iy = 1, 2
          do ix = 1, 2
             if (cKov(iy,ix) == 0.0_dp) cycle
             do k = 1, 2
                if (Xt(k,nC+iy,v) == 0.0_dp) cycle
                do j = 1, nC
                   call daxpy(nV, -ccp*cKov(iy,ix)*Xt(k,nC+iy,v), NpT(1,j,ix,k), 1, S(3,nocca+j,v), 1)
                enddo
             enddo
          enddo
       enddo
       ! shift
       do j = 1, nC
          call daxpy(nV, A_G, X(3,nocca+j,v), 1, S(3,nocca+j,v), 1)
       enddo
    enddo
    time_ext = time_ext + wall_time() - t0

    ! ---------------------------------------------------------------
    ! XC kernel of the singlet CV block (closed-shell kernel installed by
    ! the driver): (jb| f_aa + f_ab |ld) y_ld, all vectors in one pass
    ! ---------------------------------------------------------------
    if (mult == 1 .and. use_kernel .and. .not. kernel_off) then
       t0 = wall_time()
       allocate(Vcv(nC,nV,nvec))
       call xc_kernel_cv(nvec, X, Vcv)
       do v = 1, nvec
          do j = 1, nC
             do b = 1, nV
                S(2+b,nocca+j,v) = S(2+b,nocca+j,v) + Vcv(j,b,v)
             enddo
          enddo
       enddo
       deallocate(Vcv)
       time_xcs = time_xcs + wall_time() - t0
       nxcs_vecs = nxcs_vecs + nvec
    endif

  end subroutine ext_terms

!######################################################################
! ext_diagonal: CV part of the diagonal (columns nocca+1..ncol)
!######################################################################
  subroutine ext_diagonal(mult, d)

    integer(is), intent(in) :: mult
    real(dp), intent(inout) :: d(nvirb,ncol)
    integer(is)             :: j, a

    if (.not. ext_ready) call mrsf_error('ext_diagonal: ext_initialise not called')
    do j = 1, nC
       d(1:2,nocca+j) = 1.0e20_dp
       do a = 3, nvirb
          d(a,nocca+j) = diag_cv(a,j) + A_G
          if (mult == 1) d(a,nocca+j) = d(a,nocca+j) + 2.0_dp*coul_scale*jbjb(a,j)
       enddo
    enddo

  end subroutine ext_diagonal

!######################################################################
! ext_free
!######################################################################
  subroutine ext_free()

    if (allocated(Fp_emb)) deallocate(Fp_emb)
    if (allocated(Fp_cc)) deallocate(Fp_cc)
    if (allocated(Bo12)) deallocate(Bo12)
    if (allocated(diag_cv)) deallocate(diag_cv)
    if (allocated(jbjb)) deallocate(jbjb)
    if (allocated(NpT)) deallocate(NpT)
    if (allocated(Fvec)) deallocate(Fvec, Fhol, Fblk3)
    if (allocated(Wv)) deallocate(Wv, Wc)
    if (allocated(fco)) deallocate(fco, fov, Mcomb, Gcomb, Fblk)
    ext_ready = .false.
    coef_mult = 0

  end subroutine ext_free

end module mrsf_extended
