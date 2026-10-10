!**********************************************************************
! mrsf_extended: extended MRSF-TDDFT (EMRSF-TDDFT; Oh, Kim, Jung, Choi,
! Lee, ChemRxiv 2026). The response vector gains the core-to-virtual
! (CV) amplitudes y(a,i), a in P (SOMO rows masked), i in C, of the
! closed-shell configuration G = |C Cbar O1 O1bar|, stored as the extra
! columns nocca+1..ncol of the (nvirb x ncol) vector [x | y].
!
! Sigma contributions handled here (sigma_core does the Fock terms and
! the exchange sweep, which also delivers the exchange-type coupling
! terms through Scv = sum_Q T_cv B^Q_{HH}, T_cv = sum_b B^Q_ab y_bj):
!   CV block:  2 (jb|ld) y   (singlet),  + A_G y,  + K^xc y  (singlet, DFT);
!              -c_H (jl|bd) y and the F' terms are added by sigma_core
!   couplings: C y into the MRSF slots and C^T x into the CV columns,
!              scaled by c_cp, with the closed forms pinned in
!              ~/calculations/mrsf_dev/emrsf_ref.py (COEF):
!   singlet CV CSFs (t=+1); CV slot (j in C -> b in V), MRSF slot (a,i):
!     CO1 (p):  -2 (jb|O2p) + (jp|bO2) - d_pj F_{O2 b}
!     CO2 (p):  -d_pj (O2O1|O2b)
!     O1V (q):  +d_qb (jO1|O2O1)
!     O2V (q):  +2 (jb|qO1) - (jO1|qb) - d_qb F_{j O1}
!     CV (p,q): -d_pj (qO1|O2b) + d_qb (jO1|O2p)
!     OO:       +(jO1|O2b) - 2 (jb|O2O1)     (compressed OO slot)
!     G:        +sqrt2 F_{jb};  D: 0
!   triplet CV CSFs (t=-1): Coulomb-type terms vanish; O1V, O2V(F), CV(B)
!     change sign: CO1 (jp|bO2) - d F;  CO2 -d W;  O1V -d W;
!     O2V -(jO1|qb)*(-1) + d F -> +(jO1|qb) + d_qb F;  CV -A - B;  OO +K.
! F = F^DFT (eq. 17), the closed-shell KS matrix of G with the triplet
! orbitals; F' = F^DFT + delta (1-c_H)[(pp|O2O2) - (pp|O1O1)] (eq. 16).
!**********************************************************************
module mrsf_extended

  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_space
  use mrsf_xcgrid, only: xc_ready, ao_cached, xc_closed_shell_init, xc_kernel_cv

  implicit none

  real(dp), allocatable :: jbjb(:,:)      ! (nvirb,nC) (jb|jb)
  real(dp), allocatable :: Np11T(:,:), Np12T(:,:)   ! (nV,nC) (jO1|O1b), (jO1|O2b)
  ! coupling coefficients of the current multiplicity
  real(dp) :: cJ1, cK1, cF1, cW2, cW3, cJ4, cK4, cF4, cA5, cB5, cK6, cJ6, cG
  integer(is) :: coef_mult = 0
  ! kernel-free operator A0 for the inner correction equations of the
  ! Davidson solver (set around the inner sigma batches)
  logical :: kernel_off = .false.
  ! F' correction: 1 = diagonal (published), 2 = covariant VV and CC blocks
  ! of (1 - c_H)(J[D_O2] - J[D_O1]) (Part XII)
  integer(is) :: fprime_mode = 1

contains

!######################################################################
! ext_initialise: F', shift, diagonal, SOMO-integral vectors
!   fdft(nmo,nmo): F^DFT in the MO basis
!######################################################################
  subroutine ext_initialise(fdft, ccp1, use_kernel1, fprime1)

    real(dp), intent(in)    :: fdft(nmo,nmo)
    real(dp), intent(in)    :: ccp1
    logical, intent(in)     :: use_kernel1
    integer(is), intent(in) :: fprime1

    real(dp), allocatable :: corr(:), dO1(:), dO2(:), Dpart(:,:), dQ(:), Jvv(:,:), Jcc(:,:)
    integer(is)           :: p, a, j, b, blk, Ql, Q
    real(dp)              :: o2o2o1o1

    if (.not. extended) call mrsf_error('ext_initialise: extended flag not set')
    if (.not. ints_loaded) call mrsf_error('ext_initialise: integrals not loaded')
    if (.not. allocated(Bcv)) call mrsf_error('ext_initialise: particle-core block missing')
    call ext_free()

    ccp = ccp1
    use_kernel = use_kernel1
    if (use_kernel) then
       if (.not. xc_ready) call mrsf_error('ext_initialise: the grid kernel is not initialised')
       if (.not. ao_cached) call mrsf_error('ext_initialise: the extended method needs the '// &
            'AO values cached on the grid (raise mem_budget)')
       call xc_closed_shell_init()
    endif

    ! diagonal correction (1 - c_H)[(pp|O2O2) - (pp|O1O1)]
    allocate(corr(nmo), dO1(naux), dO2(naux))
    dO1 = Dall(:,iO1)
    dO2 = Dall(:,iO2)
    call dgemv('T', naux, nmo, 1.0_dp, Dall, naux, dO2, 1, 0.0_dp, corr, 1)
    call dgemv('T', naux, nmo, -1.0_dp, Dall, naux, dO1, 1, 1.0_dp, corr, 1)
    corr = (1.0_dp - chf) * corr
    o2o2o1o1 = dot_product(dO2, dO1)
    if (fprime1 /= 1 .and. fprime1 /= 2) call mrsf_error('ext_initialise: fprime must be 1 or 2')
    fprime_mode = fprime1
    if (fprime_mode == 2) then
       if (store_sp) call mrsf_error('ext_initialise: the covariant F'' needs double-precision integrals')
       allocate(dQ(naux), Jvv(nvirb,nvirb), Jcc(max(nC,1_is),max(nC,1_is)))
       dQ = (1.0_dp - chf) * (dO2 - dO1)
       call ext_coulomb_blocks(dQ, Jvv, Jcc)
    endif

    ! F' blocks (V x V embedded in the particle layout, C x C)
    allocate(Fp_emb(nvirb,nvirb), source=0.0_dp)
    allocate(Fp_cc(max(nC,1_is),max(nC,1_is)), source=0.0_dp)
    do b = 3, nvirb
       do a = 3, nvirb
          Fp_emb(a,b) = fdft(Pmap(a),Pmap(b))
          if (fprime_mode == 2) Fp_emb(a,b) = Fp_emb(a,b) + Jvv(a,b)
       enddo
       if (fprime_mode == 1) Fp_emb(b,b) = Fp_emb(b,b) + corr(Pmap(b))
    enddo
    do j = 1, nC
       do p = 1, nC
          Fp_cc(p,j) = fdft(Hmap(p),Hmap(j))
          if (fprime_mode == 2) Fp_cc(p,j) = Fp_cc(p,j) + Jcc(p,j)
       enddo
       if (fprime_mode == 1) Fp_cc(j,j) = Fp_cc(j,j) + corr(Hmap(j))
    enddo
    if (allocated(dQ)) deallocate(dQ, Jvv, Jcc)

    ! Fock-type coupling vectors / block
    allocate(Fcv(max(nV,1_is),max(nC,1_is)), fO2V(max(nV,1_is)), fCO1(max(nC,1_is)))
    do j = 1, nC
       do b = 1, nV
          Fcv(b,j) = fdft(Pmap(2+b),Hmap(j))
       enddo
       fCO1(j) = fdft(Hmap(j),iO1)
    enddo
    do b = 1, nV
       fO2V(b) = fdft(iO2,Pmap(2+b))
    enddo

    ! B^Q_{O1O2} and the SOMO-integral vectors
    allocate(Bo12(naux), w1(max(nV,1_is)), w2(max(nC,1_is)))
    do Q = 1, naux
       blk = (Q-1)/nQ + 1
       Ql  = Q - (blk-1)*nQ
       Bo12(Q) = Boo(Ql,nC+1,nC+2,blk)
    enddo
    w1 = 0.0_dp; w2 = 0.0_dp
    if (nV > 0) call dgemv('T', naux, nV, 1.0_dp, Bvo(1,1,2), naux, Bo12, 1, 0.0_dp, w1, 1)
    if (nC > 0) call dgemv('T', naux, nC, 1.0_dp, Bco(1,1,1), naux, Bo12, 1, 0.0_dp, w2, 1)

    ! transposed SOMO-pair blocks (contiguous columns for the rank-1 updates)
    allocate(Np11T(max(nV,1_is),max(nC,1_is)), Np12T(max(nV,1_is),max(nC,1_is)))
    do j = 1, nC
       do b = 1, nV
          Np11T(b,j) = Np(j,b,1,1)
          Np12T(b,j) = Np(j,b,1,2)
       enddo
    enddo

    ! shift of the CV block (eq. 13)
    A_G = FbPP(1,1) - FaHH(nC+2,nC+2) - chf * o2o2o1o1

    ! CV diagonal without the shift: F'_bb - F'_jj - c_H (jj|bb) [+ 2 (jb|jb) singlet]
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
    deallocate(corr, dO1, dO2, Dpart)

    ext_ready = .true.
    if (verbose) then
       write(6,'(/,2x,a)') 'Extended MRSF-TDDFT initialised'
       write(6,'(2x,a,i0,a,i0,a,f8.4,a,f12.6,a,l1)') 'CV columns = ', nC, ', xdim_tot = ', &
            xdim_tot, ', c_cp = ', ccp, ', A_G = ', A_G, ' Ha, kernel = ', use_kernel
       if (fprime_mode == 2) then
          write(6,'(2x,a)') "F' correction: covariant (full VV and CC blocks)"
       else
          write(6,'(2x,a)') "F' correction: diagonal (published form)"
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
! ext_free
!######################################################################
  subroutine ext_free()

    if (allocated(Fp_emb)) deallocate(Fp_emb)
    if (allocated(Fp_cc)) deallocate(Fp_cc)
    if (allocated(Fcv)) deallocate(Fcv, fO2V, fCO1)
    if (allocated(Bo12)) deallocate(Bo12, w1, w2)
    if (allocated(diag_cv)) deallocate(diag_cv)
    if (allocated(jbjb)) deallocate(jbjb)
    if (allocated(Np11T)) deallocate(Np11T, Np12T)
    ext_ready = .false.
    coef_mult = 0

  end subroutine ext_free

!######################################################################
! set_coefficients: coupling coefficients of a multiplicity
!######################################################################
  subroutine set_coefficients(mult)

    integer(is), intent(in) :: mult

    if (mult == coef_mult) return
    coef_mult = mult
    if (mult == 1) then
       cJ1 = -2.0_dp; cK1 = 1.0_dp; cF1 = -1.0_dp
       cW2 = -1.0_dp
       cW3 = 1.0_dp
       cJ4 = 2.0_dp; cK4 = -1.0_dp; cF4 = -1.0_dp
       cA5 = -1.0_dp; cB5 = 1.0_dp
       cK6 = 1.0_dp; cJ6 = -2.0_dp
       cG  = sqrt(2.0_dp)
    else
       cJ1 = 0.0_dp; cK1 = 1.0_dp; cF1 = -1.0_dp
       cW2 = -1.0_dp
       cW3 = -1.0_dp
       cJ4 = 0.0_dp; cK4 = 1.0_dp; cF4 = 1.0_dp
       cA5 = -1.0_dp; cB5 = -1.0_dp
       cK6 = 1.0_dp; cJ6 = 0.0_dp
       cG  = 0.0_dp
    endif

  end subroutine set_coefficients

!######################################################################
! ext_block_terms: exchange-type C^T x contributions evaluated inside
! the Q-block loop of the exchange sweep, right after step 2 of block
! blk (T of the MRSF columns and the hole-hole block of blk are the
! data just touched):
!   O2V K: sigma_cv(b,j) += ccp cK4 sum_Q B^Q_{jO1} T(b,Q,hO2)
!   CO1 K: sigma_cv(b,j) += ccp cK1 sum_Q B^Q_{bO2} t(Q,j),
!          t(Q,j) = sum_{p in C} B^Q_{jp} x_{p O1}
! X: compressed vectors (nvirb,ncol,nvec); S: sigma, CV columns updated
!######################################################################
  subroutine ext_block_terms(blk, nvec, X, S)

    integer(is), intent(in) :: blk, nvec
    real(dp), intent(in)    :: X(nvirb,ncol,nvec)
    real(dp), intent(inout) :: S(nvirb,ncol,nvec)

    real(dp), allocatable :: xCO1(:,:), t(:,:,:)
    integer(is)           :: Q0, nQl, v, hO1, hO2, p
    real(dp)              :: t0

    if (nC == 0 .or. nV == 0) return
    t0 = wall_time()
    Q0  = (blk-1)*nQ
    nQl = min(nQ, naux - Q0)
    hO1 = nC + 1
    hO2 = nC + 2

    ! O2V exchange-type term, transposed (SOMO-row corrections in ext_terms)
    do v = 1, nvec
       call dgemm('N','N', nvirb, nC, nQl, ccp*cK4, Twork(1,1,hO2,v), nvirb, &
            Boo(1,1,hO1,blk), nQ, 1.0_dp, S(1,nocca+1,v), nvirb)
    enddo

    ! CO1 exchange-type term, transposed
    allocate(xCO1(nC,nvec), t(nQ,nC,nvec))
    do v = 1, nvec
       do p = 1, nC
          xCO1(p,v) = X(1,p,v)
       enddo
    enddo
    call dgemm('N','N', nQ*nC, nvec, nC, 1.0_dp, Boo(1,1,1,blk), nQ*nocca, xCO1, nC, &
         0.0_dp, t, nQ*nC)
    do v = 1, nvec
       call dgemm('T','N', nV, nC, nQl, ccp*cK1, Bvo(Q0+1,1,2), naux, t(1,1,v), nQ, &
            1.0_dp, S(3,nocca+1,v), nvirb)
    enddo
    deallocate(xCO1, t)
    time_ext = time_ext + wall_time() - t0

  end subroutine ext_block_terms

!######################################################################
! ext_terms: all remaining extended contributions, added to the folded
! compressed sigma S (MRSF columns) and to the CV columns.
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
    real(dp), allocatable :: xCO1(:,:), xO2V(:,:), Vcv(:,:,:)
    integer(is)           :: v, p, q, j, hO1, hO2, Q0, nQc, nQl, cv0, b
    real(dp)              :: t0, fac, sgl

    if (.not. ext_ready) call mrsf_error('ext_terms: ext_initialise not called')
    t0 = wall_time()
    call set_coefficients(mult)
    hO1 = nC + 1
    hO2 = nC + 2
    cv0 = nocca + 1
    sgl = merge(1.0_dp, 0.0_dp, mult == 1)
    if (nC == 0 .or. nV == 0) return

    ! ---------------------------------------------------------------
    ! Coulomb pass over Bcv (one streaming pass, Q-blocked so that the
    ! block read by the first dgemm is still cached for the second):
    !   j_Y(Q) = sum_{ai} B^Q_ai y_ai   (from the CV columns)
    !   j_x(Q) = cJ1 sum_p B^Q_{pO2} x_{pO1} + cJ4 sum_q B^Q_{qO1} x_{O2 q}
    !          + cJ6 B^Q_{O1O2} x_OO       (Coulomb-type C^T x terms)
    !   Yout   = sum_Q B^Q_ai [2 sgl j_Y(Q) + ccp j_x(Q)]
    ! ---------------------------------------------------------------
    allocate(Ypk(nvirb*nC,nvec), jY(naux,nvec), jx(naux,nvec), jtot(naux,nvec))
    allocate(Yout(nvirb*nC,nvec), xCO1(nC,nvec), xO2V(nV,nvec), tmp(max(nV,nC),nvec))
    do v = 1, nvec
       call dcopy(nvirb*nC, Xt(1,cv0,v), 1, Ypk(1,v), 1)
       do p = 1, nC
          xCO1(p,v) = X(1,p,v)
       enddo
       do q = 1, nV
          xO2V(q,v) = X(2+q,hO2,v)
       enddo
    enddo
    call dgemm('N','N', naux, nvec, nC, cJ1, Bco(1,1,2), naux, xCO1, nC, 0.0_dp, jx, naux)
    call dgemm('N','N', naux, nvec, nV, cJ4, Bvo(1,1,1), naux, xO2V, nV, 1.0_dp, jx, naux)
    do v = 1, nvec
       call daxpy(naux, cJ6*X(1,hO1,v), Bo12, 1, jx(1,v), 1)
    enddo
    nQc = max(1_is, min(naux, int(4.0e6_dp / (8.0_dp*real(nvirb,dp)*real(nC,dp)), is)))
    Yout = 0.0_dp
    do Q0 = 0, naux-1, nQc
       nQl = min(nQc, naux - Q0)
       call dgemm('T','N', nQl, nvec, nvirb*nC, 1.0_dp, Bcv(1,1,Q0+1), nvirb*nC, Ypk, nvirb*nC, &
            0.0_dp, jY(Q0+1,1), naux)
       jtot(Q0+1:Q0+nQl,:) = 2.0_dp*sgl*jY(Q0+1:Q0+nQl,:) + ccp*jx(Q0+1:Q0+nQl,:)
       call dgemm('N','N', nvirb*nC, nvec, nQl, 1.0_dp, Bcv(1,1,Q0+1), nvirb*nC, jtot(Q0+1,1), naux, &
            1.0_dp, Yout, nvirb*nC)
    enddo
    do v = 1, nvec
       call daxpy(nvirb*nC, 1.0_dp, Yout(1,v), 1, S(1,cv0,v), 1)
    enddo
    ! Coulomb-type C y terms
    if (cJ1 /= 0.0_dp) then
       call dgemm('T','N', nC, nvec, naux, ccp*cJ1, Bco(1,1,2), naux, jY, naux, 0.0_dp, tmp, max(nV,nC))
       do v = 1, nvec
          do p = 1, nC
             S(1,p,v) = S(1,p,v) + tmp(p,v)
          enddo
       enddo
    endif
    if (cJ4 /= 0.0_dp) then
       call dgemm('T','N', nV, nvec, naux, ccp*cJ4, Bvo(1,1,1), naux, jY, naux, 0.0_dp, tmp, max(nV,nC))
       do v = 1, nvec
          do q = 1, nV
             S(2+q,hO2,v) = S(2+q,hO2,v) + tmp(q,v)
          enddo
       enddo
    endif
    if (cJ6 /= 0.0_dp) then
       do v = 1, nvec
          S(1,hO1,v) = S(1,hO1,v) + ccp*cJ6*dot_product(Bo12, jY(:,v))
       enddo
    endif
    deallocate(Ypk, jY, jx, jtot, Yout, xCO1, xO2V, tmp)

    do v = 1, nvec
       ! ------------------------------------------------------------
       ! exchange-type C y terms from the sweep output and the CV
       ! exchange term itself
       ! ------------------------------------------------------------
       do p = 1, nC
          S(1,p,v) = S(1,p,v) + ccp*cK1*Scv(2,p,v)
       enddo
       do q = 1, nV
          S(2+q,hO2,v) = S(2+q,hO2,v) + ccp*cK4*Scv(2+q,hO1,v)
       enddo
       S(1,hO1,v) = S(1,hO1,v) + ccp*cK6*Scv(2,hO1,v)
       do j = 1, nC
          call daxpy(nV, -chf, Scv(3,j,v), 1, S(3,nocca+j,v), 1)
       enddo
       ! ------------------------------------------------------------
       ! Fock-type and SOMO-pair C y terms
       ! ------------------------------------------------------------
       ! CO1 F: -d_pj F_{O2 b};  CO2 W: -d_pj (O2O1|O2b)
       call dgemv('T', nV, nC, ccp*cF1, Xt(3,cv0,v), nvirb, fO2V, 1, 1.0_dp, S(1,1,v), nvirb)
       call dgemv('T', nV, nC, ccp*cW2, Xt(3,cv0,v), nvirb, w1, 1, 1.0_dp, S(2,1,v), nvirb)
       ! O2V F: d_qb F_{jO1};  O1V W: d_qb (jO1|O2O1)
       call dgemv('N', nV, nC, ccp*cF4, Xt(3,cv0,v), nvirb, fCO1, 1, 1.0_dp, S(3,hO2,v), 1)
       call dgemv('N', nV, nC, ccp*cW3, Xt(3,cv0,v), nvirb, w2, 1, 1.0_dp, S(3,hO1,v), 1)
       ! CV A: -d_pj (qO1|O2b) -> Mp12 Y;  CV B: d_qb (jO1|O2p) -> Y Gp12
       call dgemm('N','N', nV, nC, nV, ccp*cA5, Mp(1,1,1,2), nV, Xt(3,cv0,v), nvirb, &
            1.0_dp, S(3,1,v), nvirb)
       call dgemm('N','N', nV, nC, nC, ccp*cB5, Xt(3,cv0,v), nvirb, Gp(1,1,1,2), nC, &
            1.0_dp, S(3,1,v), nvirb)
       ! G: sqrt2 sum_{jb} F_{jb} y_bj  (singlet only, slot (1,hO2))
       if (cG /= 0.0_dp) then
          fac = 0.0_dp
          do j = 1, nC
             fac = fac + dot_product(Fcv(:,j), Xt(3:nvirb,nocca+j,v))
          enddo
          S(1,hO2,v) = S(1,hO2,v) + ccp*cG*fac
       endif
       ! ------------------------------------------------------------
       ! C^T x terms into the CV columns (the exchange-type parts were
       ! accumulated in the sweep)
       ! ------------------------------------------------------------
       ! Fock type: CO1 F, O2V F, G
       call dger(nV, nC, ccp*cF1, fO2V, 1, X(1,1,v), nvirb, S(3,cv0,v), nvirb)
       call dger(nV, nC, ccp*cF4, X(3,hO2,v), 1, fCO1, 1, S(3,cv0,v), nvirb)
       if (cG /= 0.0_dp) then
          do j = 1, nC
             call daxpy(nV, ccp*cG*X(1,hO2,v), Fcv(1,j), 1, S(3,nocca+j,v), 1)
          enddo
       endif
       ! SOMO pair: CO2 W, O1V W, CV A, CV B, OO K
       call dger(nV, nC, ccp*cW2, w1, 1, X(2,1,v), nvirb, S(3,cv0,v), nvirb)
       call dger(nV, nC, ccp*cW3, X(3,hO1,v), 1, w2, 1, S(3,cv0,v), nvirb)
       call dgemm('T','N', nV, nC, nV, ccp*cA5, Mp(1,1,1,2), nV, X(3,1,v), nvirb, &
            1.0_dp, S(3,cv0,v), nvirb)
       call dgemm('N','T', nV, nC, nC, ccp*cB5, X(3,1,v), nvirb, Gp(1,1,1,2), nC, &
            1.0_dp, S(3,cv0,v), nvirb)
       do j = 1, nC
          call daxpy(nV, ccp*cK6*X(1,hO1,v), Np12T(1,j), 1, S(3,nocca+j,v), 1)
       enddo
       ! SOMO-row corrections of the transposed O2V exchange term (the
       ! sweep contracted the full hole column hO2 of X~, including the
       ! OO-slot entries X~(1,hO2) and X~(2,hO2))
       do j = 1, nC
          call daxpy(nV, -ccp*cK4*Xt(1,hO2,v), Np11T(1,j), 1, S(3,nocca+j,v), 1)
          call daxpy(nV, -ccp*cK4*Xt(2,hO2,v), Np12T(1,j), 1, S(3,nocca+j,v), 1)
       enddo
       ! shift A_G
       do j = 1, nC
          call daxpy(nV, A_G, X(3,nocca+j,v), 1, S(3,nocca+j,v), 1)
       enddo
    enddo
    time_ext = time_ext + wall_time() - t0

    ! ---------------------------------------------------------------
    ! XC kernel of the singlet CV block (closed-shell kernel at the G
    ! density): (jb| f_aa + f_ab |ld) y_ld, all vectors in one pass
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
          if (mult == 1) d(a,nocca+j) = d(a,nocca+j) + 2.0_dp*jbjb(a,j)
       enddo
    enddo

  end subroutine ext_diagonal

end module mrsf_extended
