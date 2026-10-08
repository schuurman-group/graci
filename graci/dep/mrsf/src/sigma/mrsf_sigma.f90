!**********************************************************************
! mrsf_sigma: the MRSF-TDDFT (TDA) sigma vector
!   sigma = U^T [ Fb X~ - X~ Fa - c (ij|ab) X~ ] + pairing(x)
! and, for the extended method, the CV columns
!   sigma_cv = F' Y - Y F' - c (jl|bd) Y + [2 (jb|ld) Y + K^xc Y] + A_G Y
!            + C^T x,   sigma += C y
! The exchange sweep runs once over the paired planes for all columns
! (MRSF and CV) of the batch; its step 2 on the CV columns with all
! hole columns of B^Q_{HH} delivers the CV exchange term and the
! exchange-type coupling terms (see mrsf_extended).
!**********************************************************************
module mrsf_sigma

  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_space
  use mrsf_extended

  implicit none

contains

!######################################################################
! ensure_work: work arrays for a batch of nvec vectors
!######################################################################
  subroutine ensure_work(nvec)

    integer(is), intent(in) :: nvec

    if (nvec > nvmax) then
       if (allocated(Twork)) deallocate(Twork)
       nvmax = nvec
       allocate(Twork(nvirb,nQ,ncol,nvmax), source=0.0_dp)
    endif
    if (store_sp) then
       if (.not. allocated(plane_scr)) allocate(plane_scr(nvirb,nvirb))
       if (.not. allocated(Boo_scr)) allocate(Boo_scr(nQ,nocca,nocca))
    endif

  end subroutine ensure_work

!######################################################################
! sigma_batch: A x for nvec compressed vectors of one (mult, irrep)
! block; irrep < 0 means no irrep projection
!######################################################################
  subroutine sigma_batch(nvec, mult, irrep, x, ax)

    integer(is), intent(in) :: nvec, mult, irrep
    real(dp), intent(in)    :: x(xdim_tot,nvec)
    real(dp), intent(out)   :: ax(xdim_tot,nvec)
    real(dp)                :: t0

    if (.not. ints_loaded) call mrsf_error('integrals not loaded')

    t0 = wall_time()
    call sigma_core(nvec, mult, x, ax)
    call apply_mask(mult, irrep, nvec, ax)
    time_sigma = time_sigma + wall_time() - t0
    nsigma_calls = nsigma_calls + 1
    nsigma_vecs  = nsigma_vecs + nvec

  end subroutine sigma_batch

!######################################################################
! sigma_core: X(nvirb,ncol,nvec) compressed -> S(nvirb,ncol,nvec)
! compressed (masks not applied)
!######################################################################
  subroutine sigma_core(nvec, mult, X, S)

    integer(is), intent(in) :: nvec, mult
    real(dp), intent(in)    :: X(nvirb,ncol,nvec)
    real(dp), intent(out)   :: S(nvirb,ncol,nvec)
    real(dp), allocatable   :: Xt(:,:,:), Scv(:,:,:)
    integer(is)             :: v, cv0

    allocate(Xt(nvirb,ncol,nvec))
    call expand(mult, nvec, X, Xt)
    if (extended) call set_coefficients(mult)

    ! Fock terms: S = Fb_PP X~ - X~ Fa_HH (MRSF columns)
    if (ncol == nocca) then
       call dgemm('N','N', nvirb, nocca*nvec, nvirb, 1.0_dp, FbPP, nvirb, Xt, nvirb, &
            0.0_dp, S, nvirb)
    else
       do v = 1, nvec
          call dgemm('N','N', nvirb, nocca, nvirb, 1.0_dp, FbPP, nvirb, Xt(1,1,v), nvirb, &
               0.0_dp, S(1,1,v), nvirb)
       enddo
    endif
    do v = 1, nvec
       call dgemm('N','N', nvirb, nocca, nocca, -1.0_dp, Xt(1,1,v), nvirb, FaHH, nocca, &
            1.0_dp, S(1,1,v), nvirb)
    enddo

    ! CV columns: F'_VV Y - Y F'_CC
    if (extended) then
       cv0 = nocca + 1
       do v = 1, nvec
          call dgemm('N','N', nvirb, nC, nvirb, 1.0_dp, Fp_emb, nvirb, Xt(1,cv0,v), nvirb, &
               0.0_dp, S(1,cv0,v), nvirb)
          call dgemm('N','N', nvirb, nC, nC, -1.0_dp, Xt(1,cv0,v), nvirb, Fp_cc, nC, &
               1.0_dp, S(1,cv0,v), nvirb)
       enddo
       allocate(Scv(nvirb,nocca,nvec), source=0.0_dp)
    endif

    ! exchange term: S -= c (ij|ab) X~ (and the CV sweep output Scv)
    if (extended) then
       call exchange_contract(nvec, Xt, S, X, Scv)
    else
       call exchange_contract(nvec, Xt, S)
    endif

    ! fold back to the compressed representation
    call fold(mult, nvec, S)

    ! pairing-strength couplings (compressed amplitudes in and out)
    call pairing_terms(mult, nvec, X, S)

    ! extended terms: CV Coulomb/shift/kernel, couplings in both directions
    if (extended) then
       call ext_terms(nvec, mult, X, Xt, S, Scv)
       deallocate(Scv)
    endif

    deallocate(Xt)

  end subroutine sigma_core

!######################################################################
! exchange_contract: S(a,i,v) -= c sum_{jb} (ij|ab) Xt(b,j,v) over the
! MRSF columns.  Paired-plane path: per plane two dsymm calls into
! T(nvirb,nQ,ncol,nvec) for all columns of the batch, then per vector
! one dgemm with the Q-blocked hole-hole integrals.  With the extended
! method (X, Scv present) the CV columns of T are contracted with all
! hole columns of B^Q_{HH} into Scv (unscaled), and the exchange-type
! C^T x terms of the block are accumulated (ext_block_terms).
!######################################################################
  subroutine exchange_contract(nvec, Xt, S, X, Scv)

    integer(is), intent(in)           :: nvec
    real(dp), intent(in)              :: Xt(nvirb,ncol,nvec)
    real(dp), intent(inout)           :: S(nvirb,ncol,nvec)
    real(dp), intent(in), optional    :: X(nvirb,ncol,nvec)
    real(dp), intent(inout), optional :: Scv(nvirb,nocca,nvec)

    integer(is)           :: blk, Q0, Ql, Q, kk, k, v, ldT, ncols
    real(dp), allocatable :: dd(:)
    real(dp)              :: t0
    logical               :: ext

    t0 = wall_time()
    call ensure_work(nvec)
    allocate(dd(nvirb))
    ldT   = nvirb * nQ
    ncols = ncol * nvec
    ext   = present(Scv)

    do blk = 1, nblk
       Q0 = (blk-1)*nQ

       ! step 1: T(a,Q,j,v) = sum_b B^Q_ab Xt(b,j,v)
       if (vv_full) then
          do Ql = 1, nQ
             Q = Q0 + Ql
             if (Q > naux) exit
             if (store_sp) then
                plane_scr = real(Bvv_sp(:,:,Q), dp)
                call dgemm('N','N', nvirb, ncols, nvirb, 1.0_dp, plane_scr, nvirb, &
                     Xt, nvirb, 0.0_dp, Twork(1,Ql,1,1), ldT)
             else
                call dgemm('N','N', nvirb, ncols, nvirb, 1.0_dp, Bvv(1,1,Q), nvirb, &
                     Xt, nvirb, 0.0_dp, Twork(1,Ql,1,1), ldT)
             endif
          enddo
       else
          do kk = 1, npblk
             k = (blk-1)*npblk + kk
             if (k > nplane) exit
             Q  = 2*k - 1
             Ql = 2*kk - 1
             if (store_sp) then
                plane_scr = real(Bvv_sp(:,:,k), dp)
                call plane_apply(ncols, plane_scr, Xt, Q, Ql, ldT, dd)
             else
                call plane_apply(ncols, Bvv(1,1,k), Xt, Q, Ql, ldT, dd)
             endif
          enddo
       endif

       ! step 2: S(a,i,v) -= c sum_{Q,j} T(a,Q,j,v) B^Q_ij
       if (store_sp) then
          Boo_scr = real(Boo_sp(:,:,:,blk), dp)
          do v = 1, nvec
             call dgemm('N','N', nvirb, nocca, nQ*nocca, -chf, Twork(1,1,1,v), nvirb, &
                  Boo_scr, nQ*nocca, 1.0_dp, S(1,1,v), nvirb)
          enddo
       else
          do v = 1, nvec
             call dgemm('N','N', nvirb, nocca, nQ*nocca, -chf, Twork(1,1,1,v), nvirb, &
                  Boo(1,1,1,blk), nQ*nocca, 1.0_dp, S(1,1,v), nvirb)
          enddo
       endif

       ! CV columns: Scv(a,i,v) += sum_{Q, j in C} T(a,Q,nocca+j,v) B^Q_ji
       ! (all holes i; the core-core sub-block of Boo is addressed through
       ! K = nQ*nC rows and ldb = nQ*nocca)
       if (ext) then
          do v = 1, nvec
             call dgemm('N','N', nvirb, nocca, nQ*nC, 1.0_dp, Twork(1,1,nocca+1,v), nvirb, &
                  Boo(1,1,1,blk), nQ*nocca, 1.0_dp, Scv(1,1,v), nvirb)
          enddo
          call ext_block_terms(blk, nvec, X, S)
       endif
    enddo

    deallocate(dd)
    time_exch = time_exch + wall_time() - t0

  end subroutine exchange_contract

!######################################################################
! plane_apply: the two dsymm calls (lower -> Q, upper -> Q+1) for one
! paired plane on ncols columns Xt(nvirb,ncols), plus the diagonal
! correction of the upper plane
!######################################################################
  subroutine plane_apply(ncols, plane, Xt, Q, Ql, ldT, dd)

    integer(is), intent(in) :: ncols, Q, Ql, ldT
    real(dp), intent(in)    :: plane(nvirb,nvirb)
    real(dp), intent(in)    :: Xt(nvirb,ncols)
    real(dp), intent(inout) :: dd(nvirb)
    integer(is)             :: a

    call dsymm('L','L', nvirb, ncols, 1.0_dp, plane, nvirb, Xt, nvirb, &
         0.0_dp, Twork(1,Ql,1,1), ldT)

    if (Q + 1 <= naux) then
       call dsymm('L','U', nvirb, ncols, 1.0_dp, plane, nvirb, Xt, nvirb, &
            0.0_dp, Twork(1,Ql+1,1,1), ldT)
       do a = 1, nvirb
          dd(a) = Dall(Q+1,Pmap(a)) - Dall(Q,Pmap(a))
       enddo
       call diag_correct(nvirb, nQ, ncols, Twork, Ql+1, dd, Xt)
    endif

  end subroutine plane_apply

!######################################################################
! diag_correct: T(:,Ql,c) += dd * Xt(:,c) for all columns c
!######################################################################
  subroutine diag_correct(n1, n2, ncols, T, Ql, dd, Xt)

    integer(is), intent(in) :: n1, n2, ncols, Ql
    real(dp), intent(inout) :: T(n1,n2,ncols)
    real(dp), intent(in)    :: dd(n1), Xt(n1,ncols)
    integer(is)             :: c

    !$omp parallel do private(c)
    do c = 1, ncols
       T(:,Ql,c) = T(:,Ql,c) + dd(:) * Xt(:,c)
    enddo
    !$omp end parallel do

  end subroutine diag_correct

!######################################################################
! pairing_terms: spin-pairing couplings between the CO and OV
! amplitudes (compressed representation), added to S
!######################################################################
  subroutine pairing_terms(mult, nvec, X, S)

    integer(is), intent(in) :: mult, nvec
    real(dp), intent(in)    :: X(nvirb,ncol,nvec)
    real(dp), intent(inout) :: S(nvirb,ncol,nvec)

    real(dp), allocatable :: tC1(:,:), tC2(:,:), oC1(:,:), oC2(:,:), oV1(:,:), oV2(:,:)
    real(dp)              :: sgn, kc, ko, kv
    integer(is)           :: v, hO1, hO2

    sgn = merge(1.0_dp, -1.0_dp, mult == 1)
    kc = sgn * spc(1)
    ko = sgn * spc(2)
    kv = sgn * spc(3)
    hO1 = nC + 1
    hO2 = nC + 2

    allocate(tC1(max(nC,1_is),nvec), tC2(max(nC,1_is),nvec))
    allocate(oC1(max(nC,1_is),nvec), oC2(max(nC,1_is),nvec))
    allocate(oV1(max(nV,1_is),nvec), oV2(max(nV,1_is),nvec))
    oC1 = 0.0_dp; oC2 = 0.0_dp; oV1 = 0.0_dp; oV2 = 0.0_dp

    do v = 1, nvec
       tC1(1:nC,v) = X(1,1:nC,v)
       tC2(1:nC,v) = X(2,1:nC,v)
    enddo

    ! core -> O1/O2 outputs
    if (nC > 0) then
       if (nV > 0) then
          call dgemm('N','N', nC, nvec, nV,  kv, Hp, nC, X(3,hO2,1), xdim_tot, 0.0_dp, oC1, nC)
          call dgemm('N','N', nC, nvec, nV, -kv, Hp, nC, X(3,hO1,1), xdim_tot, 0.0_dp, oC2, nC)
       endif
       call dgemm('N','N', nC, nvec, nC,  kc, Gp(1,1,2,2), nC, tC1, nC, 1.0_dp, oC1, nC)
       call dgemm('N','N', nC, nvec, nC, -kc, Gp(1,1,1,2), nC, tC2, nC, 1.0_dp, oC1, nC)
       call dgemm('N','N', nC, nvec, nC, -kc, Gp(1,1,2,1), nC, tC1, nC, 1.0_dp, oC2, nC)
       call dgemm('N','N', nC, nvec, nC,  kc, Gp(1,1,1,1), nC, tC2, nC, 1.0_dp, oC2, nC)
    endif

    ! O1/O2 -> virtual outputs
    if (nV > 0) then
       if (nC > 0) then
          call dgemm('T','N', nV, nvec, nC, -kv, Hp, nC, tC2, nC, 0.0_dp, oV1, nV)
          call dgemm('T','N', nV, nvec, nC,  kv, Hp, nC, tC1, nC, 0.0_dp, oV2, nV)
       endif
       call dgemm('N','N', nV, nvec, nV,  ko, Mp(1,1,2,2), nV, X(3,hO1,1), xdim_tot, 1.0_dp, oV1, nV)
       call dgemm('N','N', nV, nvec, nV, -ko, Mp(1,1,1,2), nV, X(3,hO2,1), xdim_tot, 1.0_dp, oV1, nV)
       call dgemm('N','N', nV, nvec, nV,  ko, Mp(1,1,1,1), nV, X(3,hO2,1), xdim_tot, 1.0_dp, oV2, nV)
       call dgemm('N','N', nV, nvec, nV, -ko, Mp(1,1,2,1), nV, X(3,hO1,1), xdim_tot, 1.0_dp, oV2, nV)
    endif

    do v = 1, nvec
       S(1,1:nC,v) = S(1,1:nC,v) + oC1(1:nC,v)
       S(2,1:nC,v) = S(2,1:nC,v) + oC2(1:nC,v)
       S(3:nvirb,hO1,v) = S(3:nvirb,hO1,v) + oV1(1:nV,v)
       S(3:nvirb,hO2,v) = S(3:nvirb,hO2,v) + oV2(1:nV,v)
    enddo

    deallocate(tC1, tC2, oC1, oC2, oV1, oV2)

  end subroutine pairing_terms

!######################################################################
! diagonal: exact diagonal of A for the given multiplicity (packed);
! redundant slots are set to a large value
!######################################################################
  subroutine diagonal(mult, d)

    integer(is), intent(in) :: mult
    real(dp), intent(out)   :: d(nvirb,ncol)
    real(dp)                :: s, kc, ko
    integer(is)             :: i, a, hO1, hO2

    if (.not. ints_loaded) call mrsf_error('integrals not loaded')

    s  = merge(1.0_dp, -1.0_dp, mult == 1)
    kc = s * spc(1)
    ko = s * spc(2)
    hO1 = nC + 1
    hO2 = nC + 2

    d(:,1:nocca) = diag0
    do i = 1, nC
       d(1,i) = d(1,i) + kc * Gp(i,i,2,2)
       d(2,i) = d(2,i) + kc * Gp(i,i,1,1)
    enddo
    do a = 1, nV
       d(2+a,hO1) = d(2+a,hO1) + ko * Mp(a,a,2,2)
       d(2+a,hO2) = d(2+a,hO2) + ko * Mp(a,a,1,1)
    enddo

    ! folded OO slot
    d(1,hO1) = 0.5_dp * (diag0(1,hO1) + diag0(2,hO2)) + s * chf * K12

    ! redundant slots
    d(2,hO2) = 1.0e20_dp
    if (mult == 3) then
       d(1,hO2) = 1.0e20_dp
       d(2,hO1) = 1.0e20_dp
    endif

    ! CV columns
    if (extended) call ext_diagonal(mult, d)

  end subroutine diagonal

end module mrsf_sigma
