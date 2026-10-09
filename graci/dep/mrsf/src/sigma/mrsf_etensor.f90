!**********************************************************************
! mrsf_etensor: explicit exchange tensor E_ij(a,b) = (ij|ab) for the
! active hole pairs i <= j (hole set H minus the frozen core), a,b in P,
! pre-contracted once from the DF blocks, E_ij = sum_Q B^Q_ij B^Q_ab.
! Every E_ij is a symmetric (nvirb x nvirb) matrix; two pairs share one
! square plane (pair 2k-1 in the lower triangle, pair 2k in the strict
! upper triangle, all diagonals in Ediag): the layout of the DF planes
! with the pair index in place of Q. The exchange term of the sigma
! vector, sigma(a,i) -= c sum_{j,b} (ij|ab) X(b,j), becomes one dsymm
! per (pair, column block) on the vectors in their (nvirb,ncol,nvec)
! layout: 2 nocca_act nvirb^2 ncol flops per vector instead of
! 2 naux nvirb^2 ncol for the DF sweep, and no T buffer. The CV columns
! of the extended method and the exchange-type coupling by-products
! (Scv, the O2V and CO1 block terms of mrsf_extended) come from the
! same planes.
!**********************************************************************
module mrsf_etensor

  use mrsf_constants
  use mrsf_global
  use mrsf_io

  implicit none

  logical                  :: etensor_ready = .false.
  integer(is)              :: etensor_nfc = -1
  integer(is)              :: npair_e = 0, nplane_e = 0
  integer(is), allocatable :: pair_i(:), pair_j(:)   ! local hole indices of pair p (i before j in active order)
  integer(is), allocatable :: pair_of(:,:)           ! (nocca,nocca) pair index of two active holes (0: none)
  real(dp), allocatable    :: Eplane(:,:,:)          ! (nvirb,nvirb,nplane_e)
  real(dp), allocatable    :: Ediag(:,:)             ! (nvirb,npair_e)
  real(dp), allocatable    :: Eo2(:,:,:)             ! (nV,nC,nC) E_{jp}(2+q,2) = (jp|V_q O2)
  real(dp)                 :: etensor_gb = 0.0_dp
  real(dp)                 :: time_ebuild = 0.0_dp, time_etens = 0.0_dp, flop_etens = 0.0_dp
  integer(is)              :: netens_vecs = 0

contains

!######################################################################
! etensor_bytes: memory of the tensor for the current active hole set
!######################################################################
  function etensor_bytes() result(bytes)

    real(dp)    :: bytes
    integer(is) :: n, np

    n  = nocca_act
    np = n*(n+1)/2
    bytes = 8.0_dp * (real((np+1)/2,dp) * real(nvirb,dp)**2 + real(nvirb,dp)*real(np,dp))
    if (extended) bytes = bytes + 8.0_dp*real(max(nV,1_is),dp)*real(max(nC,1_is),dp)**2

  end function etensor_bytes

!######################################################################
! etensor_free
!######################################################################
  subroutine etensor_free()

    if (allocated(pair_i)) deallocate(pair_i, pair_j)
    if (allocated(pair_of)) deallocate(pair_of)
    if (allocated(Eplane)) deallocate(Eplane)
    if (allocated(Ediag)) deallocate(Ediag)
    if (allocated(Eo2)) deallocate(Eo2)
    etensor_ready = .false.
    etensor_nfc   = -1
    npair_e  = 0
    nplane_e = 0
    etensor_gb = 0.0_dp

  end subroutine etensor_free

!######################################################################
! etensor_build: E from the loaded DF blocks (planes, Boo, Dall), both
! symmetries used: packed lower triangles (a >= b) of the planes times
! the hole-pair columns (i <= j) of B^Q_HH, in batches of pairs with
! the Q loop inside (flops ~ nocca_act^2 naux nvirb^2 / 2)
!######################################################################
  subroutine etensor_build()

    integer(is), parameter :: nqs_max = 64
    integer(is)            :: n, ii, jj, p, pp, k, a, b, q, j, npab, npb, nb, pb0, pb1
    integer(is)            :: blk, Q0, nQl, qs0, nqs, Ql, pr
    integer(is), allocatable :: off(:)
    real(dp), allocatable  :: Etmp(:,:), Apk(:,:), Bho(:,:)
    real(dp)               :: t0

    if (.not. ints_loaded) call mrsf_error('etensor_build: integrals not loaded')
    t0 = wall_time()
    call etensor_free()

    ! active hole pairs i <= j (active order: active core, O1, O2)
    n = nocca_act
    npair_e  = n*(n+1)/2
    nplane_e = (npair_e + 1)/2
    allocate(pair_i(npair_e), pair_j(npair_e), pair_of(nocca,nocca))
    pair_of = 0
    p = 0
    do jj = 1, n
       do ii = 1, jj
          p = p + 1
          pair_i(p) = hact(ii)
          pair_j(p) = hact(jj)
          pair_of(hact(ii),hact(jj)) = p
          pair_of(hact(jj),hact(ii)) = p
       enddo
    enddo

    allocate(Eplane(nvirb,nvirb,nplane_e), source=0.0_dp)
    allocate(Ediag(nvirb,npair_e), source=0.0_dp)

    ! packed lower-triangle column offsets (column-major, a >= b)
    npab = nvirb*(nvirb+1)/2
    allocate(off(nvirb))
    off(1) = 0
    do b = 2, nvirb
       off(b) = off(b-1) + nvirb - b + 2
    enddo

    ! pair batches: Etmp at most ~0.5 GB
    npb = max(1_is, min(npair_e, int(6.0e7_dp/real(npab,dp), is)))
    allocate(Etmp(npab,npb), Apk(npab,nqs_max), Bho(nqs_max,npb))

    do pb0 = 1, npair_e, npb
       pb1 = min(pb0 + npb - 1, npair_e)
       nb  = pb1 - pb0 + 1
       Etmp(:,1:nb) = 0.0_dp
       do blk = 1, nblk
          Q0  = (blk-1)*nQ
          nQl = min(nQ, naux - Q0)
          do qs0 = 1, nQl, nqs_max
             nqs = min(nqs_max, nQl - qs0 + 1)
             call expand_block(blk, qs0, nqs, off, npab, Apk)
             do pp = 1, nb
                p = pb0 + pp - 1
                if (store_sp) then
                   do Ql = 1, nqs
                      Bho(Ql,pp) = real(Boo_sp(qs0+Ql-1,pair_i(p),pair_j(p),blk), dp)
                   enddo
                else
                   do Ql = 1, nqs
                      Bho(Ql,pp) = Boo(qs0+Ql-1,pair_i(p),pair_j(p),blk)
                   enddo
                endif
             enddo
             call dgemm('N','N', npab, nb, nqs, 1.0_dp, Apk, npab, Bho, nqs_max, &
                  1.0_dp, Etmp, npab)
          enddo
       enddo
       ! scatter the batch into the planes
       do pp = 1, nb
          p = pb0 + pp - 1
          k = (p+1)/2
          do b = 1, nvirb
             Ediag(b,p) = Etmp(off(b)+1,pp)
          enddo
          if (mod(p,2_is) == 1) then
             do b = 1, nvirb
                Eplane(b:nvirb,b,k) = Etmp(off(b)+1:off(b)+nvirb-b+1,pp)
             enddo
          else
             do b = 1, nvirb
                do a = b+1, nvirb
                   Eplane(b,a,k) = Etmp(off(b)+a-b+1,pp)
                enddo
             enddo
          endif
       enddo
    enddo
    deallocate(Etmp, Apk, Bho, off)

    ! (jp|V_q O2) for the CO1 block term of the extended method
    if (extended) then
       allocate(Eo2(max(nV,1_is),max(nC,1_is),max(nC,1_is)), source=0.0_dp)
       do p = 1, nC
          do j = 1, nC
             pr = pair_of(j,p)
             if (pr == 0) cycle
             k = (pr+1)/2
             if (mod(pr,2_is) == 1) then
                do q = 1, nV
                   Eo2(q,j,p) = Eplane(2+q,2,k)
                enddo
             else
                do q = 1, nV
                   Eo2(q,j,p) = Eplane(2,2+q,k)
                enddo
             endif
          enddo
       enddo
    endif

    etensor_ready = .true.
    etensor_nfc   = nfc
    etensor_gb    = etensor_bytes()/1.0e9_dp
    time_ebuild   = time_ebuild + wall_time() - t0
    if (verbose) write(6,'(2x,a,i0,a,f8.3,a,f8.2,a)') 'exchange tensor (ij|ab) built for ', &
         npair_e, ' active hole pairs: ', etensor_gb, ' GB, ', wall_time() - t0, ' s'

  end subroutine etensor_build

!######################################################################
! expand_block: packed lower triangles (a >= b, column-major) of the
! planes B^Q for Q = Q0+qs0 .. Q0+qs0+nqs-1 of block blk into A
!######################################################################
  subroutine expand_block(blk, qs0, nqs, off, npab, Apk)

    integer(is), intent(in) :: blk, qs0, nqs, npab
    integer(is), intent(in) :: off(nvirb)
    real(dp), intent(out)   :: Apk(npab,nqs)
    integer(is)             :: Q0, Ql, Q, k, a, b
    real(dp), allocatable   :: pl(:,:)

    Q0 = (blk-1)*nQ
    allocate(pl(nvirb,nvirb))
    if (vv_full) then
       do Ql = 1, nqs
          Q = Q0 + qs0 + Ql - 1
          if (store_sp) then
             pl = real(Bvv_sp(:,:,Q), dp)
          else
             pl = Bvv(:,:,Q)
          endif
          do b = 1, nvirb
             Apk(off(b)+1:off(b)+nvirb-b+1,Ql) = pl(b:nvirb,b)
          enddo
       enddo
    else
       ! paired planes: Q = 2k-1 lower, Q+1 strict upper; qs0 is odd
       do Ql = 1, nqs, 2
          Q = Q0 + qs0 + Ql - 1
          k = (Q+1)/2
          if (store_sp) then
             pl = real(Bvv_sp(:,:,k), dp)
          else
             pl = Bvv(:,:,k)
          endif
          do b = 1, nvirb
             Apk(off(b)+1:off(b)+nvirb-b+1,Ql) = pl(b:nvirb,b)
          enddo
          if (Ql+1 <= nqs .and. Q+1 <= naux) then
             do b = 1, nvirb
                Apk(off(b)+1,Ql+1) = Dall(Q+1,Pmap(b))
                do a = b+1, nvirb
                   Apk(off(b)+a-b+1,Ql+1) = pl(b,a)
                enddo
             enddo
          else if (Ql+1 <= nqs) then
             Apk(:,Ql+1) = 0.0_dp
          endif
       enddo
    endif
    deallocate(pl)

  end subroutine expand_block

!######################################################################
! etensor_contract: exchange term of the sigma vectors from the tensor
!   S(:,i,v) -= c sum_j E_ij Xt(:,j,v)                (MRSF columns)
! and, for the extended method (Scv, X present):
!   Scv(:,i,v) += sum_{j in C} E_ij Y_j(:,v)          (unscaled)
!   S(:,nocca+j,v) += c_o2v E_{j,O1} Xt(:,hO2,v)       (O2V block term)
!   S(3:,nocca+j,v) += c_co1 sum_p (jp|V O2) X(1,p,v)  (CO1 block term)
! Xt: expanded vectors; X: compressed vectors (O1 row of the CO1 term)
!######################################################################
  subroutine etensor_contract(nvec, Xt, S, c_o2v, c_co1, X, Scv)

    integer(is), intent(in)           :: nvec
    real(dp), intent(in)              :: Xt(nvirb,ncol,nvec)
    real(dp), intent(inout)           :: S(nvirb,ncol,nvec)
    real(dp), intent(in)              :: c_o2v, c_co1
    real(dp), intent(in), optional    :: X(nvirb,ncol,nvec)
    real(dp), intent(inout), optional :: Scv(nvirb,nocca,nvec)

    integer(is)           :: k, p1, p2, j, v, p, ldv
    real(dp), allocatable :: dd(:), xCO1(:,:)
    real(dp)              :: t0

    if (.not. etensor_ready) call mrsf_error('etensor_contract: tensor not built')
    t0 = wall_time()
    allocate(dd(nvirb))
    ldv = nvirb*ncol

    do k = 1, nplane_e
       p1 = 2*k - 1
       call pair_apply(nvec, Eplane(1,1,k), 'L', p1, .false., dd, Xt, S, c_o2v, Scv)
       p2 = 2*k
       if (p2 <= npair_e) then
          dd = Ediag(:,p2) - Ediag(:,p1)
          call pair_apply(nvec, Eplane(1,1,k), 'U', p2, .true., dd, Xt, S, c_o2v, Scv)
       endif
    enddo

    ! CO1 block term
    if (present(Scv) .and. present(X) .and. c_co1 /= 0.0_dp .and. nC > 0 .and. nV > 0) then
       allocate(xCO1(nC,nvec))
       do v = 1, nvec
          do p = 1, nC
             xCO1(p,v) = X(1,p,v)
          enddo
       enddo
       do j = 1, nC
          if (frozen_hole(j)) cycle
          call dgemm('N','N', nV, nvec, nC, c_co1, Eo2(1,j,1), nV*nC, xCO1, nC, &
               1.0_dp, S(3,nocca+j,1), ldv)
       enddo
       flop_etens = flop_etens + 2.0_dp*real(nV,dp)*real(nC,dp)*real(nC_act,dp)*real(nvec,dp)
       deallocate(xCO1)
    endif

    deallocate(dd)
    time_etens  = time_etens + wall_time() - t0
    netens_vecs = netens_vecs + nvec

  end subroutine etensor_contract

!######################################################################
! pair_apply: all products of one pair (i,j) stored in plane (uplo)
!######################################################################
  subroutine pair_apply(nvec, plane, uplo, p, upper, dd, Xt, S, c_o2v, Scv)

    integer(is), intent(in)           :: nvec, p
    real(dp), intent(in)              :: plane(nvirb,nvirb)
    character(len=1), intent(in)      :: uplo
    logical, intent(in)               :: upper
    real(dp), intent(in)              :: dd(nvirb)
    real(dp), intent(in)              :: Xt(nvirb,ncol,nvec)
    real(dp), intent(inout)           :: S(nvirb,ncol,nvec)
    real(dp), intent(in)              :: c_o2v
    real(dp), intent(inout), optional :: Scv(nvirb,nocca,nvec)

    integer(is) :: i, j, ldv, ldc, hO1, hO2, nprod

    i   = pair_i(p)
    j   = pair_j(p)
    ldv = nvirb*ncol
    ldc = nvirb*nocca
    hO1 = nC + 1
    hO2 = nC + 2
    nprod = 0

    ! MRSF columns
    call dsymm('L', uplo, nvirb, nvec, -chf, plane, nvirb, Xt(1,j,1), ldv, 1.0_dp, S(1,i,1), ldv)
    if (upper) call diag_update(nvec, -chf, dd, Xt(1,j,1), ldv, S(1,i,1), ldv)
    nprod = nprod + 1
    if (i /= j) then
       call dsymm('L', uplo, nvirb, nvec, -chf, plane, nvirb, Xt(1,i,1), ldv, 1.0_dp, S(1,j,1), ldv)
       if (upper) call diag_update(nvec, -chf, dd, Xt(1,i,1), ldv, S(1,j,1), ldv)
       nprod = nprod + 1
    endif

    if (present(Scv)) then
       ! CV columns: Scv(:,i) += E_ij Y_j (j in C), Scv(:,j) += E_ij Y_i (i in C)
       if (j <= nC) then
          call dsymm('L', uplo, nvirb, nvec, 1.0_dp, plane, nvirb, Xt(1,nocca+j,1), ldv, &
               1.0_dp, Scv(1,i,1), ldc)
          if (upper) call diag_update(nvec, 1.0_dp, dd, Xt(1,nocca+j,1), ldv, Scv(1,i,1), ldc)
          nprod = nprod + 1
       endif
       if (i /= j .and. i <= nC) then
          call dsymm('L', uplo, nvirb, nvec, 1.0_dp, plane, nvirb, Xt(1,nocca+i,1), ldv, &
               1.0_dp, Scv(1,j,1), ldc)
          if (upper) call diag_update(nvec, 1.0_dp, dd, Xt(1,nocca+i,1), ldv, Scv(1,j,1), ldc)
          nprod = nprod + 1
       endif
       ! O2V block term: pair (c, O1), c in C
       if (j == hO1 .and. i <= nC .and. c_o2v /= 0.0_dp) then
          call dsymm('L', uplo, nvirb, nvec, c_o2v, plane, nvirb, Xt(1,hO2,1), ldv, &
               1.0_dp, S(1,nocca+i,1), ldv)
          if (upper) call diag_update(nvec, c_o2v, dd, Xt(1,hO2,1), ldv, S(1,nocca+i,1), ldv)
          nprod = nprod + 1
       endif
    endif
    flop_etens = flop_etens + 2.0_dp*real(nprod,dp)*real(nvirb,dp)**2*real(nvec,dp)

  end subroutine pair_apply

!######################################################################
! diag_update: Sc(:,v) += alpha dd(:) Xc(:,v), columns with leading
! dimensions ldx (input) and lds (output)
!######################################################################
  subroutine diag_update(nvec, alpha, dd, Xc, ldx, Sc, lds)

    integer(is), intent(in) :: nvec, ldx, lds
    real(dp), intent(in)    :: alpha, dd(nvirb)
    real(dp), intent(in)    :: Xc(ldx,*)
    real(dp), intent(inout) :: Sc(lds,*)
    integer(is)             :: v

    do v = 1, nvec
       Sc(1:nvirb,v) = Sc(1:nvirb,v) + alpha * dd(1:nvirb) * Xc(1:nvirb,v)
    enddo

  end subroutine diag_update

!######################################################################
! etensor_kx: S(:,i,v) += scale sum_j E_ij Z(:,j,v) on (nvirb,nocca,
! nvec) arrays (Z-vector operator and stability check; full hole set)
!######################################################################
  subroutine etensor_kx(nvec, scale, Z, S)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: scale
    real(dp), intent(in)    :: Z(nvirb,nocca,nvec)
    real(dp), intent(inout) :: S(nvirb,nocca,nvec)

    integer(is)           :: k, p, pp, i, j, ld
    real(dp), allocatable :: dd(:)
    real(dp)              :: t0

    if (.not. etensor_ready) call mrsf_error('etensor_kx: tensor not built')
    if (etensor_nfc /= 0) call mrsf_error('etensor_kx: the tensor must cover the full hole set')
    t0 = wall_time()
    allocate(dd(nvirb))
    ld = nvirb*nocca
    do k = 1, nplane_e
       do pp = 1, 2
          p = 2*k - 2 + pp
          if (p > npair_e) exit
          i = pair_i(p)
          j = pair_j(p)
          if (pp == 2) dd = Ediag(:,p) - Ediag(:,p-1)
          call dsymm('L', merge('L','U',pp == 1), nvirb, nvec, scale, Eplane(1,1,k), nvirb, &
               Z(1,j,1), ld, 1.0_dp, S(1,i,1), ld)
          if (pp == 2) call diag_update(nvec, scale, dd, Z(1,j,1), ld, S(1,i,1), ld)
          if (i /= j) then
             call dsymm('L', merge('L','U',pp == 1), nvirb, nvec, scale, Eplane(1,1,k), nvirb, &
                  Z(1,i,1), ld, 1.0_dp, S(1,j,1), ld)
             if (pp == 2) call diag_update(nvec, scale, dd, Z(1,i,1), ld, S(1,j,1), ld)
             flop_etens = flop_etens + 2.0_dp*real(nvirb,dp)**2*real(nvec,dp)
          endif
          flop_etens = flop_etens + 2.0_dp*real(nvirb,dp)**2*real(nvec,dp)
       enddo
    enddo
    deallocate(dd)
    time_etens  = time_etens + wall_time() - t0
    netens_vecs = netens_vecs + nvec

  end subroutine etensor_kx

end module mrsf_etensor
