!**********************************************************************
! scf_df: Coulomb and exchange matrices from the packed density-fitting
! tensor of PySCF (lower-Cholesky convention, C-ordered (naux, npair)
! array, i.e. the Fortran (npair, naux) matrix whose columns are the
! upper-packed AO matrices B^Q, pair index c(c-1)/2 + r for r <= c). The
! tensor stays in Python's memory; the library keeps its address.
!   J = sum_Q B^Q (B^Q . D)                 (two dgemv passes)
!   Y^Q = B^Q X (dsymm per Q, OpenMP over Q), stored in the slab
!         Y(nao, naux, ka) so that every column set is one matrix,
!   K over the first kb columns and over the columns kb+1..ka by dsyrk
!**********************************************************************
module scf_df

  use scf_constants
  use scf_io
  use iso_c_binding
  !$ use omp_lib

  implicit none

  type(c_ptr)           :: bptr = c_null_ptr
  integer(is)           :: naux_d = 0, npair_d = 0, nao_d = 0
  real(dp), allocatable :: Yslab(:)
  real(dp)              :: time_j = 0.0_dp, time_k = 0.0_dp, flop_k = 0.0_dp, time_y = 0.0_dp

  external :: dgemv, dsymm, dsyrk

contains

  subroutine df_set(nao, naux, cptr)
    integer(is), intent(in) :: nao, naux
    type(c_ptr), intent(in) :: cptr
    nao_d = nao
    naux_d = naux
    npair_d = nao*(nao+1)/2
    bptr = cptr
    if (allocated(Yslab)) deallocate(Yslab)
  end subroutine df_set

  subroutine df_free()
    bptr = c_null_ptr
    if (allocated(Yslab)) deallocate(Yslab)
    naux_d = 0
  end subroutine df_free

!######################################################################
! df_jk: J of the total density D; wantk = 0: no K; 1: K1 = K of the
! first kb columns of X, K2 = K of the columns kb+1..ka (ROKS); 2: K1 =
! K of all ka columns (closed shell or equal spins)
!######################################################################
  subroutine df_jk(nao, D, ka, kb, X, wantk, J, K1, K2)

    integer(is), intent(in) :: nao, ka, kb, wantk
    real(dp), intent(in)    :: D(nao,nao), X(nao,ka)
    real(dp), intent(out)   :: J(nao,nao), K1(nao,nao), K2(nao,nao)

    real(dp), pointer     :: B(:,:)
    real(dp), allocatable :: dpk(:), jq(:), jpk(:), P(:,:)
    integer(is) :: r, c, off, q, naux, npair
    integer(ib) :: need
    real(dp)    :: t0, t1, tj, acc

    naux = naux_d
    npair = npair_d
    call c_f_pointer(bptr, B, [npair, naux])

    ! Coulomb: packed density with doubled off-diagonal elements; with
    ! exchange, j^Q is formed in the same pass over the tensor as Y^Q
    t0 = wall_time()
    allocate(dpk(npair), jq(naux), jpk(npair))
    do c = 1, nao
       off = c*(c-1)/2
       do r = 1, c-1
          dpk(off+r) = 2.0_dp*D(r,c)
       enddo
       dpk(off+c) = D(c,c)
    enddo
    K1 = 0.0_dp
    K2 = 0.0_dp
    if (wantk == 0 .or. ka == 0) then
       call dgemv('T', npair, naux, 1.0_dp, B, npair, dpk, 1, 0.0_dp, jq, 1)
    else
       ! exchange: Y^Q = B^Q X for all Q (and j^Q), then dsyrk over the column sets
       need = int(nao,ib)*int(naux,ib)*int(ka,ib)
       if (allocated(Yslab)) then
          if (size(Yslab,kind=ib) < need) deallocate(Yslab)
       endif
       if (.not. allocated(Yslab)) allocate(Yslab(need))
       !$omp parallel private(P, q, c, off, r, acc)
       call blas_sequential()
       allocate(P(nao,nao))
       !$omp do schedule(dynamic)
       do q = 1, naux
          acc = 0.0_dp
          do c = 1, nao
             off = c*(c-1)/2
             do r = 1, c
                P(r,c) = B(off+r,q)
                acc = acc + B(off+r,q)*dpk(off+r)
             enddo
          enddo
          jq(q) = acc
          call dsymm('L', 'U', nao, ka, 1.0_dp, P, nao, X, nao, 0.0_dp, Yslab(int(q-1,ib)*nao+1), nao*naux)
       enddo
       !$omp end do
       deallocate(P)
       call blas_restore()
       !$omp end parallel
       time_y = time_y + wall_time() - t0
    endif
    t1 = wall_time()
    call dgemv('N', npair, naux, 1.0_dp, B, npair, jq, 1, 0.0_dp, jpk, 1)
    do c = 1, nao
       off = c*(c-1)/2
       do r = 1, c
          J(r,c) = jpk(off+r)
          J(c,r) = jpk(off+r)
       enddo
    enddo
    deallocate(dpk, jq, jpk)
    tj = wall_time() - t1
    time_j = time_j + tj
    if (wantk == 0 .or. ka == 0) then
       time_j = time_j + t1 - t0
       return
    endif
    if (wantk == 1) then
       if (kb > 0) call dsyrk('U', 'N', nao, naux*kb, 1.0_dp, Yslab, nao, 0.0_dp, K1, nao)
       if (ka > kb) call dsyrk('U', 'N', nao, naux*(ka-kb), 1.0_dp, Yslab(int(nao,ib)*naux*kb+1), &
            nao, 0.0_dp, K2, nao)
       call fill_lower(nao, K1)
       call fill_lower(nao, K2)
    else
       call dsyrk('U', 'N', nao, naux*ka, 1.0_dp, Yslab, nao, 0.0_dp, K1, nao)
       call fill_lower(nao, K1)
    endif
    flop_k = flop_k + 2.0_dp*real(nao,dp)**2*real(ka,dp)*real(naux,dp) &
         + real(nao,dp)**2*real(naux,dp)*real(ka,dp)
    time_k = time_k + wall_time() - t0 - tj

  end subroutine df_jk

  subroutine fill_lower(n, A)
    integer(is), intent(in) :: n
    real(dp), intent(inout) :: A(n,n)
    integer(is) :: i, j
    do j = 1, n
       do i = j+1, n
          A(i,j) = A(j,i)
       enddo
    enddo
  end subroutine fill_lower

end module scf_df
