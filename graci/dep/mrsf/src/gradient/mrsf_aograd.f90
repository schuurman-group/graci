!**********************************************************************
! mrsf_aograd: AO-side contraction of the density-fitted two-electron
! gradient with hole-width factor pairs, one auxiliary block at a time.
! For a family with Gamma^P = L^P R^T + R L^P^T (L^P: nao x nk per aux
! function P, R: nao x nk) the per-AO and per-aux vectors are
!   t(mu,x) += sum_{nu,P} (d_x mu nu|P) Gamma^P_{mu nu}
!            = sum_{P,o} L^P(mu,o) [sum_nu (d_x mu nu|P) R(nu,o)]      (term A)
!            + sum_o R(mu,o) [sum_{nu,P} (d_x mu nu|P) L^P(nu,o)]       (term B)
!   u(P,x)  += 2 sum_{mu,o} L^P(mu,o) [sum_nu (mu nu|d_x P) R(nu,o)]
! The derivative integrals arrive in libcint's buffer layout: ip1(mu,nu,P,x)
! = (d_x mu nu|P) (derivative on mu) and the packed lower triangle
! ip2p(ij,P,x) = (mu nu|d_x P), ij = mu(mu-1)/2 + nu, mu >= nu (PySCF
! aosym s2ij). The half-transforms of term A and of the aux derivative
! are one small dgemm per (P, x) inside an OpenMP loop (BLAS runs
! single-threaded there), with (mu nu|d_x P) unpacked into a
! thread-private buffer; term B is one dgemm per x with N = nao*nP.
!**********************************************************************
module mrsf_aograd

  use mrsf_constants
  use mrsf_io

  implicit none

  real(dp), allocatable, private :: M1T(:), N1T(:), M2(:), LfT(:)
  integer(is), private           :: n_m1 = 0, n_n1 = 0, n_m2 = 0, n_lft = 0
  real(dp)                       :: time_ao = 0.0_dp, tpart_ao(4) = 0.0_dp
  integer(is)                    :: nao_calls = 0

contains

  subroutine aograd_block(nao, nk, nP, nf, ip1, ip2p, Lf, Rf, t, u)

    integer(is), intent(in) :: nao, nk, nP, nf
    real(dp), intent(in)    :: ip1(nao,nao,nP,3), ip2p(nao*(nao+1)/2,nP,3)
    real(dp), intent(in)    :: Lf(nao,nk,nP,nf), Rf(nao,nk,nf)
    real(dp), intent(inout) :: t(nao,3,nf)
    real(dp), intent(out)   :: u(nP,3,nf)

    real(dp) :: tstart

    tstart = wall_time()
    call ensure(M1T, n_m1, 3*nk*nao*nP*nf)
    call ensure(N1T, n_n1, 3*nk*nao*nP*nf)
    call ensure(M2, n_m2, 3*nao*nk)
    call ensure(LfT, n_lft, nao*nP*nk)
    call half_transforms(nao, nk, nP, nf, ip1, ip2p, Rf, M1T, N1T)
    call contract(nao, nk, nP, nf, ip1, Lf, Rf, t, u, M1T, N1T, M2, LfT)
    time_ao = time_ao + wall_time() - tstart
    nao_calls = nao_calls + 1

  end subroutine aograd_block

  subroutine ensure(arr, n, nneed)
    real(dp), allocatable, intent(inout) :: arr(:)
    integer(is), intent(inout) :: n
    integer(is), intent(in)    :: nneed
    if (nneed > n) then
       if (allocated(arr)) deallocate(arr)
       allocate(arr(nneed))
       n = nneed
    endif
  end subroutine ensure

  subroutine aograd_free()
    if (allocated(M1T)) deallocate(M1T)
    if (allocated(N1T)) deallocate(N1T)
    if (allocated(M2)) deallocate(M2)
    if (allocated(LfT)) deallocate(LfT)
    n_m1 = 0; n_n1 = 0; n_m2 = 0; n_lft = 0
  end subroutine aograd_free

!######################################################################
! half_transforms: M1T(o,mu,P,x,f) = sum_nu R_f(nu,o) (d_x mu nu|P),
! N1T(o,mu,P,x,f) = sum_nu R_f(nu,o) (mu nu|d_x P)
!######################################################################
  subroutine half_transforms(nao, nk, nP, nf, ip1, ip2p, Rf, M1T, N1T)

    integer(is), intent(in) :: nao, nk, nP, nf
    real(dp), intent(in)    :: ip1(nao,nao,nP,3), ip2p(nao*(nao+1)/2,nP,3), Rf(nao,nk,nf)
    real(dp), intent(out)   :: M1T(nk,nao,nP,3,nf), N1T(nk,nao,nP,3,nf)

    real(dp), allocatable :: buf(:,:)
    integer(is) :: x, P, f, mu, nu, ij
    real(dp) :: t0

    t0 = wall_time()
    !$omp parallel private(buf,P,f,mu,nu,ij)
    allocate(buf(nao,nao))
    !$omp do collapse(2) schedule(dynamic,1)
    do x = 1, 3
       do P = 1, nP
          ij = 0
          do mu = 1, nao
             do nu = 1, mu
                ij = ij + 1
                buf(nu,mu) = ip2p(ij,P,x)
             enddo
          enddo
          do mu = 1, nao
             do nu = mu+1, nao
                buf(nu,mu) = buf(mu,nu)
             enddo
          enddo
          do f = 1, nf
             call dgemm('T','T', nk, nao, nao, 1.0_dp, Rf(1,1,f), nao, ip1(1,1,P,x), nao, &
                  0.0_dp, M1T(1,1,P,x,f), nk)
             call dgemm('T','N', nk, nao, nao, 1.0_dp, Rf(1,1,f), nao, buf, nao, &
                  0.0_dp, N1T(1,1,P,x,f), nk)
          enddo
       enddo
    enddo
    !$omp end do
    deallocate(buf)
    !$omp end parallel
    tpart_ao(1) = tpart_ao(1) + wall_time() - t0

  end subroutine half_transforms

!######################################################################
! contract: term A, term B and the aux-derivative reductions
!######################################################################
  subroutine contract(nao, nk, nP, nf, ip1, Lf, Rf, t, u, M1T, N1T, M2, LfT)

    integer(is), intent(in) :: nao, nk, nP, nf
    real(dp), intent(in)    :: ip1(nao,nao,nP,3), Lf(nao,nk,nP,nf), Rf(nao,nk,nf)
    real(dp), intent(in)    :: M1T(nk,nao,nP,3,nf), N1T(nk,nao,nP,3,nf)
    real(dp), intent(inout) :: t(nao,3,nf)
    real(dp), intent(out)   :: u(nP,3,nf)
    real(dp), intent(inout) :: M2(nao,nk,3), LfT(nao,nP,nk)

    integer(is) :: f, x, P, mu, nu, o
    real(dp) :: acc, t0, t1

    u = 0.0_dp
    do f = 1, nf
       t0 = wall_time()
       ! term A and the aux derivative from the half-transforms
       !$omp parallel do private(mu,P,o,acc) collapse(2)
       do x = 1, 3
          do mu = 1, nao
             acc = 0.0_dp
             do P = 1, nP
                do o = 1, nk
                   acc = acc + Lf(mu,o,P,f)*M1T(o,mu,P,x,f)
                enddo
             enddo
             t(mu,x,f) = t(mu,x,f) + acc
          enddo
       enddo
       !$omp end parallel do
       !$omp parallel do private(P,mu,o,acc) collapse(2)
       do x = 1, 3
          do P = 1, nP
             acc = 0.0_dp
             do mu = 1, nao
                do o = 1, nk
                   acc = acc + Lf(mu,o,P,f)*N1T(o,mu,P,x,f)
                enddo
             enddo
             u(P,x,f) = 2.0_dp*acc
          enddo
       enddo
       !$omp end parallel do
       ! term B: M2(mu,o,x) = sum_{nu,P} (d_x mu nu|P) L^P(nu,o)
       !$omp parallel do private(P,nu) collapse(2)
       do o = 1, nk
          do P = 1, nP
             do nu = 1, nao
                LfT(nu,P,o) = Lf(nu,o,P,f)
             enddo
          enddo
       enddo
       !$omp end parallel do
       t1 = wall_time(); tpart_ao(3) = tpart_ao(3) + t1 - t0; t0 = t1
       do x = 1, 3
          call dgemm('N','N', nao, nk, nao*nP, 1.0_dp, ip1(1,1,1,x), nao, LfT, nao*nP, &
               0.0_dp, M2(1,1,x), nao)
       enddo
       t1 = wall_time(); tpart_ao(2) = tpart_ao(2) + t1 - t0; t0 = t1
       !$omp parallel do private(mu,o,acc) collapse(2)
       do x = 1, 3
          do mu = 1, nao
             acc = 0.0_dp
             do o = 1, nk
                acc = acc + Rf(mu,o,f)*M2(mu,o,x)
             enddo
             t(mu,x,f) = t(mu,x,f) + acc
          enddo
       enddo
       !$omp end parallel do
       tpart_ao(3) = tpart_ao(3) + wall_time() - t0
    enddo

  end subroutine contract

  subroutine aograd_timings()
    write(6,'(/,2x,a)') 'AO derivative contraction timings (s)'
    write(6,'(2x,a,f10.3)') 'half-transforms (per-P)  : ', tpart_ao(1)
    write(6,'(2x,a,f10.3)') 'term B dgemm             : ', tpart_ao(2)
    write(6,'(2x,a,f10.3)') 'reductions/transposes    : ', tpart_ao(3)
    write(6,'(2x,a,f10.3,a,i0,a)') 'total                    : ', time_ao, ' (', nao_calls, ' blocks)'
    flush(6)
  end subroutine aograd_timings

end module mrsf_aograd

!**********************************************************************
! C-bound entry points
!**********************************************************************
module mrsf_aograd_interface

  use iso_c_binding
  use mrsf_constants
  use mrsf_aograd

  implicit none

contains

  subroutine mrsf_aograd_block_c(nao1, nk1, nP1, nf1, ip1, ip2p, Lf, Rf, t, u) &
       bind(c, name='mrsf_aograd_block')
    integer(is), intent(in) :: nao1, nk1, nP1, nf1
    real(dp), intent(in)    :: ip1(*), ip2p(*), Lf(*), Rf(*)
    real(dp), intent(inout) :: t(*)
    real(dp), intent(out)   :: u(*)
    call aograd_block(nao1, nk1, nP1, nf1, ip1, ip2p, Lf, Rf, t, u)
  end subroutine mrsf_aograd_block_c

  subroutine mrsf_aograd_free_c() bind(c, name='mrsf_aograd_free')
    call aograd_free()
  end subroutine mrsf_aograd_free_c

  subroutine mrsf_aograd_timings_c() bind(c, name='mrsf_aograd_timings')
    call aograd_timings()
  end subroutine mrsf_aograd_timings_c

end module mrsf_aograd_interface
