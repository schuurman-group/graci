!**********************************************************************
! mrsf_density: unrelaxed state and transition density matrices from
! MRSF amplitudes (spin-summed, MO basis, rho(p,q) = <bra|E_pq|ket>).
! Stateless: the orbital classes are rebuilt from the occupation vector.
!**********************************************************************
module mrsf_density

  use mrsf_constants
  use mrsf_io
  use mrsf_space, only: classify_orbitals

  implicit none

contains

!######################################################################
! tdm_pairs: rho(:,:,ip) = <bra I_ip|E_pq|ket J_ip> (+ reference
! occupations if add_ref), I from xb, J from xk
!######################################################################
  subroutine tdm_pairs(nmo1, occ1, mult, xdim1, npairs, nb, nk, ipairs, xb, xk, rho, add_ref)

    integer(is), intent(in) :: nmo1, mult, xdim1, npairs, nb, nk
    real(dp), intent(in)    :: occ1(nmo1)
    integer(is), intent(in) :: ipairs(2,npairs)
    real(dp), intent(in)    :: xb(xdim1,nb), xk(xdim1,nk)
    real(dp), intent(out)   :: rho(nmo1,nmo1,npairs)
    logical, intent(in)     :: add_ref

    integer(is), allocatable :: hmap1(:), pmap1(:)
    real(dp), allocatable    :: XI(:,:), XJ(:,:), PP(:,:), HH(:,:), Y(:,:), W(:,:)
    integer(is)              :: nc1, nv1, io1, io2, nocca1, nvirb1, ip, I, J
    integer(is)              :: a, b, i1, j1, p, hO1, hO2
    real(dp), parameter      :: f = sqrt2 - 1.0_dp

    call classify_orbitals(nmo1, occ1, nc1, nv1, hmap1, pmap1, io1, io2)
    nocca1 = nc1 + 2
    nvirb1 = nv1 + 2
    if (xdim1 /= nocca1*nvirb1) call mrsf_error('tdm_pairs: inconsistent xdim')
    hO1 = nc1 + 1
    hO2 = nc1 + 2

    allocate(XI(nvirb1,nocca1), XJ(nvirb1,nocca1), PP(nvirb1,nvirb1), HH(nocca1,nocca1))
    allocate(Y(nvirb1,nvirb1), W(nocca1,nocca1))

    rho = 0.0_dp

    do ip = 1, npairs
       I = ipairs(1,ip)
       J = ipairs(2,ip)
       if (I < 1 .or. I > nb .or. J < 1 .or. J > nk) &
            call mrsf_error('tdm_pairs: state index out of range')
       call expand_one(nc1, nvirb1, nocca1, mult, xb(:,I), XI)
       call expand_one(nc1, nvirb1, nocca1, mult, xk(:,J), XJ)

       ! particle-particle block: PP(a,b) = sum_k XI(a,k) XJ(b,k), with
       ! the sqrt2 factor on OO x non-OO cross terms (k in {O1,O2})
       call dgemm('N','T', nvirb1, nvirb1, nocca1, 1.0_dp, XI, nvirb1, XJ, nvirb1, 0.0_dp, PP, nvirb1)
       call dgemm('N','T', nvirb1, nvirb1, 2_is, 1.0_dp, XI(1,hO1), nvirb1, XJ(1,hO1), nvirb1, 0.0_dp, Y, nvirb1)
       PP(1:2,3:nvirb1) = PP(1:2,3:nvirb1) + f * Y(1:2,3:nvirb1)
       PP(3:nvirb1,1:2) = PP(3:nvirb1,1:2) + f * Y(3:nvirb1,1:2)

       ! hole-hole block: HH(i,j) = -sum_k XJ(k,i) XI(k,j), cross terms k in {O1,O2}
       call dgemm('T','N', nocca1, nocca1, nvirb1, -1.0_dp, XJ, nvirb1, XI, nvirb1, 0.0_dp, HH, nocca1)
       call dgemm('T','N', nocca1, nocca1, 2_is, 1.0_dp, XJ, nvirb1, XI, nvirb1, 0.0_dp, W, nocca1)
       HH(1:nc1,hO1:hO2) = HH(1:nc1,hO1:hO2) - f * W(1:nc1,hO1:hO2)
       HH(hO1:hO2,1:nc1) = HH(hO1:hO2,1:nc1) - f * W(hO1:hO2,1:nc1)

       do b = 1, nvirb1
          do a = 1, nvirb1
             rho(pmap1(a),pmap1(b),ip) = rho(pmap1(a),pmap1(b),ip) + PP(a,b)
          enddo
       enddo
       do j1 = 1, nocca1
          do i1 = 1, nocca1
             rho(hmap1(i1),hmap1(j1),ip) = rho(hmap1(i1),hmap1(j1),ip) + HH(i1,j1)
          enddo
       enddo

       if (add_ref) then
          do p = 1, nmo1
             rho(p,p,ip) = rho(p,p,ip) + occ1(p)
          enddo
       endif
    enddo

    deallocate(XI, XJ, PP, HH, Y, W, hmap1, pmap1)

  end subroutine tdm_pairs

!######################################################################
! expand_one: compressed vector -> expanded X(nvirb,nocca) (local dims)
!######################################################################
  subroutine expand_one(nc1, nvirb1, nocca1, mult, xin, Xe)

    integer(is), intent(in) :: nc1, nvirb1, nocca1, mult
    real(dp), intent(in)    :: xin(nvirb1,nocca1)
    real(dp), intent(out)   :: Xe(nvirb1,nocca1)
    real(dp)                :: xoo

    Xe = xin
    xoo = xin(1,nc1+1)
    Xe(1,nc1+1) = xoo * isqrt2
    if (mult == 1) then
       Xe(2,nc1+2) = -xoo * isqrt2
    else
       Xe(2,nc1+2) = xoo * isqrt2
       Xe(1,nc1+2) = 0.0_dp
       Xe(2,nc1+1) = 0.0_dp
    endif

  end subroutine expand_one

end module mrsf_density
