!**********************************************************************
! mrsf_density: unrelaxed state and transition density matrices from
! MRSF and extended-MRSF amplitudes (spin-summed, MO basis,
! rho(p,q) = <bra|E_pq|ket>). Stateless: the orbital classes are rebuilt
! from the occupation vector.
!
! A vector is the (nvirb x ncol) matrix [x | y]: x the compressed MRSF
! amplitudes (nocca columns), y the core-to-virtual amplitudes of the
! extended method (nC columns, SOMO rows zero; absent when ncol = nocca).
! With the expanded amplitudes X~ (OO slots unfolded), for bra I and ket J:
!   rho = (x_I.x_J) D_ref + (y_I.y_J) D_G
!       + particle block  X~_I X~_J^T  (MRSF and CV columns in one product;
!                          sqrt2 on the OO x non-OO cross terms)
!       - hole block      [X~_J^T X~_I] on the MRSF columns (sqrt2 terms)
!       - core block      [X~_J^T X~_I] on the CV columns
!       + MRSF <-> CV cross terms: G slot (l,d), CO1 (O2,d), O2V (l,O1)
!         and their transposes
! D_ref = diag(2 C, 1 O1, 1 O2, 0 V), D_G = diag(2 C, 2 O1, 0 O2, 0 V).
! The cross-term coefficients were pinned against the determinant brute
! force (~/calculations/mrsf_dev/emrsf_ref.py, COEF_TDM).
!
! Efficiency: every state is expanded once; per pair one dgemm writes
! the whole particle block (in place for the usual MO order C,O1,O2,V)
! and one ncol x ncol Gram dgemm gives the hole blocks, the norm weights
! and the O2V cross terms; each output element is written once (plus
! O(nmo) accumulations). OpenMP over pairs (disjoint contiguous output
! slices, pairs sorted by bra state) with sequential BLAS inside.
!**********************************************************************
module mrsf_density

  use iso_c_binding, only: c_int
  use mrsf_constants
  use mrsf_io
  use mrsf_space, only: classify_orbitals
!$ use omp_lib

  implicit none

  real(dp), parameter, private :: f = sqrt2 - 1.0_dp

#ifdef USE_MKL
  interface
     function mkl_threads_local(nt) bind(c, name='MKL_Set_Num_Threads_Local') result(prev)
       import :: c_int
       integer(c_int), value :: nt
       integer(c_int)        :: prev
     end function mkl_threads_local
  end interface
#endif

contains

!######################################################################
! tdm_pairs: rho(:,:,ip) = <bra I_ip|E_pq|ket J_ip>, reference terms
! included (state densities for I = J with normalised vectors)
!######################################################################
  subroutine tdm_pairs(nmo1, occ1, mult, ncol1, xdim1, npairs, nb, nk, ipairs, xb, xk, rho)

    integer(is), intent(in) :: nmo1, mult, ncol1, xdim1, npairs, nb, nk
    real(dp), intent(in)    :: occ1(nmo1)
    integer(is), intent(in) :: ipairs(2,npairs)
    real(dp), intent(in)    :: xb(xdim1,nb), xk(xdim1,nk)
    real(dp), intent(out)   :: rho(nmo1,nmo1,npairs)

    integer(is), allocatable :: hmap1(:), pmap1(:), perm(:), off(:)
    real(dp), allocatable    :: Xbe(:,:,:), Xke(:,:,:), occG(:), Gm(:,:), PP(:)
    integer(is)              :: nc1, nv1, io1, io2, nocca1, nvirb1, ip, k, I, J, nthr
    integer(c_int)           :: prev
    logical                  :: ext, fast, par
    real(dp)                 :: cG, cCO1, cO2V

    call classify_orbitals(nmo1, occ1, nc1, nv1, hmap1, pmap1, io1, io2)
    nocca1 = nc1 + 2
    nvirb1 = nv1 + 2
    if (ncol1 /= nocca1 .and. ncol1 /= nocca1 + nc1) call mrsf_error('tdm_pairs: inconsistent ncol')
    if (xdim1 /= nvirb1*ncol1) call mrsf_error('tdm_pairs: inconsistent xdim')
    ext = ncol1 > nocca1

    ! MRSF <-> CV cross-term coefficients
    if (mult == 1) then
       cG = sqrt2;  cCO1 = -1.0_dp; cO2V = -1.0_dp
    else
       cG = 0.0_dp; cCO1 = -1.0_dp; cO2V = 1.0_dp
    endif

    ! usual MO order (C, O1, O2, V contiguous): particle block in place
    fast = .true.
    do k = 1, nocca1
       if (hmap1(k) /= k) fast = .false.
    enddo
    do k = 1, nvirb1
       if (pmap1(k) /= nc1 + k) fast = .false.
    enddo

    ! occupations of the closed-shell configuration G
    allocate(occG(nmo1))
    occG = 0.0_dp
    do k = 1, nc1
       occG(hmap1(k)) = 2.0_dp
    enddo
    occG(io1) = 2.0_dp

    ! every state expanded once
    allocate(Xbe(nvirb1,ncol1,nb), Xke(nvirb1,ncol1,nk))
    do I = 1, nb
       call expand_ext(nc1, nvirb1, nocca1, ncol1, mult, xb(1,I), Xbe(1,1,I))
    enddo
    do J = 1, nk
       call expand_ext(nc1, nvirb1, nocca1, ncol1, mult, xk(1,J), Xke(1,1,J))
    enddo

    ! pairs sorted by bra state (stable counting sort): consecutive pairs
    ! of a thread reuse the bra expansion from cache
    allocate(off(nb+1), perm(npairs))
    off = 0
    do ip = 1, npairs
       I = ipairs(1,ip)
       J = ipairs(2,ip)
       if (I < 1 .or. I > nb .or. J < 1 .or. J > nk) &
            call mrsf_error('tdm_pairs: state index out of range')
       off(I+1) = off(I+1) + 1
    enddo
    do I = 1, nb
       off(I+1) = off(I+1) + off(I)
    enddo
    do ip = 1, npairs
       I = ipairs(1,ip)
       off(I) = off(I) + 1
       perm(off(I)) = ip
    enddo

    ! OpenMP over pairs with sequential BLAS inside; a single pair (or
    ! fewer pairs than threads) uses threaded BLAS instead
    nthr = 1
!$  nthr = omp_get_max_threads()
    par = (nthr > 1) .and. (npairs >= nthr)

    !$omp parallel if(par) default(shared) private(k, ip, I, J, Gm, PP, prev)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(1_c_int)
#endif
    allocate(Gm(ncol1,ncol1))
    if (fast) then
       allocate(PP(1))
    else
       allocate(PP(nvirb1*nvirb1))
    endif
    !$omp do schedule(static)
    do k = 1, npairs
       ip = perm(k)
       I  = ipairs(1,ip)
       J  = ipairs(2,ip)
       call pair_tdm(nmo1, nc1, nv1, nocca1, nvirb1, ncol1, ext, fast, io1, hmap1, pmap1, &
            occ1, occG, cG, cCO1, cO2V, Xbe(1,1,I), Xke(1,1,J), Gm, PP, rho(1,1,ip))
    enddo
    !$omp end do
    deallocate(Gm, PP)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(0_c_int)
#endif
    !$omp end parallel

    deallocate(hmap1, pmap1, perm, off, Xbe, Xke, occG)

  end subroutine tdm_pairs

!######################################################################
! pair_tdm: rho_ip = <I|E_pq|J> for one pair of expanded vectors
!######################################################################
  subroutine pair_tdm(nmo1, nc1, nv1, nocca1, nvirb1, ncol1, ext, fast, io1, hmap1, pmap1, &
       occ1, occG, cG, cCO1, cO2V, XI, XJ, Gm, PP, rho_ip)

    integer(is), intent(in) :: nmo1, nc1, nv1, nocca1, nvirb1, ncol1, io1
    logical, intent(in)     :: ext, fast
    integer(is), intent(in) :: hmap1(nocca1), pmap1(nvirb1)
    real(dp), intent(in)    :: occ1(nmo1), occG(nmo1), cG, cCO1, cO2V
    real(dp), intent(in)    :: XI(nvirb1,ncol1), XJ(nvirb1,ncol1)
    real(dp), intent(inout) :: Gm(ncol1,ncol1), PP(*)
    real(dp), intent(inout) :: rho_ip(nmo1,nmo1)

    integer(is) :: i, j, l, d, p, hO1, hO2
    real(dp)    :: xx, yy, sI, sJ

    hO1 = nc1 + 1
    hO2 = nc1 + 2

    ! particle block, written once: in place for the usual MO order (the
    ! rest of the slice is zeroed first), else through a scratch block
    if (fast) then
       rho_ip(1:nc1,:) = 0.0_dp
       rho_ip(nc1+1:nmo1,1:nc1) = 0.0_dp
       call pblock(nv1, nocca1, nvirb1, ncol1, ext, cCO1, XI, XJ, rho_ip(nc1+1,nc1+1), nmo1)
    else
       rho_ip = 0.0_dp
       call pblock(nv1, nocca1, nvirb1, ncol1, ext, cCO1, XI, XJ, PP, nvirb1)
       call scatter_pp(nmo1, nvirb1, pmap1, PP, rho_ip)
    endif

    ! Gram matrix Gm(r,s) = sum_a XJ(a,r) XI(a,s)
    call dgemm('T','N', ncol1, ncol1, nvirb1, 1.0_dp, XJ, nvirb1, XI, nvirb1, 0.0_dp, Gm, ncol1)
    xx = 0.0_dp
    do i = 1, nocca1
       xx = xx + Gm(i,i)
    enddo
    yy = 0.0_dp
    do l = nocca1+1, ncol1
       yy = yy + Gm(l,l)
    enddo

    ! sqrt2 terms of the MRSF hole block (C x O and O x C, rows 1:2 only)
    if (nc1 > 0) then
       call dgemm('T','N', nc1, 2, 2, f, XJ(1,1), nvirb1, XI(1,hO1), nvirb1, 1.0_dp, Gm(1,hO1), ncol1)
       call dgemm('T','N', 2, nc1, 2, f, XJ(1,hO1), nvirb1, XI(1,1), nvirb1, 1.0_dp, Gm(hO1,1), ncol1)
    endif

    ! MRSF hole block
    do j = 1, nocca1
       do i = 1, nocca1
          rho_ip(hmap1(i),hmap1(j)) = rho_ip(hmap1(i),hmap1(j)) - Gm(i,j)
       enddo
    enddo

    if (ext .and. nc1 > 0) then
       ! core block of the CV part
       do l = 1, nc1
          do j = 1, nc1
             rho_ip(hmap1(j),hmap1(l)) = rho_ip(hmap1(j),hmap1(l)) - Gm(nocca1+j,nocca1+l)
          enddo
       enddo
       ! O2V cross terms: (l, O1) and (O1, l)
       do l = 1, nc1
          rho_ip(hmap1(l),io1) = rho_ip(hmap1(l),io1) + cO2V*Gm(nocca1+l,hO2)
          rho_ip(io1,hmap1(l)) = rho_ip(io1,hmap1(l)) + cO2V*Gm(hO2,nocca1+l)
       enddo
       ! G-slot cross terms (singlets): (l, d) and (d, l)
       if (cG /= 0.0_dp .and. nv1 > 0) then
          sI = cG * XI(1,hO2)
          sJ = cG * XJ(1,hO2)
          do d = 1, nv1
             p = pmap1(2+d)
             do l = 1, nc1
                rho_ip(hmap1(l),p) = rho_ip(hmap1(l),p) + sI*XJ(2+d,nocca1+l)
             enddo
          enddo
          do l = 1, nc1
             p = hmap1(l)
             do d = 1, nv1
                rho_ip(pmap1(2+d),p) = rho_ip(pmap1(2+d),p) + sJ*XI(2+d,nocca1+l)
             enddo
          enddo
       endif
    endif

    ! norm-weighted reference occupations
    do p = 1, nmo1
       rho_ip(p,p) = rho_ip(p,p) + xx*occ1(p) + yy*occG(p)
    enddo

  end subroutine pair_tdm

!######################################################################
! pblock: particle block T(a,b) of one pair (leading dimension ldT)
!######################################################################
  subroutine pblock(nv1, nocca1, nvirb1, ncol1, ext, cCO1, XI, XJ, T, ldT)

    integer(is), intent(in) :: nv1, nocca1, nvirb1, ncol1, ldT
    logical, intent(in)     :: ext
    real(dp), intent(in)    :: cCO1
    real(dp), intent(in)    :: XI(nvirb1,ncol1), XJ(nvirb1,ncol1)
    real(dp), intent(inout) :: T(ldT,*)

    integer(is) :: nc1, hO1

    nc1 = nocca1 - 2
    hO1 = nc1 + 1

    ! all columns (MRSF and CV) in one product
    call dgemm('N','T', nvirb1, nvirb1, ncol1, 1.0_dp, XI, nvirb1, XJ, nvirb1, 0.0_dp, T, ldT)

    if (nv1 > 0) then
       ! sqrt2 on the OO x non-OO cross terms (MRSF hole columns O1, O2)
       call dgemm('N','T', 2, nv1, 2, f, XI(1,hO1), nvirb1, XJ(3,hO1), nvirb1, 1.0_dp, T(1,3), ldT)
       call dgemm('N','T', nv1, 2, 2, f, XI(3,hO1), nvirb1, XJ(1,hO1), nvirb1, 1.0_dp, T(3,1), ldT)
       ! CO1 cross terms: (O2, d) and (d, O2)
       if (ext .and. nc1 > 0) then
          call dgemv('N', nv1, nc1, cCO1, XJ(3,nocca1+1), nvirb1, XI(1,1), nvirb1, 1.0_dp, T(2,3), ldT)
          call dgemv('N', nv1, nc1, cCO1, XI(3,nocca1+1), nvirb1, XJ(1,1), nvirb1, 1.0_dp, T(3,2), 1)
       endif
    endif

  end subroutine pblock

!######################################################################
! scatter_pp: rho(pmap(a), pmap(b)) += PP(a,b) (general MO order)
!######################################################################
  subroutine scatter_pp(nmo1, nvirb1, pmap1, PP, rho_ip)

    integer(is), intent(in) :: nmo1, nvirb1, pmap1(nvirb1)
    real(dp), intent(in)    :: PP(nvirb1,nvirb1)
    real(dp), intent(inout) :: rho_ip(nmo1,nmo1)
    integer(is)             :: a, b

    do b = 1, nvirb1
       do a = 1, nvirb1
          rho_ip(pmap1(a),pmap1(b)) = rho_ip(pmap1(a),pmap1(b)) + PP(a,b)
       enddo
    enddo

  end subroutine scatter_pp

!######################################################################
! expand_one: compressed MRSF vector -> expanded X(nvirb,nocca) (used by
! the overlap module)
!######################################################################
  subroutine expand_one(nc1, nvirb1, nocca1, mult, xin, Xe)

    integer(is), intent(in) :: nc1, nvirb1, nocca1, mult
    real(dp), intent(in)    :: xin(nvirb1,nocca1)
    real(dp), intent(out)   :: Xe(nvirb1,nocca1)

    call expand_ext(nc1, nvirb1, nocca1, nocca1, mult, xin, Xe)

  end subroutine expand_one

!######################################################################
! expand_ext: compressed vector -> expanded X(nvirb,ncol): OO slots
! unfolded on the MRSF columns, SOMO rows of the CV columns zeroed
!######################################################################
  subroutine expand_ext(nc1, nvirb1, nocca1, ncol1, mult, xin, Xe)

    integer(is), intent(in) :: nc1, nvirb1, nocca1, ncol1, mult
    real(dp), intent(in)    :: xin(nvirb1,ncol1)
    real(dp), intent(out)   :: Xe(nvirb1,ncol1)
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
    if (ncol1 > nocca1) Xe(1:2,nocca1+1:ncol1) = 0.0_dp

  end subroutine expand_ext

end module mrsf_density
