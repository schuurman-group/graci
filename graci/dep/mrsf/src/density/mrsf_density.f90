!**********************************************************************
! mrsf_density: unrelaxed state and transition density matrices from
! MRSF and extended-MRSF amplitudes (spin-summed, MO basis,
! rho(p,q) = <bra|E_pq|ket>). Stateless: the orbital classes are rebuilt
! from the occupation vector; the parent of the extension enters through
! its weights.
!
! A vector is the (nvirb x ncol) matrix [x | y]: x the compressed MRSF
! amplitudes (nocca columns), y the core-to-virtual amplitudes of the
! extended method (nC columns, SOMO rows zero; absent when ncol = nocca).
! The C->V configurations are the singles of the parent
!   Pi = w_G G + w_D D + w_S S   (G = |C^2 O1^2|, D = |C^2 O2^2|,
!   S = the open-shell singlet SOMO pair; the published parent is
!   w = (1, 0, 0)), weights wpi(3).
! With the expanded amplitudes X~ (OO slots unfolded), for bra I and ket J:
!   rho = (x_I.x_J) D_ref + (y_I.y_J) D_Pi
!       + particle block  X~_I X~_J^T  (MRSF and CV columns in one product;
!                          sqrt2 on the OO x non-OO cross terms)
!       - hole block      [X~_J^T X~_I] on the MRSF columns (sqrt2 terms)
!       - core block      [X~_J^T X~_I] on the CV columns
!       + MRSF <-> CV cross terms, sum over the parents X of w_X R_X:
!           slot terms   (l,d): G slot (parent G), D slot (D), OO slot (S)
!           CO slots     (O_r,d): CO1 -> O2 (G), CO2 -> O1 (D),
!                                 CO1 -> O1 and CO2 -> O2 (S)
!           OV slots     (l,O_r): O2V -> O1 (G), O1V -> O2 (D),
!                                 O1V -> O1 and O2V -> O2 (S)
!         and their transposes (coefficients in density_tables)
! D_ref = diag(2 C, 1 O1, 1 O2, 0 V); D_Pi = 2 P_C + the SOMO block
!   [[2wG^2+wS^2, sqrt2 wS (wD-wG)], [., 2wD^2+wS^2]]  (= D_G for w = e_G).
! The cross-term coefficients were pinned against the determinant brute
! force (~/calculations/mrsf_dev/emrsf_ref.py COEF_TDM for the published
! parent; inv/emrsf_pi_density_ref.py, inv/coef_pi_density.py for the
! three parents). Unnormalised weights are handled exactly: the CV
! columns are scaled by |w| and the tables use w/|w|.
!
! Efficiency: every state is expanded once; per pair one dgemm writes
! the whole particle block (in place for the usual MO order C,O1,O2,V)
! and one ncol x ncol Gram dgemm gives the hole blocks, the norm weights
! and the OV cross terms; the cross terms are O(nV nC) dgemv/rank-1
! updates, only the ones with a nonzero contracted coefficient are
! issued (the published parent executes exactly the Part VII calls);
! each output element is written once (plus O(nmo) accumulations).
! OpenMP over pairs (disjoint contiguous output slices, pairs sorted by
! bra state) with sequential BLAS inside.
!**********************************************************************
module mrsf_density_mod

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
! included (state densities for I = J with normalised vectors);
! wpi(3) = parent weights (w_G, w_D, w_S) of the extension ((1,0,0) for
! the published parent and for standard MRSF vectors)
!######################################################################
  subroutine tdm_pairs(nmo1, occ1, mult, ncol1, xdim1, npairs, nb, nk, wpi, ipairs, xb, xk, rho)

    integer(is), intent(in) :: nmo1, mult, ncol1, xdim1, npairs, nb, nk
    real(dp), intent(in)    :: occ1(nmo1), wpi(3)
    integer(is), intent(in) :: ipairs(2,npairs)
    real(dp), intent(in)    :: xb(xdim1,nb), xk(xdim1,nk)
    real(dp), intent(out)   :: rho(nmo1,nmo1,npairs)

    integer(is), allocatable :: hmap1(:), pmap1(:), perm(:), off(:)
    real(dp), allocatable    :: Xbe(:,:,:), Xke(:,:,:), dpi(:), Gm(:,:), PP(:)
    integer(is)              :: nc1, nv1, io1, io2, nocca1, nvirb1, ip, k, I, J, nthr
    integer(c_int)           :: prev
    logical                  :: ext, fast, par
    real(dp)                 :: w(3), wn, ysc, dpi12, cslot(3), cCO(2,2), cOV(2,2)

    call classify_orbitals(nmo1, occ1, nc1, nv1, hmap1, pmap1, io1, io2)
    nocca1 = nc1 + 2
    nvirb1 = nv1 + 2
    if (ncol1 /= nocca1 .and. ncol1 /= nocca1 + nc1) call mrsf_error('tdm_pairs: inconsistent ncol')
    if (xdim1 /= nvirb1*ncol1) call mrsf_error('tdm_pairs: inconsistent xdim')
    ext = ncol1 > nocca1

    ! parent weights: normalised tables, the norm on the CV columns
    wn = sqrt(wpi(1)**2 + wpi(2)**2 + wpi(3)**2)
    if (wn > 0.0_dp) then
       w = wpi / wn
       ysc = wn
    else
       w = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
       ysc = 1.0_dp
    endif
    call density_coefs(mult, w, cslot, cCO, cOV)

    ! usual MO order (C, O1, O2, V contiguous): particle block in place
    fast = .true.
    do k = 1, nocca1
       if (hmap1(k) /= k) fast = .false.
    enddo
    do k = 1, nvirb1
       if (pmap1(k) /= nc1 + k) fast = .false.
    enddo

    ! reference density of the parent Pi: diagonal plus the SOMO off-diagonal
    allocate(dpi(nmo1))
    dpi = 0.0_dp
    do k = 1, nc1
       dpi(hmap1(k)) = 2.0_dp
    enddo
    dpi(io1) = 2.0_dp*w(1)*w(1) + w(3)*w(3)
    dpi(io2) = 2.0_dp*w(2)*w(2) + w(3)*w(3)
    dpi12 = sqrt2*w(3)*(w(2) - w(1))

    ! every state expanded once
    allocate(Xbe(nvirb1,ncol1,nb), Xke(nvirb1,ncol1,nk))
    do I = 1, nb
       call expand_ext(nc1, nvirb1, nocca1, ncol1, mult, ysc, xb(1,I), Xbe(1,1,I))
    enddo
    do J = 1, nk
       call expand_ext(nc1, nvirb1, nocca1, ncol1, mult, ysc, xk(1,J), Xke(1,1,J))
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
       call pair_tdm(nmo1, nc1, nv1, nocca1, nvirb1, ncol1, ext, fast, io1, io2, hmap1, pmap1, &
            occ1, dpi, dpi12, cslot, cCO, cOV, Xbe(1,1,I), Xke(1,1,J), Gm, PP, rho(1,1,ip))
    enddo
    !$omp end do
    deallocate(Gm, PP)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(0_c_int)
#endif
    !$omp end parallel

    deallocate(hmap1, pmap1, perm, off, Xbe, Xke, dpi)

  end subroutine tdm_pairs

!######################################################################
! density_tables: cross-term coefficients of one parent X (1 G, 2 D,
! 3 S) and multiplicity (pinned, plan Part XIX section 1):
!   tslot(3): the (l,d) term of the G slot, the D slot and the compressed
!             OO amplitude (singlets only)
!   tCO(r,s): particle row O_r from the CO slot row s (1 = CO1, 2 = CO2)
!   tOV(r,s): hole column O_r from the OV slot column s (1 = O1V, 2 = O2V)
! The entries equal the Fock-type entries of the parent coupling tables
! (mrsf_extended parent_tables), as the one-body part of the generalised
! Slater-Condon identity requires.
!######################################################################
  subroutine density_tables(mult, X, tslot, tCO, tOV)

    integer(is), intent(in) :: mult, X
    real(dp), intent(out)   :: tslot(3), tCO(2,2), tOV(2,2)
    real(dp), parameter     :: h2 = 0.7071067811865476_dp
    logical                 :: sgl

    sgl = (mult == 1)
    tslot = 0.0_dp; tCO = 0.0_dp; tOV = 0.0_dp
    select case (X)
    case (1)   ! parent G
       if (sgl) tslot(1) = sqrt2
       tCO(2,1) = -1.0_dp
       tOV(1,2) = merge(-1.0_dp, 1.0_dp, sgl)
    case (2)   ! parent D (the O1 <-> O2 mirror of G)
       if (sgl) tslot(2) = sqrt2
       tCO(1,2) = -1.0_dp
       tOV(2,1) = merge(-1.0_dp, 1.0_dp, sgl)
    case (3)   ! parent S (open-shell singlet pair)
       if (sgl) tslot(3) = sqrt2
       tCO(1,1) = -h2
       tCO(2,2) = h2
       tOV(1,1) = merge(-h2, h2, sgl)
       tOV(2,2) = merge(h2, -h2, sgl)
    end select

  end subroutine density_tables

!######################################################################
! density_coefs: the tables contracted with the (normalised) parent
! weights; zero weights contribute nothing (no 0 * t terms), so the
! published parent reproduces the Part VII coefficients exactly
!######################################################################
  subroutine density_coefs(mult, w, cslot, cCO, cOV)

    integer(is), intent(in) :: mult
    real(dp), intent(in)    :: w(3)
    real(dp), intent(out)   :: cslot(3), cCO(2,2), cOV(2,2)
    real(dp)                :: tslot(3), tCO(2,2), tOV(2,2)
    integer(is)             :: X

    cslot = 0.0_dp; cCO = 0.0_dp; cOV = 0.0_dp
    do X = 1, 3
       if (w(X) == 0.0_dp) cycle
       call density_tables(mult, X, tslot, tCO, tOV)
       cslot = cslot + w(X)*tslot
       cCO = cCO + w(X)*tCO
       cOV = cOV + w(X)*tOV
    enddo

  end subroutine density_coefs

!######################################################################
! pair_tdm: rho_ip = <I|E_pq|J> for one pair of expanded vectors
!######################################################################
  subroutine pair_tdm(nmo1, nc1, nv1, nocca1, nvirb1, ncol1, ext, fast, io1, io2, hmap1, pmap1, &
       occ1, dpi, dpi12, cslot, cCO, cOV, XI, XJ, Gm, PP, rho_ip)

    integer(is), intent(in) :: nmo1, nc1, nv1, nocca1, nvirb1, ncol1, io1, io2
    logical, intent(in)     :: ext, fast
    integer(is), intent(in) :: hmap1(nocca1), pmap1(nvirb1)
    real(dp), intent(in)    :: occ1(nmo1), dpi(nmo1), dpi12, cslot(3), cCO(2,2), cOV(2,2)
    real(dp), intent(in)    :: XI(nvirb1,ncol1), XJ(nvirb1,ncol1)
    real(dp), intent(inout) :: Gm(ncol1,ncol1), PP(*)
    real(dp), intent(inout) :: rho_ip(nmo1,nmo1)

    integer(is) :: i, j, l, d, p, r, s, hO1, hO2, ior(2), hOs(2)
    real(dp)    :: xx, yy, sI, sJ

    hO1 = nc1 + 1
    hO2 = nc1 + 2
    ior = (/ io1, io2 /)
    hOs = (/ hO1, hO2 /)

    ! particle block, written once: in place for the usual MO order (the
    ! rest of the slice is zeroed first), else through a scratch block
    if (fast) then
       rho_ip(1:nc1,:) = 0.0_dp
       rho_ip(nc1+1:nmo1,1:nc1) = 0.0_dp
       call pblock(nv1, nocca1, nvirb1, ncol1, ext, cCO, XI, XJ, rho_ip(nc1+1,nc1+1), nmo1)
    else
       rho_ip = 0.0_dp
       call pblock(nv1, nocca1, nvirb1, ncol1, ext, cCO, XI, XJ, PP, nvirb1)
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

    ! sqrt2 terms of the MRSF hole block (C x O and O x C, rows 1:2 only);
    ! the CV rows/columns nocca+l of Gm are not touched (the OV cross
    ! terms below read Gm(nocca+l, hO_s) = sum_d y_J(d,l) X~_I(2+d, hO_s))
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
       ! OV cross terms: (l, O_r) and (O_r, l) from the slot columns hO_s
       do s = 1, 2
          do r = 1, 2
             if (cOV(r,s) == 0.0_dp) cycle
             do l = 1, nc1
                rho_ip(hmap1(l),ior(r)) = rho_ip(hmap1(l),ior(r)) + cOV(r,s)*Gm(nocca1+l,hOs(s))
                rho_ip(ior(r),hmap1(l)) = rho_ip(ior(r),hmap1(l)) + cOV(r,s)*Gm(hOs(s),nocca1+l)
             enddo
          enddo
       enddo
       ! slot cross terms (singlets): (l, d) and (d, l); the OO slot
       ! enters with its compressed amplitude sqrt2 X~(1,hO1)
       if (nv1 > 0 .and. any(cslot /= 0.0_dp)) then
          sI = 0.0_dp
          sJ = 0.0_dp
          if (cslot(1) /= 0.0_dp) then
             sI = sI + cslot(1)*XI(1,hO2)
             sJ = sJ + cslot(1)*XJ(1,hO2)
          endif
          if (cslot(2) /= 0.0_dp) then
             sI = sI + cslot(2)*XI(2,hO1)
             sJ = sJ + cslot(2)*XJ(2,hO1)
          endif
          if (cslot(3) /= 0.0_dp) then
             sI = sI + cslot(3)*sqrt2*XI(1,hO1)
             sJ = sJ + cslot(3)*sqrt2*XJ(1,hO1)
          endif
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

    ! norm-weighted reference densities
    do p = 1, nmo1
       rho_ip(p,p) = rho_ip(p,p) + xx*occ1(p) + yy*dpi(p)
    enddo
    if (ext .and. dpi12 /= 0.0_dp) then
       rho_ip(io1,io2) = rho_ip(io1,io2) + yy*dpi12
       rho_ip(io2,io1) = rho_ip(io2,io1) + yy*dpi12
    endif

  end subroutine pair_tdm

!######################################################################
! pblock: particle block T(a,b) of one pair (leading dimension ldT)
!######################################################################
  subroutine pblock(nv1, nocca1, nvirb1, ncol1, ext, cCO, XI, XJ, T, ldT)

    integer(is), intent(in) :: nv1, nocca1, nvirb1, ncol1, ldT
    logical, intent(in)     :: ext
    real(dp), intent(in)    :: cCO(2,2)
    real(dp), intent(in)    :: XI(nvirb1,ncol1), XJ(nvirb1,ncol1)
    real(dp), intent(inout) :: T(ldT,*)

    integer(is) :: nc1, hO1, r, s

    nc1 = nocca1 - 2
    hO1 = nc1 + 1

    ! all columns (MRSF and CV) in one product
    call dgemm('N','T', nvirb1, nvirb1, ncol1, 1.0_dp, XI, nvirb1, XJ, nvirb1, 0.0_dp, T, ldT)

    if (nv1 > 0) then
       ! sqrt2 on the OO x non-OO cross terms (MRSF hole columns O1, O2)
       call dgemm('N','T', 2, nv1, 2, f, XI(1,hO1), nvirb1, XJ(3,hO1), nvirb1, 1.0_dp, T(1,3), ldT)
       call dgemm('N','T', nv1, 2, 2, f, XI(3,hO1), nvirb1, XJ(1,hO1), nvirb1, 1.0_dp, T(3,1), ldT)
       ! CO cross terms: (O_r, d) and (d, O_r) from the slot rows s
       if (ext .and. nc1 > 0) then
          do s = 1, 2
             do r = 1, 2
                if (cCO(r,s) == 0.0_dp) cycle
                call dgemv('N', nv1, nc1, cCO(r,s), XJ(3,nocca1+1), nvirb1, XI(s,1), nvirb1, 1.0_dp, T(r,3), ldT)
                call dgemv('N', nv1, nc1, cCO(r,s), XI(3,nocca1+1), nvirb1, XJ(s,1), nvirb1, 1.0_dp, T(3,r), 1)
             enddo
          enddo
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

    call expand_ext(nc1, nvirb1, nocca1, nocca1, mult, 1.0_dp, xin, Xe)

  end subroutine expand_one

!######################################################################
! expand_ext: compressed vector -> expanded X(nvirb,ncol): OO slots
! unfolded on the MRSF columns, SOMO rows of the CV columns zeroed, CV
! columns scaled by ysc (the norm of the parent weights)
!######################################################################
  subroutine expand_ext(nc1, nvirb1, nocca1, ncol1, mult, ysc, xin, Xe)

    integer(is), intent(in) :: nc1, nvirb1, nocca1, ncol1, mult
    real(dp), intent(in)    :: ysc
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
    if (ncol1 > nocca1) then
       Xe(1:2,nocca1+1:ncol1) = 0.0_dp
       if (ysc /= 1.0_dp) Xe(3:nvirb1,nocca1+1:ncol1) = ysc * Xe(3:nvirb1,nocca1+1:ncol1)
    endif

  end subroutine expand_ext

end module mrsf_density_mod
