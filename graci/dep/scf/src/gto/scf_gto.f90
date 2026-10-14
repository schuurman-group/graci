!**********************************************************************
! scf_gto: Gaussian basis functions and their first derivatives on grid
! points, in PySCF/libcint conventions:
!   - shell data copied from PySCF's _atm, _bas, _env (0-based pointers
!     into env; the contraction coefficients in env already carry the
!     radial normalisation),
!   - contracted radial part R_k(r) = sum_p c_kp exp(-a_p r^2) per
!     contraction k (the primitives of a generally contracted shell are
!     evaluated once),
!   - Cartesian components x^lx y^ly z^lz in libcint order (lx from l
!     down to 0, then ly from l-lx down to 0),
!   - spherical functions from PySCF's cart2sph tables (which contain the
!     s and p factors); Cartesian functions carry the s and p factors only,
!   - AOs of a shell ordered contraction-major ([k=1: functions][k=2: ...]).
! Validated against pyscf eval_gto (values and gradients, s to g,
! spherical and Cartesian) through the NumPy reference gto_ref.py.
! The screening keeps a shell on a block of points when an upper bound
! of its values and gradients there exceeds scr_cut (more permissive
! than PySCF's screen index at the same cutoff).
!**********************************************************************
module scf_gto

  use scf_constants

  implicit none

  integer(is) :: natm_b = 0, nbas_b = 0, nao_b = 0, lmax_b = 0, npmax_b = 0, ncmax_b = 0
  logical     :: cart_b = .false.
  real(dp), allocatable    :: env_b(:)
  integer(is), allocatable :: aoloc_b(:)          ! (nbas+1) 0-based AO offsets of the shells
  real(dp), allocatable    :: sh_xyz(:,:)         ! (3,nbas) shell centres
  integer(is), allocatable :: sh_l(:), sh_np(:), sh_nc(:), sh_pe(:), sh_pc(:)
  real(dp), allocatable    :: sh_tmax(:)          ! largest column sum of |angular transform|
  real(dp), allocatable    :: c2s(:)              ! cart->sph tables, column-major (ncart, nsph) per l
  integer(is)              :: c2s_off(0:15) = 0
  real(dp)                 :: scr_cut = 1.0e-15_dp

contains

  pure function ncart(l) result(n)
    integer(is), intent(in) :: l
    integer(is)             :: n
    n = (l+1)*(l+2)/2
  end function ncart

  function nfun(l) result(n)
    integer(is), intent(in) :: l
    integer(is)             :: n
    if (cart_b) then
       n = (l+1)*(l+2)/2
    else
       n = 2*l + 1
    endif
  end function nfun

  pure function cart_fac(l) result(f)
    integer(is), intent(in) :: l
    real(dp)                :: f
    f = 1.0_dp
    if (l == 0) f = fac_s
    if (l == 1) f = fac_p
  end function cart_fac

!######################################################################
! gto_set: copy the basis (PySCF _atm (6,natm), _bas (8,nbas), _env,
! ao_loc) and the cart->sph tables (c2sm(c2soff(l)+1:) = cart2sph(l),
! column-major, l = 0..lmax); ierr 1: l > 15, 2: inconsistent ao_loc
!######################################################################
  subroutine gto_set(natm, nbas, nenv, atm, bas, env, nao, aoloc, cart, lmax, nc2s, c2sm, c2soff, ierr)

    integer(is), intent(in)  :: natm, nbas, nenv, nao, lmax, nc2s
    integer(is), intent(in)  :: atm(6,natm), bas(8,nbas), aoloc(nbas+1), c2soff(0:lmax)
    real(dp), intent(in)     :: env(nenv), c2sm(nc2s)
    logical, intent(in)      :: cart
    integer(is), intent(out) :: ierr

    integer(is) :: ish, ia, l, nc, nf, m, n
    real(dp)    :: t

    ierr = 0
    call gto_free()
    if (lmax > 15) then
       ierr = 1
       return
    endif
    natm_b = natm; nbas_b = nbas; nao_b = nao; lmax_b = lmax; cart_b = cart
    allocate(env_b(nenv), aoloc_b(nbas+1), c2s(max(nc2s,1)))
    env_b = env
    aoloc_b = aoloc
    c2s(1:nc2s) = c2sm
    c2s_off(0:lmax) = c2soff(0:lmax)
    allocate(sh_xyz(3,nbas), sh_l(nbas), sh_np(nbas), sh_nc(nbas), sh_pe(nbas), sh_pc(nbas), sh_tmax(nbas))
    npmax_b = 0
    ncmax_b = 0
    n = 0
    do ish = 1, nbas
       ia = bas(1,ish) + 1
       sh_l(ish)  = bas(2,ish)
       sh_np(ish) = bas(3,ish)
       sh_nc(ish) = bas(4,ish)
       sh_pe(ish) = bas(6,ish)
       sh_pc(ish) = bas(7,ish)
       sh_xyz(:,ish) = env(atm(2,ia)+1:atm(2,ia)+3)
       npmax_b = max(npmax_b, sh_np(ish))
       ncmax_b = max(ncmax_b, sh_nc(ish))
       l  = sh_l(ish)
       nc = ncart(l)
       if (cart) then
          sh_tmax(ish) = cart_fac(l)
       else
          t = 0.0_dp
          do m = 1, 2*l+1
             t = max(t, sum(abs(c2s(c2s_off(l)+(m-1)*nc+1:c2s_off(l)+m*nc))))
          enddo
          sh_tmax(ish) = t
       endif
       nf = nfun(l)
       if (aoloc(ish+1) - aoloc(ish) /= nf*sh_nc(ish)) ierr = 2
       n = n + nf*sh_nc(ish)
    enddo
    if (n /= nao) ierr = 2

  end subroutine gto_set

  subroutine gto_free()
    if (allocated(env_b)) deallocate(env_b)
    if (allocated(aoloc_b)) deallocate(aoloc_b)
    if (allocated(c2s)) deallocate(c2s)
    if (allocated(sh_xyz)) deallocate(sh_xyz, sh_l, sh_np, sh_nc, sh_pe, sh_pc, sh_tmax)
    nbas_b = 0
    nao_b = 0
  end subroutine gto_free

!######################################################################
! gto_eval: AO values (and, for deriv = 1, the x, y, z derivatives) of
! the shells shl(1:nsh) at the np points (x, y, z); ao(ld, ncol, 1+3*deriv)
! receives the functions of the listed shells in the listed order
!######################################################################
  subroutine gto_eval(np, x, y, z, nsh, shl, deriv, ld, ncol, ao)

    integer(is), intent(in) :: np, nsh, shl(nsh), deriv, ld, ncol
    real(dp), intent(in)    :: x(np), y(np), z(np)
    real(dp), intent(out)   :: ao(ld,ncol,*)

    real(dp), allocatable :: dx(:), dy(:), dz(:), r2(:), e(:), rad(:,:), drad(:,:)
    real(dp), allocatable :: xp(:,:), yp(:,:), zp(:,:), mon(:,:), dmon(:,:,:), cv(:,:)
    integer(is) :: s, ish, l, npr, nct, nc, nf, col, p, k, g, ic, lx, ly, lz, m, comp, ncm
    real(dp)    :: a, c, t

    ncm = ncart(lmax_b)
    allocate(dx(np), dy(np), dz(np), r2(np), e(np), rad(np,ncmax_b), drad(np,ncmax_b))
    allocate(xp(np,0:lmax_b), yp(np,0:lmax_b), zp(np,0:lmax_b), mon(np,ncm), dmon(np,ncm,3), cv(np,ncm))

    col = 0
    do s = 1, nsh
       ish = shl(s)
       l   = sh_l(ish)
       npr = sh_np(ish)
       nct = sh_nc(ish)
       nc  = ncart(l)
       nf  = nfun(l)
       do g = 1, np
          dx(g) = x(g) - sh_xyz(1,ish)
          dy(g) = y(g) - sh_xyz(2,ish)
          dz(g) = z(g) - sh_xyz(3,ish)
          r2(g) = dx(g)*dx(g) + dy(g)*dy(g) + dz(g)*dz(g)
       enddo
       ! contracted radial parts and their r-derivative factors
       rad(:,1:nct) = 0.0_dp
       if (deriv > 0) drad(:,1:nct) = 0.0_dp
       do p = 1, npr
          a = env_b(sh_pe(ish)+p)
          do g = 1, np
             e(g) = exp(-min(a*r2(g), 700.0_dp))
          enddo
          do k = 1, nct
             c = env_b(sh_pc(ish)+(k-1)*npr+p)
             do g = 1, np
                rad(g,k) = rad(g,k) + c*e(g)
             enddo
             if (deriv > 0) then
                t = -2.0_dp*a*c
                do g = 1, np
                   drad(g,k) = drad(g,k) + t*e(g)
                enddo
             endif
          enddo
       enddo
       ! Cartesian monomials and their derivatives (libcint order)
       do g = 1, np
          xp(g,0) = 1.0_dp
          yp(g,0) = 1.0_dp
          zp(g,0) = 1.0_dp
       enddo
       do m = 1, l
          do g = 1, np
             xp(g,m) = xp(g,m-1)*dx(g)
             yp(g,m) = yp(g,m-1)*dy(g)
             zp(g,m) = zp(g,m-1)*dz(g)
          enddo
       enddo
       ic = 0
       do lx = l, 0, -1
          do ly = l-lx, 0, -1
             lz = l - lx - ly
             ic = ic + 1
             do g = 1, np
                mon(g,ic) = xp(g,lx)*yp(g,ly)*zp(g,lz)
             enddo
             if (deriv > 0) then
                if (lx > 0) then
                   do g = 1, np
                      dmon(g,ic,1) = real(lx,dp)*xp(g,lx-1)*yp(g,ly)*zp(g,lz)
                   enddo
                else
                   dmon(:,ic,1) = 0.0_dp
                endif
                if (ly > 0) then
                   do g = 1, np
                      dmon(g,ic,2) = real(ly,dp)*xp(g,lx)*yp(g,ly-1)*zp(g,lz)
                   enddo
                else
                   dmon(:,ic,2) = 0.0_dp
                endif
                if (lz > 0) then
                   do g = 1, np
                      dmon(g,ic,3) = real(lz,dp)*xp(g,lx)*yp(g,ly)*zp(g,lz-1)
                   enddo
                else
                   dmon(:,ic,3) = 0.0_dp
                endif
             endif
          enddo
       enddo
       ! per contraction: Cartesian values, then the angular transform
       do k = 1, nct
          do ic = 1, nc
             do g = 1, np
                cv(g,ic) = mon(g,ic)*rad(g,k)
             enddo
          enddo
          call ang_transform(np, l, nc, nf, cv, ld, ao(1,col+(k-1)*nf+1,1))
          if (deriv > 0) then
             do comp = 1, 3
                do ic = 1, nc
                   if (comp == 1) then
                      do g = 1, np
                         cv(g,ic) = dmon(g,ic,1)*rad(g,k) + mon(g,ic)*dx(g)*drad(g,k)
                      enddo
                   else if (comp == 2) then
                      do g = 1, np
                         cv(g,ic) = dmon(g,ic,2)*rad(g,k) + mon(g,ic)*dy(g)*drad(g,k)
                      enddo
                   else
                      do g = 1, np
                         cv(g,ic) = dmon(g,ic,3)*rad(g,k) + mon(g,ic)*dz(g)*drad(g,k)
                      enddo
                   endif
                enddo
                call ang_transform(np, l, nc, nf, cv, ld, ao(1,col+(k-1)*nf+1,1+comp))
             enddo
          endif
       enddo
       col = col + nct*nf
    enddo

    deallocate(dx, dy, dz, r2, e, rad, drad, xp, yp, zp, mon, dmon, cv)

  contains

    ! out(1:np, 1:nf) = cv(1:np, 1:nc) T, T = cart2sph(l) or the s/p factor
    subroutine ang_transform(np1, l1, nc1, nf1, cv1, ld1, out)
      integer(is), intent(in) :: np1, l1, nc1, nf1, ld1
      real(dp), intent(in)    :: cv1(np1,*)
      real(dp), intent(inout) :: out(ld1,*)
      integer(is) :: mm, cc, gg, off
      real(dp)    :: f
      if (cart_b) then
         f = cart_fac(l1)
         do cc = 1, nc1
            do gg = 1, np1
               out(gg,cc) = f*cv1(gg,cc)
            enddo
         enddo
      else
         off = c2s_off(l1)
         do mm = 1, nf1
            do gg = 1, np1
               out(gg,mm) = 0.0_dp
            enddo
            do cc = 1, nc1
               f = c2s(off+(mm-1)*nc1+cc)
               if (f /= 0.0_dp) then
                  do gg = 1, np1
                     out(gg,mm) = out(gg,mm) + f*cv1(gg,cc)
                  enddo
               endif
            enddo
         enddo
      endif
    end subroutine ang_transform

  end subroutine gto_eval

!######################################################################
! gto_screen: shells whose values or gradients can exceed scr_cut on
! the np points: bound = tmax(l) sum_p max_k|c_kp| [h_l + l h_{l-1} +
! 2 a_p h_{l+1}], h_k = max_{r >= rmin} r^k exp(-a_p r^2), rmin = the
! smallest distance of the points from the shell centre
!######################################################################
  subroutine gto_screen(np, x, y, z, keep)

    integer(is), intent(in) :: np
    real(dp), intent(in)    :: x(np), y(np), z(np)
    logical, intent(out)    :: keep(nbas_b)

    integer(is) :: ish, g, p, k, l, npr, nct
    real(dp)    :: rmin2, rmin, d2, a, cm, b

    do ish = 1, nbas_b
       rmin2 = huge(1.0_dp)
       do g = 1, np
          d2 = (x(g)-sh_xyz(1,ish))**2 + (y(g)-sh_xyz(2,ish))**2 + (z(g)-sh_xyz(3,ish))**2
          rmin2 = min(rmin2, d2)
       enddo
       rmin = sqrt(rmin2)
       l   = sh_l(ish)
       npr = sh_np(ish)
       nct = sh_nc(ish)
       b = 0.0_dp
       do p = 1, npr
          a  = env_b(sh_pe(ish)+p)
          cm = 0.0_dp
          do k = 1, nct
             cm = max(cm, abs(env_b(sh_pc(ish)+(k-1)*npr+p)))
          enddo
          b = b + cm*(hmax(l, a, rmin) + 2.0_dp*a*hmax(l+1, a, rmin))
          if (l > 0) b = b + cm*real(l,dp)*hmax(l-1, a, rmin)
       enddo
       keep(ish) = b*sh_tmax(ish) > scr_cut
    enddo

  contains

    pure function hmax(kk, aa, r0) result(h)
      integer(is), intent(in) :: kk
      real(dp), intent(in)    :: aa, r0
      real(dp)                :: h, rs
      if (kk == 0) then
         h = exp(-min(aa*r0*r0, 700.0_dp))
      else
         rs = sqrt(real(kk,dp)/(2.0_dp*aa))
         if (r0 <= rs) then
            h = rs**kk*exp(-0.5_dp*real(kk,dp))
         else
            h = r0**kk*exp(-min(aa*r0*r0, 700.0_dp))
         endif
      endif
    end function hmax

  end subroutine gto_screen

end module scf_gto
