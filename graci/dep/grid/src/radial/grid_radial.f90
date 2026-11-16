!**********************************************************************
! grid_radial: radial quadratures r(n), dr(n) (PySCF's dft.radi
! conventions and point order; the radial weight is 4 pi r^2 dr)
!   rad_treutler         Treutler-Ahlrichs M4 with the per-element xi
!   rad_gauss_chebyshev  Krack-Koster Gauss-Chebyshev (second kind)
!   rad_mura_knowles     Mura-Knowles log3
!   rad_becke            Becke's Chebyshev mapping with Bragg radii
!   rad_delley           Delley log2
!**********************************************************************
module grid_radial

  use grid_constants
  use grid_tables

  implicit none

contains

  subroutine radial_grid(scheme, n, z, r, dr, ierr)

    integer(is), intent(in)  :: scheme, n, z
    real(dp), intent(out)    :: r(n), dr(n)
    integer(is), intent(out) :: ierr
    real(dp), allocatable    :: x1(:), xi(:)
    real(dp)                 :: xif, step, ln2, x, fac, far, rm, t, w, r_outer, rfac
    integer(is)              :: i

    ierr = 0
    if (n < 1) then
       ierr = 1
       return
    endif

    select case (scheme)

    case (rad_treutler)
       xif = 1.0_dp
       if (z >= 0 .and. z <= zmax_xi) xif = treutler_xi(z+1)
       step = pi/real(n+1, dp)
       ln2 = xif/log(2.0_dp)
       do i = 1, n
          x = cos(real(i, dp)*step)
          ! PySCF fills index i-1 and reverses: r(n+1-i)
          r(n+1-i) = -ln2*(1.0_dp + x)**0.6_dp*log((1.0_dp - x)/2.0_dp)
          dr(n+1-i) = step*sin(real(i, dp)*step)*ln2*(1.0_dp + x)**0.6_dp &
               *(-0.6_dp/(1.0_dp + x)*log((1.0_dp - x)/2.0_dp) + 1.0_dp/(1.0_dp - x))
       enddo

    case (rad_gauss_chebyshev)
       allocate(x1(n), xi(n))
       ln2 = 1.0_dp/log(2.0_dp)
       fac = 16.0_dp/3.0_dp/real(n+1, dp)
       do i = 1, n
          x1(i) = real(i, dp)*pi/real(n+1, dp)
          xi(i) = real(n - 1 - (i-1)*2, dp)/real(n+1, dp) &
               + (1.0_dp + 2.0_dp/3.0_dp*sin(x1(i))**2)*sin(2.0_dp*x1(i))/pi
       enddo
       do i = 1, n
          r(i) = (xi(i) - xi(n+1-i))/2.0_dp
       enddo
       xi = r
       do i = 1, n
          r(i) = 1.0_dp - log(1.0_dp + xi(i))*ln2
          dr(i) = fac*sin(x1(i))**4*ln2/(1.0_dp + xi(i))
       enddo
       deallocate(x1, xi)

    case (rad_mura_knowles)
       far = 5.2_dp
       if (z == 3 .or. z == 4 .or. z == 11 .or. z == 12 .or. z == 19 .or. z == 20) far = 7.0_dp
       do i = 1, n
          x = (real(i-1, dp) + 0.5_dp)/real(n, dp)
          r(i) = -far*log(1.0_dp - x**3)
          dr(i) = far*3.0_dp*x*x/((1.0_dp - x**3)*real(n, dp))
       enddo

    case (rad_becke)
       if (z < 0 .or. z > zmax_bragg) then
          ierr = 2
          return
       endif
       if (z == 1) then
          rm = bragg_radii(z+1)
       else
          rm = bragg_radii(z+1)*0.5_dp
       endif
       do i = 1, n
          t = cos(real(i, dp)*pi/real(n+1, dp))
          w = pi/real(n+1, dp)*sin(real(i, dp)*pi/real(n+1, dp))
          r(i) = (1.0_dp + t)/(1.0_dp - t)*rm
          dr(i) = w*2.0_dp/(1.0_dp - t)**2*rm
       enddo

    case (rad_delley)
       r_outer = 12.0_dp
       step = 1.0_dp/real(n+1, dp)
       rfac = r_outer/log(1.0_dp - (real(n, dp)*step)**2)
       do i = 1, n
          x = real(i, dp)*step
          r(i) = rfac*log(1.0_dp - x**2)
          dr(i) = rfac*(-2.0_dp*real(i, dp)*step**2)/(1.0_dp - x**2)
       enddo

    case default
       ierr = 3
    end select

  end subroutine radial_grid

end module grid_radial
