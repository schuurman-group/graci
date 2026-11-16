!**********************************************************************
! grid_step: smoothed step functions and the atomic-size adjustment.
! PySCF sign convention: mu_ij = (|r-R_i| - |r-R_j|)/d_ij,
! nu_ij = mu_ij + a_ij (1 - mu_ij^2), s_ij = S(-nu_ij) where S is the
! paper's step rising from 0 at -1 to 1 at +1:
!   step_b3:  s = (1 - p3(nu))/2, p(x) = (3x - x^3)/2 iterated 3 times
!   step_gmu: s = (1-nu)^mu / ((1-nu)^mu + (1+nu)^mu)
! and the cell of nucleus i is the product of its s_ij.
!**********************************************************************
module grid_step

  use grid_constants
  use grid_tables

  implicit none

contains

!######################################################################
! adjust_table: a(i,j) = (rad_j/rad_i - rad_i/rad_j)/4 clipped to
! [-1/2, 1/2] with rad = sqrt(Bragg) (Treutler) or Bragg (Becke);
! zero for adjust_none
!######################################################################
  subroutine adjust_table(natm, Z, adjust, a, ierr)

    integer(is), intent(in)  :: natm, Z(natm), adjust
    real(dp), intent(out)    :: a(natm,natm)
    integer(is), intent(out) :: ierr
    real(dp), allocatable    :: rad(:)
    integer(is)              :: i, j

    ierr = 0
    a = 0.0_dp
    if (adjust == adjust_none) return
    if (adjust /= adjust_treutler .and. adjust /= adjust_becke) then
       ierr = 41
       return
    endif
    allocate(rad(natm))
    do i = 1, natm
       if (Z(i) < 0 .or. Z(i) > zmax_bragg) then
          ierr = 42
          return
       endif
       if (adjust == adjust_treutler) then
          rad(i) = sqrt(bragg_radii(Z(i)+1)) + 1.0e-200_dp
       else
          rad(i) = bragg_radii(Z(i)+1) + 1.0e-200_dp
       endif
    enddo
    do j = 1, natm
       do i = 1, natm
          a(i,j) = 0.25_dp*(rad(j)/rad(i) - rad(i)/rad(j))
          if (a(i,j) < -0.5_dp) a(i,j) = -0.5_dp
          if (a(i,j) > 0.5_dp) a(i,j) = 0.5_dp
       enddo
    enddo
    deallocate(rad)

  end subroutine adjust_table

!######################################################################
! step_check: validity of the step code and parameter
!######################################################################
  subroutine step_check(step, mu, ierr)

    integer(is), intent(in)  :: step
    real(dp), intent(in)     :: mu
    integer(is), intent(out) :: ierr

    ierr = 0
    if (step == step_b3) return
    if (step == step_gmu) then
       if (mu < 1.0_dp) ierr = 44
       return
    endif
    ierr = 43

  end subroutine step_check

!######################################################################
! step_eval: s(nu)
!######################################################################
  pure function step_eval(step, mu, nu) result(s)

    integer(is), intent(in) :: step
    real(dp), intent(in)    :: mu, nu
    real(dp)                :: s, p, x, aa, bb

    if (step == step_b3) then
       p = nu
       p = (3.0_dp - p*p)*p*0.5_dp
       p = (3.0_dp - p*p)*p*0.5_dp
       p = (3.0_dp - p*p)*p*0.5_dp
       s = 0.5_dp*(1.0_dp - p)
    else
       x = min(max(nu, -1.0_dp), 1.0_dp)
       aa = (1.0_dp - x)**mu
       bb = (1.0_dp + x)**mu
       s = aa/(aa + bb)
    endif

  end function step_eval

!######################################################################
! step_eval_der: s(nu) and ds/dnu
!######################################################################
  pure subroutine step_eval_der(step, mu, nu, s, ds)

    integer(is), intent(in) :: step
    real(dp), intent(in)    :: mu, nu
    real(dp), intent(out)   :: s, ds
    real(dp)                :: p1, p2, p3, x, aa, bb

    if (step == step_b3) then
       p1 = (3.0_dp - nu*nu)*nu*0.5_dp
       p2 = (3.0_dp - p1*p1)*p1*0.5_dp
       p3 = (3.0_dp - p2*p2)*p2*0.5_dp
       s = 0.5_dp*(1.0_dp - p3)
       ds = -27.0_dp/16.0_dp*(1.0_dp - nu*nu)*(1.0_dp - p1*p1)*(1.0_dp - p2*p2)
    else
       x = min(max(nu, -1.0_dp), 1.0_dp)
       aa = (1.0_dp - x)**mu
       bb = (1.0_dp + x)**mu
       s = aa/(aa + bb)
       if (x <= -1.0_dp .or. x >= 1.0_dp) then
          ds = 0.0_dp
       else
          ds = -2.0_dp*mu*(1.0_dp - x)**(mu - 1.0_dp)*(1.0_dp + x)**(mu - 1.0_dp)/(aa + bb)**2
       endif
    endif

  end subroutine step_eval_der

end module grid_step
