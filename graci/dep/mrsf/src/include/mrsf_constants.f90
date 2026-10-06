!**********************************************************************
! mrsf_constants: kind parameters and numerical constants
!**********************************************************************
module mrsf_constants

  implicit none

  integer, parameter :: is = selected_int_kind(8)
  integer, parameter :: ib = selected_int_kind(18)
  integer, parameter :: sp = selected_real_kind(6)
  integer, parameter :: dp = selected_real_kind(15)

  real(dp), parameter :: eh2ev  = 27.211386245988_dp
  real(dp), parameter :: sqrt2  = 1.4142135623730951_dp
  real(dp), parameter :: isqrt2 = 0.7071067811865476_dp

end module mrsf_constants
