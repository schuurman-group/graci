!**********************************************************************
! scf_constants: kind parameters and numerical constants
!**********************************************************************
module scf_constants

  implicit none

  integer, parameter :: is = selected_int_kind(8)
  integer, parameter :: ib = selected_int_kind(18)
  integer, parameter :: dp = selected_real_kind(15)

  ! factors of the s and p functions in libcint's convention
  real(dp), parameter :: fac_s = 0.282094791773878143_dp
  real(dp), parameter :: fac_p = 0.488602511902919921_dp

end module scf_constants
