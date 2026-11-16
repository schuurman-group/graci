!**********************************************************************
! grid_constants: kind parameters and numerical constants of libgrid
!**********************************************************************
module grid_constants

  implicit none

  integer, parameter :: is = selected_int_kind(8)
  integer, parameter :: ib = selected_int_kind(18)
  integer, parameter :: dp = selected_real_kind(15)

  real(dp), parameter :: pi = 3.14159265358979323846_dp
  ! Angstrom per bohr (PySCF's value, used for the radii tables)
  real(dp), parameter :: bohr = 0.52917721092_dp

  ! integer codes of the schemes (shared with the Python layer)
  integer(is), parameter :: rad_treutler = 1, rad_gauss_chebyshev = 2, rad_mura_knowles = 3, &
                            rad_becke = 4, rad_delley = 5
  integer(is), parameter :: prune_none = 0, prune_nwchem = 1, prune_sg1 = 2, prune_treutler = 3
  integer(is), parameter :: adjust_none = 0, adjust_treutler = 1, adjust_becke = 2
  integer(is), parameter :: step_b3 = 1, step_gmu = 2
  integer(is), parameter :: mode_becke = 1, mode_delaunay = 2

end module grid_constants
