module constants

  save

  ! Short integer (32-bit)
  integer, parameter :: is=selected_int_kind(8)

  ! Long integer (64-bit)
  integer, parameter :: il=selected_int_kind(16)

  ! Bit-kind integer (64-bit) [used for representing bit strings]
  integer, parameter :: ib=selected_int_kind(16)

  ! Single precision real (32-bit)
  integer, parameter :: sp=selected_real_kind(4)

  ! Double precision real (64-bit)
  integer, parameter :: dp=selected_real_kind(8)

  ! Quadruple precision real (128-bit)
  integer, parameter :: qp=selected_real_kind(24)

  ! Real zero
  real(dp), parameter :: rzero=0._dp

  ! Real one
  real(dp), parameter :: rone=1._dp

  ! Complex zero
  complex(dp), parameter :: czero=(0._dp,0._dp)

  ! Complex one
  complex(dp), parameter :: cone=(1._dp,0._dp)

  ! Complex i
  complex(dp), parameter :: ci=(0._dp,1._dp)

  ! pi
  real(dp), parameter :: pi=3.1415926535897932_dp

  ! Number of 64-bit integers needed to represent each alpha and beta
  ! string in a determinant
  integer(is) :: n_int
  
end module constants
