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

  ! Conversion factors
  real(dp), parameter :: ang2bohr=1.88972612d0
  real(dp), parameter :: invcm2ev=1.23985e-4_dp
  real(dp), parameter :: eh2ev=27.2113845d0
  real(dp), parameter :: mu2me=1822.888486d0
  real(dp), parameter :: au2wn=219474.63068d0

  ! Length of the integers used to represent the dets, confs and SOPs
  integer(is), parameter :: n_bits=64
  
  ! Number of length n_bits integers needed to represent each component
  ! of a det, conf or SOP
  integer(is) :: n_int

  ! Convenient spin labels
  integer(is), parameter :: ialpha=1
  integer(is), parameter :: ibeta=2
  
end module constants
