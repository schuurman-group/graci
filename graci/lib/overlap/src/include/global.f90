!**********************************************************************
! Globally accessible variables
!**********************************************************************
module global

  use constants
  
  implicit none

  save

  ! Number of bra and ket MOs
  integer(is) :: nmoB,nmoK
  
end module global
