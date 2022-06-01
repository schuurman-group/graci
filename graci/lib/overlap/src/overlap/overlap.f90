!**********************************************************************
! Public-facing interface to the overlap library
!**********************************************************************
subroutine overlap(nmoB1,nmoK1,n_intB1,n_intK1,ndetB1,ndetK1,nrootsB1,&
     nrootsK1,detB,detK,vecB,vecK)

  use constants
  use global
  
  implicit none

  ! No. MOs
  integer(is), intent(in)    :: nmoB1,nmoK1

  ! No. (n_bits)-bit integers needed to represent each alpha/beta
  ! string
  integer(is), intent(in)    :: n_intB1,n_intK1
  
  ! No. determinants
  integer(is), intent(in)    :: ndetB1,ndetK1

  ! No. roots
  integer(is), intent(in)    :: nrootsB1,nrootsK1
  
  ! Determinant bit strings
  integer(ib), intent(inout) :: detB(n_intB1,2,ndetB1)
  integer(ib), intent(inout) :: detK(n_intK1,2,ndetK1)

  ! Eigenvectors
  real(dp), intent(inout)    :: vecB(ndetB1,nrootsB1)
  real(dp), intent(inout)    :: vecK(ndetK1,nrootsK1)

  
  
  STOP
  
  return
  
end subroutine overlap
