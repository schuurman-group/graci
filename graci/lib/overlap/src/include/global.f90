!**********************************************************************
! Globally accessible variables
!**********************************************************************
module global

  use constants
  
  implicit none

  save

  ! Number of bra and ket MOs
  integer(is) :: nmoB,nmoK

  ! Number of bra and ket electrons
  integer(is) :: nelB,nelK,nel_alphaB,nel_alphaK,nel_betaB,nel_betaK
  
  ! Determinant bit strings
  integer(is)              :: n_intB,n_intK
  integer(is)              :: ndetB,ndetK
  integer(ib), allocatable :: detB(:,:,:), detK(:,:,:)

  ! Eigenvectors
  integer(is)           :: nrootsB,nrootsK
  real(dp), allocatable :: vecB(:,:),vecK(:,:)

  ! Number of unique alpha and beta strings
  integer(is)           :: nalphaB,nbetaB,nalphaK,nbetaK

  ! Unique alpha and beta strings
  integer(ib), allocatable :: alphaB(:,:),betaB(:,:),&
                              alphaK(:,:),betaK(:,:)
  
end module global
