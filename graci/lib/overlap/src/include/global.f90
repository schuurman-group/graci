!**********************************************************************
! Globally accessible variables
!**********************************************************************
module global

  use constants
  
  implicit none

  save

  ! Verbose output
  logical     :: verbose
  
  ! Number of bra and ket MOs
  integer(is) :: nmoB,nmoK

  ! Number of bra and ket electrons
  integer(is) :: nelB,nelK,nel_alphaB,nel_alphaK,nel_betaB,nel_betaK

  ! Bra-ket MO overlaps
  real(dp), allocatable :: smo(:,:)
  
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

  ! Alpha string offsets
  integer(is), allocatable :: offsetB(:),offsetK(:)

  ! Detetrminant-to-beta-string mapping
  integer(is), allocatable :: det2betaB(:),det2betaK(:)

  ! Unique beta factors
  real(dp), allocatable    :: betafac(:,:)

  ! Determinant thresholds
  real(dp), parameter      :: hthrsh=1e-6_dp
  real(dp), parameter      :: fthrsh=1e-6_dp
  
end module global
