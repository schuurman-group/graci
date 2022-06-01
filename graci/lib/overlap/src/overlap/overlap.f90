!**********************************************************************
! Public-facing interface to the overlap library
!**********************************************************************
subroutine overlap(nmoB1,nmoK1,n_intB1,n_intK1,ndetB1,ndetK1,nrootsB1,&
     nrootsK1,detB1,detK1,vecB1,vecK1,normthrsh)

  use constants
  use global
  use detfuncs
  use detsort
  
  implicit none

  ! No. MOs
  integer(is), intent(in) :: nmoB1,nmoK1

  ! No. (n_bits)-bit integers needed to represent each alpha/beta
  ! string
  integer(is), intent(in) :: n_intB1,n_intK1
  
  ! No. determinants
  integer(is), intent(in) :: ndetB1,ndetK1

  ! No. roots
  integer(is), intent(in) :: nrootsB1,nrootsK1
  
  ! Untruncated determinant bit strings
  integer(ib), intent(in) :: detB1(n_intB1,2,ndetB1)
  integer(ib), intent(in) :: detK1(n_intK1,2,ndetK1)

  ! Untruncated eigenvectors
  real(dp), intent(in)    :: vecB1(ndetB1,nrootsB1)
  real(dp), intent(in)    :: vecK1(ndetK1,nrootsK1)

  ! Norm-based truncation threshold
  real(dp), intent(in)    :: normthrsh

!----------------------------------------------------------------------
! Make sure that all globally accessible allocatable arrays are
! not allocated
!----------------------------------------------------------------------
  if (allocated(detB)) deallocate(detB)
  if (allocated(detK)) deallocate(detK)
  if (allocated(vecB)) deallocate(vecB)
  if (allocated(vecK)) deallocate(vecK)
  
!----------------------------------------------------------------------
! Set some globally accessible variables
!----------------------------------------------------------------------
  ! No. MOs
  nmoB=nmoB1
  nmoK=nmoK1

  ! No. integers needed to represent each alpha- and beta-string
  n_intB=n_intB1
  n_intK=n_intK1

  ! No. states
  nrootsB=nrootsB1
  nrootsK=nrootsK1

  ! No. electrons
  call get_nel(n_intB,detB1(:,:,1),nelB,nel_alphaB,nel_betaB)
  call get_nel(n_intK,detK1(:,:,1),nelK,nel_alphaK,nel_betaK)

!----------------------------------------------------------------------
! Truncate the wave functions
!----------------------------------------------------------------------
  ! Bra
  call truncate_wave_functions(n_intB,ndetB1,nrootsB,detB1,vecB1,&
       normthrsh,ndetB,detB,vecB)

  ! Ket
  call truncate_wave_functions(n_intK,ndetK1,nrootsK,detK1,vecK1,&
       normthrsh,ndetK,detK,vecK)

!----------------------------------------------------------------------
! Get the unique alpha and beta strings
!----------------------------------------------------------------------
  ! Bra
  call unique_strings(n_intB,ndetB,detB)
  
  STOP
  
  return
  
end subroutine overlap
