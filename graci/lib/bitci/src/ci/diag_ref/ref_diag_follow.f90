!**********************************************************************
! Top level routines for the calculation and selection of selected
! roots of the MRCI reference space Hamiltonian based on an overlap
! criterion
!**********************************************************************

!######################################################################
! ref_diag_mrci_follow: Diagonalisation of the reference space
!                       Hamiltonian followed by the selection of a
!                       a subset of the roots based on their overlaps
!                       with a set of previously computed ones
!######################################################################
#ifdef CBINDING
subroutine ref_diag_mrci_follow(irrep,nroots,confscr,n_intR0,ndetR0,&
     nrootsR0,detR0,vecR0,nmo0,smat,nconf,vecscr) &
     bind(c,name="ref_diag_mrci_follow")
#else
subroutine ref_diag_mrci_follow(irrep,nroots,confscr,n_intR0,ndetR0,&
     nrootsR0,detR0,vecR0,nmo0,smat,nconf,vecscr)
#endif

  use constants
  use bitglobal
  use hbuild_double
  use ref_guess
  use full_diag
  use gendav
  use utils
  use iomod
  use conftype
  
  implicit none

  ! Irrep number and the requested number of roots
  integer(is), intent(in)    :: irrep
  integer(is), intent(inout) :: nroots
  
  ! Array of reference configuration scratch file numbers
  integer(is), intent(in)    :: confscr(0:nirrep-1)

  ! Eigenvectors to follow, expressed in a Slater determinant basis
  integer(is), intent(in)    :: n_intR0,ndetR0,nrootsR0
  integer(ib), intent(in)    :: detR0(n_intR0,2,ndetR0)
  real(dp), intent(in)       :: vecR0(ndetR0,nrootsR0)

  ! MO overlaps
  integer(is), intent(in)    :: nmo0
  real(dp), intent(in)       :: smat(nmo0,nmo)
  
  ! Array of numbers of reference configurations
  integer(is), intent(in)    :: nconf(0:nirrep-1)

  ! Eigenvector scratch file index
  integer(is), intent(out)   :: vecscr

  ! Everything else
  integer(is)                :: i,j


  print*,''
  do i=1,nmo
     print*,i,dot_product(smat(:,i),smat(:,i))
  enddo
  
  stop
  
  return
  
end subroutine ref_diag_mrci_follow
