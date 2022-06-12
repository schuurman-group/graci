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
subroutine ref_diag_mrci_follow(irrep,nroots,confscr,nconf,vecscr) &
     bind(c,name="ref_diag_mrci_follow")
#else
subroutine ref_diag_mrci_follow(irrep,nroots,confscr,nconf,vecscr)
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

  ! Array of numbers of referecne configurations
  integer(is), intent(in)    :: nconf(0:nirrep-1)

  ! Eigenvector scratch file index
  integer(is), intent(out)   :: vecscr

  print*,''
  print*,'Here?'
  stop
  
  return
  
end subroutine ref_diag_mrci_follow
