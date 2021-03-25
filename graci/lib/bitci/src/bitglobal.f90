module bitglobal

  use constants
  
  implicit none

  save

  !
  ! Name of the scratch directory
  !
  character(len=15) :: scratchdir
  
  !
  ! Dimensions
  !
  integer(is) :: nmo,nel,nel_alpha,nel_beta

  !
  ! Spin multiplicity
  !
  integer(is) :: imult
  
  !
  ! Determinants
  !
  integer(ib), allocatable :: det0(:,:)

  !
  ! MO information
  !
  integer(is), allocatable :: mosym(:)
  real(dp), allocatable    :: moen(:)
  
  !
  ! Symmetry
  !
  integer(is)      :: ipg             ! Index of the point group
  integer(is)      :: nirrep(8)       ! Number of irreps in the point group
  character(len=3) :: pgroup          ! Point group character label
  character(len=3) :: irreplbl(0:7,8) ! Labels of the irreps in each point group
  character(len=3) :: pglbls(8)       ! Labels of every point group
  
end module bitglobal
