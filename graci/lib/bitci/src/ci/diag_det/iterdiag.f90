!**********************************************************************
! Gateway routines for the Davidson diagonalisation of a selected CI
! Hamiltonian matrix
!**********************************************************************

!######################################################################
! hdiag_spinpure: Iterative diagonalisation of a selected CI
!                 Hamiltonian matrix using spin purification
!######################################################################
#ifdef CBINDING
subroutine hdiag_spinpure(irrep,detscr,vecscr,nroots) &
     bind(c,name="hdiag_spinpure")
#else
subroutine hdiag_spinpure(irrep,detscr,vecscr,nroots)
#endif

  use constants
  use bitglobal
  use utils
  use iomod
  use detsort
  use sigma
  use guessvecs
  use timing
  
  implicit none

  ! Input irrep and requested number of roots
  integer(is), intent(in)    :: irrep
  integer(is), intent(inout) :: nroots
  
  ! Determinant and eigenvector scratch file number
  integer(is), intent(in)    :: detscr
  integer(is), intent(out)   :: vecscr

  ! Determinant arrays
  integer(is)                :: ndet,offdima,offdimb
  integer(ib), allocatable   :: da(:,:,:),db(:,:,:)
  integer(is)                :: nsym(0:nirrep-1)
  integer(is)                :: nsym_sum(0:nirrep-1)
  integer(is)                :: nuniquea(0:nirrep-1),&
                                nuniqueb(0:nirrep-1)
  integer(is), allocatable   :: offseta(:,:),offsetb(:,:)
  integer(is), allocatable   :: mapab(:)

  ! On-diagonal elements of the Hamiltonian matrix
  real(dp), allocatable      :: hdiag(:)

  ! Guess vectors
  integer(is)                :: guessdim
  
  ! Timing variables
  real(dp)                   :: tcpu_start,tcpu_end,twall_start,&
                                twall_end
  
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)
  
!----------------------------------------------------------------------
! Read in the determinants and associated information from disk
!----------------------------------------------------------------------
  call read_det_file_sorted(detscr,ndet,offdima,offdimb,nsym,da,db,&
       nuniquea,nuniqueb,offseta,offsetb,mapab)

!----------------------------------------------------------------------
! Return if there are no determinants for the current irrep
!----------------------------------------------------------------------
  if (nsym(irrep) == 0) then
     write(6,'(/,x,a)') 'No MRCI determinants of '&
          //trim(irreplbl(irrep,ipg))//' symmetry'
     vecscr=-1
     deallocate(da)
     deallocate(db)
     deallocate(offseta)
     deallocate(offsetb)
     return
  endif

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(hdiag(nsym(irrep)))
  hdiag=0.0d0
  
!----------------------------------------------------------------------
! Calculate the on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
  call h_ondiag(irrep,hdiag,ndet,offdima,nsym,da,nuniquea,offseta)

!----------------------------------------------------------------------
! Determine the guess vector space
!----------------------------------------------------------------------
  ! Requested guess vector space dimension. Note that this will
  ! be increased in guess_space as needed to obtain spin-pure
  ! guess vectors
  guessdim=min(400,nsym(irrep))

  ! Generate the guess space
  call guess_space(hdiag,nsym,irrep,guessdim,da,ndet,offseta,offdima)
  
!----------------------------------------------------------------------
! Deallocate arrays  
!----------------------------------------------------------------------
  deallocate(da)
  deallocate(db)
  deallocate(offseta)
  deallocate(offsetb)
  deallocate(mapab)
  deallocate(hdiag)

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine hdiag_spinpure
