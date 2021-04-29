!**********************************************************************
! Calculation of 1-RDMs for MRCI wavefunctions
!**********************************************************************

!######################################################################
! density_mrci: Controls the calculation of MRCI 1-electron reduced
!               density matrices
!######################################################################
#ifdef CBINDING
subroutine density_mrci(irrep,nroots,iroots,dmat,confscr,vecscr) &
     bind(c,name='density_mrci')
#else
subroutine density_mrci(irrep,nroots,iroots,dmat,confscr,vecscr)
#endif

  use constants
  use bitglobal
  use iomod
  use conftype
  use rdm
  
  implicit none

  ! Irrep and no. roots
  integer(is), intent(in) :: irrep,nroots

  ! Indices of the states for which 1-RDMs are requested
  integer(is), intent(in) :: iroots(nroots)
  
  ! Density matrix
  real(dp), intent(out)   :: dmat(nmo,nmo,nroots)

  ! MRCI configuration and eigenvector scratch file numbers
  integer(is), intent(in) :: confscr(0:nirrep-1),vecscr(0:nirrep-1)
  
  ! MRCI configuration derived type
  type(mrcfg)             :: cfg

  ! Eigenpairs
  real(dp), allocatable   :: vec(:,:),ener(:)
  
  ! Everything else
  integer(is)             :: i

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  write(6,'(/,52a)') ('-',i=1,52)
  write(6,'(3(x,a))') 'Density matrix calculation for the',&
       trim(irreplbl(irrep,ipg)),'subspace'
  write(6,'(52a)') ('-',i=1,52)
  
!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr(irrep))

  write(6,'(/,x,a,x,i0)') 'CSF basis dimension:',cfg%csfdim

!----------------------------------------------------------------------
! Read in the eigenvectors
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(vec(cfg%csfdim,nroots))
  allocate(ener(nroots))
  
  ! Read in the eigenvectors
  call read_some_eigenpairs(vecscr(irrep),vec,ener,cfg%csfdim,&
       nroots,iroots)
  
!----------------------------------------------------------------------
! Compute the 1-RDMs
!----------------------------------------------------------------------
  call rdm_mrci(cfg,cfg%csfdim,nroots,vec,dmat)
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  call cfg%finalise
  deallocate(vec)
  deallocate(ener)
  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
    
end subroutine density_mrci
