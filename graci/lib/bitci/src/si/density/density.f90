!**********************************************************************
! Calculation of 1-RDMs and 1-TDMs for MRCI wavefunctions
!**********************************************************************

!######################################################################
! density_mrci: Controls the calculation of MRCI 1-electron reduced
!               density matrices
!######################################################################
#ifdef CBINDING
subroutine density_mrci(irrep,nroots,iroots,dmat,conffile_in,vecfile_in) &
     bind(c,name='density_mrci')
#else
subroutine density_mrci(irrep,nroots,iroots,dmat,conffile_in,vecfile_in)
#endif

  use constants
  use bitglobal
  use iomod
  use iso_c_binding, only: C_CHAR
  use conftype
  use rdm
  
  implicit none

  ! MRCI configuration and eigenvector file names
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: conffile_in(*),vecfile_in(*)
  character(len=255)                 :: conffile,vecfile
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: conffile_in,vecfile_in
  character(len=255)                 :: conffile,vecfile
#endif

  ! Irrep and no. roots
  integer(is), intent(in) :: irrep,nroots

  ! Indices of the states for which 1-RDMs are requested
  integer(is), intent(in) :: iroots(nroots)
  
  ! Density matrix
  real(dp), intent(out)   :: dmat(nmo,nmo,nroots)

  ! MRCI configuration derived type
  type(mrcfg)             :: cfg

  ! Scratch file numbers
  integer(is)             :: confscr,vecscr

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
! If C bindings are on, then convert the MRCI configuration and
! eigenvector file names from the C char type to the Fortran character
! type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(conffile_in)
  call c2fstr(conffile_in,conffile,length)
  length=cstrlen(vecfile_in)
  call c2fstr(vecfile_in,vecfile,length)
#else
  conffile=adjustl(trim(conffile_in))
  vecfile=adjustl(trim(conffile_in))
#endif

!----------------------------------------------------------------------
! Register the configuration and eigenvector scratch files
!----------------------------------------------------------------------
  call register_scratch_file(confscr,conffile)
  call register_scratch_file(vecscr,vecfile)

!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr)

  write(6,'(/,x,a,x,i0)') 'CSF basis dimension:',cfg%csfdim

!----------------------------------------------------------------------
! Read in the eigenvectors
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(vec(cfg%csfdim,nroots))
  allocate(ener(nroots))
  
  ! Read in the eigenvectors
  call read_some_eigenpairs(vecscr,vec,ener,cfg%csfdim,nroots,iroots)
  
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
