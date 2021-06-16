!**********************************************************************
! Extraction of the determinant expansions of the MRCI wave functions
!**********************************************************************

!######################################################################
! wf_mrci: Controls the calculation and writing to an HDF5 file of a
!          a requested subset of the MRCI wave functions for a single
!          irrep
!######################################################################
#ifdef CBINDING
subroutine wf_mrci(irrep,nroots,iroots,confscr,vecscr,h5file_in,ndet) &
     bind(c,name='wf_mrci')
#else
subroutine wf_mrci(irrep,nroots,iroots,confscr,vecscr,h5file_in,ndet)
#endif

  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use iomod
  use conftype
  use csf2det
  
  implicit none

  ! Output HDF5 file name
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: h5file_in
  character(len=255)                 :: h5file
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: h5file_in
  character(len=255)                 :: h5file
#endif
  
  ! Irrep and no. roots
  integer(is), intent(in)  :: irrep,nroots

  ! Indices of the states for which wave functions are requested
  integer(is), intent(in) :: iroots(nroots)
  
  ! MRCI configuration and eigenvector scratch file numbers
  integer(is), intent(in)  :: confscr(0:nirrep-1),vecscr(0:nirrep-1)

  ! No. determinants
  integer(is), intent(out) :: ndet
  
  ! MRCI configuration derived type
  type(mrcfg)              :: cfg

  ! Eigenpairs
  real(dp), allocatable    :: vec(:,:),ener(:)

  ! Everything else
  integer(is)              :: i
  
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  write(6,'(/,52a)') ('-',i=1,52)
  write(6,'(3(x,a))') 'Wave function extraction for the',&
       trim(irreplbl(irrep,ipg)),'subspace'
  write(6,'(52a)') ('-',i=1,52)
  
!----------------------------------------------------------------------
! If C bindings are on, then convert the output HDF5 filename from the
! C char type to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(h5file_in)
  call c2fstr(h5file_in,h5file,length)
#else
  h5file=adjustl(trim(h5file_in))
#endif

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
! Compute the determinant expansion of the MRCI wave functions and
! save them to disk
!----------------------------------------------------------------------
  call csf2det_mrci(cfg,nroots,cfg%csfdim,vec,ndet)
  
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
  
end subroutine wf_mrci
