#ifdef CBINDING
subroutine truncate_mrci_wf(irrep,nroots,confscr,vecscr,thrsh) &
     bind(c,name="truncate_mrci_wf")
#else
subroutine truncate_mrci_wf(irrep,nroots,confscr,vecscr,thrsh)
#endif

  use constants
  use bitglobal
  use conftype
  use pspace
  use iomod
    
  implicit none

  ! Irrep number
  integer(is), intent(in)  :: irrep

  ! Number of roots requested
  integer(is), intent(in)  :: nroots

  ! Array of MRCI configuration scratch file numbers
  integer(is), intent(in)  :: confscr

  ! Eigenpair scratch file number
  integer(is), intent(in)  :: vecscr

  ! Wave function truncation threshold
  real(dp), intent(in)     :: thrsh

  ! Number of surviving MRCI configurations
  integer(is)              :: nconf
  
  ! MRCI configuration derived type
  type(mrcfg)              :: cfg,cfg_new

  ! Wave function vectors
  real(dp), allocatable    :: vec(:,:)
  
  ! Surviving configuration flags
  integer(is), allocatable :: i1I(:),i2I(:),i1E(:),i2E(:),i1I1E(:)

  ! Everything else
  integer(is)              :: i
  integer(is), allocatable :: iroots(:)
  real(dp), allocatable    :: ener(:)
  
!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(vec(cfg%csfdim,nroots))
  vec=0.0d0

  allocate(ener(nroots))
  ener=0.0d0

  allocate(iroots(nroots))
  iroots=0
  
  allocate(i1I(cfg%n1I), i2I(cfg%n2I), i1E(cfg%n1E), i2E(cfg%n2E),&
       i1I1E(cfg%n1I1E))
  i1I=0; i2I=0; i1E=0; i2E=0; i1I1E=0
  
!----------------------------------------------------------------------
! Read in the wave function vectors
!----------------------------------------------------------------------
! Note that read_some_eigenpairs has to be called here in case extra
! roots were used in the MR-ENPT2 calculation
!----------------------------------------------------------------------
  ! Roots of interest
  do i=1,nroots
     iroots(i)=i
  enddo

  ! Read the wave functions from disk
  call read_some_eigenpairs(vecscr,vec,ener,cfg%csfdim,nroots,iroots)

!----------------------------------------------------------------------
! Get the indices of the configurations that generate CSFs with
! coefficients above threshold
!----------------------------------------------------------------------
  call pspace_conf_indices(cfg,thrsh,vec,cfg%csfdim,cfg%confdim,&
       nroots,nroots,iroots,i1I,i2I,i1E,i2E,i1I1E,cfg%n1I,cfg%n2I,&
       cfg%n1E,cfg%n2E,cfg%n1I1E)

!----------------------------------------------------------------------
! Fill in the surviving configuration information
!----------------------------------------------------------------------
  call pspace_set_new_confs(cfg,cfg_new,i1I,i2I,i1E,i2E,i1I1E,cfg%n1I,&
       cfg%n2I,cfg%n1E,cfg%n2E,cfg%n1I1E,confscr,nconf)

!----------------------------------------------------------------------
! Truncate the wave functions
!----------------------------------------------------------------------
  print*,''
  print*,'The wf vector truncation routines needs to be written...'
  stop
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(vec, ener, iroots, i1I, i2I, i1E, i2E, i1I1E)
    
  return
  
end subroutine truncate_mrci_wf
