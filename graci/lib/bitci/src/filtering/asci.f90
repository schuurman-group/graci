!**********************************************************************
! Additional filtering of the MRCI configurations based on their
! interaction with the reference space configurations
!**********************************************************************

!######################################################################
! filter_asci: Filtering out of MRCI configurations based on the
!              Adaptive Sampling Configuration Interaction (ASCI)
!              selection algorithm. For a description of this, see
!              J. Chem. Phys., 145, 044112 (2016).
!######################################################################
#ifdef CBINDING
subroutine filter_asci(Athrsh,irrep,nroots,confscr,vec0scr,nconf) &
     bind(c,name="filter_asci")
#else
subroutine filter_asci(Athrsh,irrep,nroots,confscr,vec0scr,nconf)
#endif
    
  use constants
  use bitglobal
  use conftype
  use hii
  use asci_avec
  use asci_trim
  
  implicit none

  ! Configuration selection threshold
  real(dp), intent(in)       :: Athrsh
  
  ! Irrep number
  integer(is), intent(in)    :: irrep

  ! Number of roots requested
  integer(is), intent(in)    :: nroots(0:nirrep-1)

  ! MRCI space configuration scratch file numbers
  integer(is), intent(in)    :: confscr(0:nirrep-1)

  ! Reference space eigenvector scratch file numbers
  integer(is), intent(in)    :: vec0scr(0:nirrep-1)
  
  ! Number of MRCI configurations
  integer(is), intent(inout) :: nconf(0:nirrep-1)

  ! MRCI configuration derived types
  type(mrcfg)                :: cfg,cfg_new
  
  ! On-diagonal Hamiltonian matrix elements
  real(dp), allocatable      :: hdiag(:)

  ! On-diagonal Hamiltonian matrix elements
  ! averaged over spin couplings
  real(dp), allocatable      :: averageii(:)

  ! A-vectors
  real(dp), allocatable      :: Avec(:,:)

  ! Surviving configuration flags
  integer(is), allocatable   :: i1I(:),i2I(:),i1E(:),i2E(:),i1I1E(:)
  
  ! Everything else
  integer(is)                :: i
  
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  ! Section header
  write(6,'(/,52a)') ('-',i=1,52)
  write(6,'(3(x,a))') 'ASCI configuration filtering for the',&
       trim(irreplbl(irrep,ipg)),'subspace'
  write(6,'(52a)') ('-',i=1,52)

  ! Configuration selection threshold
  write(6,'(/,x,a,x,ES10.4)') 'Selection theshold:',Athrsh

  ! Original no. configurations
  write(6,'(/,x,a,x,i0)') 'Original number of MRCI configurations:',&
       nconf(irrep)
  
!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr(irrep))

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(hdiag(cfg%csfdim))
  hdiag=0.0d0

  allocate(averageii(cfg%confdim))
  averageii=0.0d0

  allocate(Avec(cfg%csfdim,nroots(irrep)))
  Avec=0.0d0

  allocate(i1I(cfg%n1I), i2I(cfg%n2I), i1E(cfg%n1E), i2E(cfg%n2E),&
       i1I1E(cfg%n1I1E))
  i1I=0; i2I=0; i1E=0; i2E=0; i1I1E=0
  
!----------------------------------------------------------------------
! Compute the on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
  call hmat_diagonal(hdiag,cfg%csfdim,averageii,cfg%confdim,cfg)

!----------------------------------------------------------------------
! Compute the A-vectors
!
! A_w,omega,I = <w omega|H|Psi^0_I>/(H_{w omega,w omega} - E^0_I),
!
! where {|Psi^0_I>, E^0_I} is the set of reference space eigenpairs.
!----------------------------------------------------------------------
  call asci_avectors(cfg,hdiag,averageii,cfg%csfdim,cfg%confdim,&
       vec0scr(irrep),Avec,nroots(irrep))

!----------------------------------------------------------------------
! Get the indices of the configurations that generate CSFs with
! A-vector elements above threshold
!----------------------------------------------------------------------
  call asci_conf_indices(cfg,Athrsh,Avec,cfg%csfdim,cfg%confdim,&
       nroots(irrep),i1I,i2I,i1E,i2E,i1I1E,&
       cfg%n1I,cfg%n2I,cfg%n1E,cfg%n2E,cfg%n1I1E)

!----------------------------------------------------------------------
! Fill in the surviving configuration information
!----------------------------------------------------------------------
  call asci_set_new_confs(cfg,cfg_new,i1I,i2I,i1E,i2E,i1I1E,&
       cfg%n1I,cfg%n2I,cfg%n1E,cfg%n2E,cfg%n1I1E,confscr(irrep),&
       nconf(irrep))

!----------------------------------------------------------------------
! Output the new no. configurations
!----------------------------------------------------------------------
  write(6,'(/,x,a,x,i0)') 'New number of MRCI configurations:',&
       nconf(irrep)
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(hdiag)
  deallocate(averageii)
  deallocate(Avec)
  deallocate(i1I, i2I, i1E, i2E, i1I1E)
  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)

  return
    
end subroutine filter_asci

!######################################################################
