!**********************************************************************
! Additional pruning of the MRCI configurations based on their
! contributions to the ENPT2 1st-order corrected wave functions
!**********************************************************************

!######################################################################
! mrci_prune: Pruning of the MRCI configuration space using the
!             Epstein-Nesbet 1st-order wave function corrections as
!             the pruning metric.
!######################################################################
#ifdef CBINDING
subroutine mrci_prune(Athrsh,irrep,nroots,confscr,vec0scr,nconf) &
     bind(c,name="mrci_prune")
#else
subroutine mrci_prune(Athrsh,irrep,nroots,confscr,vec0scr,nconf)
#endif
    
  use constants
  use bitglobal
  use conftype
  use hii
  use epstein_nesbet
  use trimconfs
  
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

  ! ENPT2 energy and wavefunction corrections
  real(dp), allocatable      :: Avec(:,:)
  real(dp), allocatable      :: E2(:)
  
  ! Surviving configuration flags
  integer(is), allocatable   :: i1I(:),i2I(:),i1E(:),i2E(:),i1I1E(:)
  
  ! Everything else
  integer(is)                :: i
  
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  ! Section header
  write(6,'(/,52a)') ('-',i=1,52)
  write(6,'(3(x,a))') 'Pruning of the',trim(irreplbl(irrep,ipg)),&
       'subspace'
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

  allocate(E2(nroots(irrep)))
  E2=0.0d0
  
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
  call enpt2(cfg,hdiag,averageii,cfg%csfdim,cfg%confdim,&
       vec0scr(irrep),Avec,E2,nroots(irrep))

!----------------------------------------------------------------------
! Get the indices of the configurations that generate CSFs with
! A-vector elements above threshold
!----------------------------------------------------------------------
  call trim_conf_indices(cfg,Athrsh,Avec,cfg%csfdim,cfg%confdim,&
       nroots(irrep),i1I,i2I,i1E,i2E,i1I1E,&
       cfg%n1I,cfg%n2I,cfg%n1E,cfg%n2E,cfg%n1I1E)

!----------------------------------------------------------------------
! Fill in the surviving configuration information
!----------------------------------------------------------------------
  call trim_set_new_confs(cfg,cfg_new,i1I,i2I,i1E,i2E,i1I1E,cfg%n1I,&
       cfg%n2I,cfg%n1E,cfg%n2E,cfg%n1I1E,confscr(irrep),nconf(irrep))

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
    
end subroutine mrci_prune

!######################################################################
