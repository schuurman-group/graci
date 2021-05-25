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
subroutine mrci_prune(Athrsh,irrep,nroots,nextra,confscr,vec0scr,nconf) &
     bind(c,name="mrci_prune")
#else
subroutine mrci_prune(Athrsh,irrep,nroots,nextra,confscr,vec0scr,nconf)
#endif
    
  use constants
  use bitglobal
  use conftype
  use hii
  use epstein_nesbet
  use trimconfs
  use utils
  use iomod
  use timing
  
  implicit none

  ! Configuration selection threshold
  real(dp), intent(in)       :: Athrsh
  
  ! Irrep number
  integer(is), intent(in)    :: irrep

  ! Number of roots requested
  integer(is), intent(in)    :: nroots(0:nirrep-1)

  ! Number of extra roots to include in the ENPT2 calculation
  integer(is), intent(in)    :: nextra
    
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

  ! Zeroth-order energies
  real(dp), allocatable      :: E0(:)

  ! 2nd-order corrected energies
  real(dp), allocatable      :: EPT2(:)
  
  ! Surviving configuration flags
  integer(is), allocatable   :: i1I(:),i2I(:),i1E(:),i2E(:),i1I1E(:)
  
  ! Everything else
  integer(is)                :: i,nvec
  integer(is), allocatable   :: indx(:)

  ! Timing variables
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                                twall_end

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)
    
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
  ! Total number of roots for which the ENPT2 corrections will be
  ! calculated
  nvec=nroots(irrep)+nextra

  allocate(hdiag(cfg%csfdim))
  hdiag=0.0d0

  allocate(averageii(cfg%confdim))
  averageii=0.0d0

  allocate(Avec(cfg%csfdim,nvec))
  Avec=0.0d0

  allocate(E0(nvec))
  E0=0.0d0
  
  allocate(E2(nvec))
  E2=0.0d0

  allocate(EPT2(nvec))
  EPT2=0.0d0

  allocate(indx(nvec))
  indx=0
   
  allocate(i1I(cfg%n1I), i2I(cfg%n2I), i1E(cfg%n1E), i2E(cfg%n2E),&
       i1I1E(cfg%n1I1E))
  i1I=0; i2I=0; i1E=0; i2E=0; i1I1E=0

!----------------------------------------------------------------------
! Read in the zeroth-order (i.e., ref space) energies
!----------------------------------------------------------------------
  call read_energies(vec0scr(irrep),nvec,E0)

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
       vec0scr(irrep),Avec,E2,nvec)

!----------------------------------------------------------------------
! Sort the 2nd-order corrected energies
!----------------------------------------------------------------------
  EPT2=E0+E2
  
  call dsortindxa1('A',nvec,EPT2,indx)

!----------------------------------------------------------------------
! Get the indices of the configurations that generate CSFs with
! A-vector elements above threshold
!----------------------------------------------------------------------
  call trim_conf_indices(cfg,Athrsh,Avec,cfg%csfdim,cfg%confdim,&
       nroots(irrep),nvec,indx(1:nroots(irrep)),i1I,i2I,i1E,i2E,i1I1E,&
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
! Output the ENPT2 2nd-order corrected excitation energies
!----------------------------------------------------------------------
  write(6,'(/,x,a,/)') 'ENPT2 energies:'
  do i=1,nroots(irrep)
     write(6,'(3x,a,x,i3,a,2(2x,F12.6),x,a)') &
          'State' ,i,':',EPT2(indx(i)), &
          (EPT2(indx(i))-EPT2(indx(1)))*eh2ev,'eV' 
  enddo
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(hdiag)
  deallocate(averageii)
  deallocate(Avec)
  deallocate(E0)
  deallocate(E2)
  deallocate(EPT2)
  deallocate(indx)
  deallocate(i1I, i2I, i1E, i2E, i1I1E)

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'mrci_prune')
  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)

  return
    
end subroutine mrci_prune

