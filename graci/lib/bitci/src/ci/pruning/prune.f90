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
subroutine mrci_prune(Athrsh,irrep,nroots,nextra,confscr,vec0scr,&
     nconf,eqscr) bind(c,name="mrci_prune")
#else
subroutine mrci_prune(Athrsh,irrep,nroots,nextra,confscr,vec0scr,&
     nconf,eqscr)
#endif
    
  use constants
  use bitglobal
  use conftype
  use confinfo
  use hii
  use epstein_nesbet
  use pspace
  use qspace
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

  ! Q-space energy correction scratch file numbers
  integer(is), intent(out)   :: eqscr(0:nirrep-1)
  
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
  real(dp), allocatable      :: E2Q(:)
  
  ! Zeroth-order energies
  real(dp), allocatable      :: E0(:)

  ! 2nd-order corrected energies
  real(dp), allocatable      :: EPT2(:)
  
  ! Surviving configuration flags
  integer(is), allocatable   :: i1I(:),i2I(:),i1E(:),i2E(:),i1I1E(:)

  ! Active MO information
  integer(is)                :: nactive
  integer(is)                :: active(nmo)
  
  ! Everything else
  integer(is)                :: i,nvec
  integer(is), allocatable   :: indx(:)
  integer(is)                :: iscratch
  character(len=60)          :: eqfile
  character(len=2)           :: amult,airrep
  
  ! Timing variables
  real(dp)                   :: tcpu_start,tcpu_end,twall_start,&
                                twall_end
  
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  ! Section header
  if (verbose) then
     write(6,'(/,52a)') ('-',i=1,52)
     write(6,'(3(x,a))') 'Pruning of the',trim(irreplbl(irrep,ipg)),&
          'subspace'
     write(6,'(52a)') ('-',i=1,52)
     
     ! Configuration selection threshold
     write(6,'(/,x,a,x,ES10.4)') 'Selection theshold:',Athrsh

     ! Original no. configurations
     write(6,'(/,x,a,x,i0)') 'Original number of MRCI configurations:',&
          nconf(irrep)
  endif
     
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

  allocate(E2Q(nvec))
  E2Q=0.0d0
  
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
  call enpt2(irrep,cfg,hdiag,averageii,cfg%csfdim,cfg%confdim,&
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
  call pspace_conf_indices(cfg,Athrsh,Avec,cfg%csfdim,cfg%confdim,&
       nroots(irrep),nvec,indx(1:nroots(irrep)),i1I,i2I,i1E,i2E,i1I1E,&
       cfg%n1I,cfg%n2I,cfg%n1E,cfg%n2E,cfg%n1I1E)

!----------------------------------------------------------------------
! Calculate the contributions of the discarded CSFs to the 2nd-order
! energy corrections
!----------------------------------------------------------------------
! Important: this has to be done before calling pspace_set_new_confs
!----------------------------------------------------------------------
  call qspace_e2(cfg,cfg%csfdim,nvec,Avec,hdiag,E0,i1I,i2I,i1E,i2E,&
       i1I1E,cfg%n1I,cfg%n2I,cfg%n1E,cfg%n2E,cfg%n1I1E,E2Q)
  
!----------------------------------------------------------------------
! Fill in the surviving configuration information
!----------------------------------------------------------------------
  call pspace_set_new_confs(cfg,cfg_new,i1I,i2I,i1E,i2E,i1I1E,cfg%n1I,&
       cfg%n2I,cfg%n1E,cfg%n2E,cfg%n1I1E,confscr(irrep),nconf(irrep))

!----------------------------------------------------------------------
! Output the new no. configurations
!----------------------------------------------------------------------
  if (verbose) &
       write(6,'(/,x,a,x,i0)') 'New number of MRCI configurations:',&
       nconf(irrep)

!----------------------------------------------------------------------
! Determine and output the new no. active MOs
!----------------------------------------------------------------------
  call get_active_mos(cfg,nactive,active)

  if (verbose) &
       write(6,'(/,x,a,x,i0)') 'New number of active MOs:',nactive
  
!----------------------------------------------------------------------
! Output the ENPT2 2nd-order corrected excitation energies
!----------------------------------------------------------------------
  ! Q-space contributions to the ENPT2 energies
  if (verbose) then
     write(6,'(/,x,a,/)') 'Q-space energy corrections:'
     do i=1,nroots(irrep)
        write(6,'(3x,a,x,i3,a,2x,F12.6)') &
             'State' ,i,':',E2Q(indx(i))
     enddo
  endif
     
!----------------------------------------------------------------------
! Write the Q-space contributions to the ENPT2 energies to disk
!----------------------------------------------------------------------
  ! Register the scratch file
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name('e2q.mult'//trim(amult)//&
       '.sym'//trim(airrep),eqfile)
  call register_scratch_file(eqscr(irrep),eqfile)

  ! Open the scratch file
  iscratch=scrunit(eqscr(irrep))
  open(iscratch,file=scrname(eqscr(irrep)),form='unformatted',&
       status='unknown')

  ! No. vectors
  write(iscratch) nvec

  ! Q-space energy corrections
  write(iscratch) E2Q

  ! Close the scratch file
  close(iscratch)
  
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
    if (verbose) &
         call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'mrci_prune')
  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)

  return
    
end subroutine mrci_prune

!######################################################################
! retrieve_qcorr: Returns the Q-space energy corrections for a pruned
!                 MRCI calculation. Note that this is simply a publicly
!                 available gateway to the get_qcorr subroutine of the
!                 qspace module.
!######################################################################
#ifdef CBINDING
subroutine retrieve_qcorr(irrep,vecscr,vec0scr,confscr,eqscr,nroots,&
     nextra,qcorr,max_overlap) bind(c,name="retrieve_qcorr")
#else
subroutine retrieve_qcorr(irrep,vecscr,vec0scr,confscr,eqscr,nroots,&
     nextra,qcorr,max_overlap)
#endif
  
  use constants
  use bitglobal
  use qspace
  
  implicit none
  
  ! Irrep number
  integer(is), intent(in)  :: irrep
  
  ! MRCI and ref space eigenvector scratch file numbers
  integer(is), intent(in)  :: vecscr,vec0scr
  
  ! MRCI configuration scratch file number
  integer(is), intent(in)  :: confscr
  
  ! Q-space energy correction scratch file number
  integer(is), intent(in)  :: eqscr

  ! No. MRCI roots
  integer(is), intent(in)  :: nroots

  ! No. extra ref space roots
  integer(is), intent(in)  :: nextra

  ! Q-space energy corrections
  real(dp), intent(out)    :: qcorr(nroots)
  
  ! <Ref|MRCI> overlaps
  real(dp), intent(out)    :: max_overlap(nroots)
  
  call get_qcorr(irrep,vecscr,vec0scr,confscr,eqscr,nroots,&
       nextra,qcorr,max_overlap)
  
  return
  
end subroutine retrieve_qcorr
