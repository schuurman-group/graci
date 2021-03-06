!**********************************************************************
! Routines for the calculation of the ENPT2 energy and wave function
! corrections:
!
! E_w,omega = -|<w omega|H|Psi^0_I>|^2/(H_{w omega,w omega} - E^0_I),
!
! A_w,omega = <w omega|H|Psi^0_I>/(H_{w omega,w omega} - E^0_I),
!
! where {|Psi^0_I>, E^0_I} is the set of reference space eigenpairs.
!**********************************************************************
module epstein_nesbet

  use constants

  implicit none

  ! Temporary Hij array
  integer(is), private           :: harr2dim
  real(dp), allocatable, private :: harr2(:)

contains

!######################################################################
! enpt2: Computes a batch of ENPT2 energy and wave function corrections
!######################################################################
  subroutine enpt2(irrep,cfg,hdiag,averageii,csfdim,confdim,vec0scr,&
       Avec,E2,nroots,shift,dspscr,multistate,EQD,mix)

    use constants
    use bitglobal
    use conftype
    use iomod

    implicit none

    ! Irrep
    integer(is)                     :: irrep
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in)         :: cfg

    ! No. CSFs and configurations
    integer(is), intent(in)         :: csfdim,confdim
    
    ! On-diagonal Hamiltonian matrix elements and their
    ! spin-coupling averaged values
    real(dp), intent(in)            :: hdiag(csfdim),averageii(confdim)

    ! Reference space eigenpair scratch file number
    integer(is), intent(in)         :: vec0scr

    ! Number of roots
    integer(is), intent(in)         :: nroots

    ! ENPT2 wave function and energy corrections
    real(dp), intent(out)           :: Avec(csfdim,nroots)
    real(dp), intent(out)           :: E2(nroots)

    ! ISA shift
    real(dp), intent(in)            :: shift
    
    ! Multistate flag
    logical, intent(in)             :: multistate

    ! Damped strong perturber scratch file number
    integer, intent(out)            :: dspscr
    
    ! Multistate energies and mixing coefficients
    real(dp), optional, intent(out) :: EQD(nroots),mix(nroots,nroots)
    
    ! Reference space eigenpairs
    integer(is)                     :: refdim
    integer(is), allocatable        :: iroots(:)
    real(dp), allocatable           :: e0(:),vec0(:,:)

    ! Damped strong perturbers
    integer(is)                     :: ndsp
    integer(is), allocatable        :: idsp(:)

    ! I/O
    integer(is)                     :: iscratch
    character(len=60)               :: dspfile
    character(len=2)                :: amult,airrep
        
    ! Everything else
    integer(is)                     :: i

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Number of reference space CSFs
    refdim=cfg%csfs0h(cfg%n0h+1)-1

    ! Reference space eigenvalues
    allocate(e0(nroots))
    e0=0.0d0

    ! Reference space eigenvectors
    allocate(iroots(nroots))
    allocate(vec0(refdim,nroots))
    iroots=0
    vec0=0.0d0

    ! Hij working array
    harr2dim=maxval(ncsfs(0:nomax))**2    
    allocate(harr2(harr2dim))

    ! Indices of the damped strong perturbers
    allocate(idsp(csfdim))
    idsp=0
    
!----------------------------------------------------------------------
! Reference space eigenpairs
!----------------------------------------------------------------------
! Note that we have to use read_some_eigenpairs here as different
! numbers of extra ref space eigenvectors may be used for
! different tasks. We will, however, assume for now that the 1st
! nroots ref space eigenvectors are needed.
!----------------------------------------------------------------------
    ! Read in the eigenpairs
    do i=1,nroots
       iroots(i)=i
    enddo
    call read_some_eigenpairs(vec0scr,vec0,e0,refdim,nroots,iroots)

    ! Subtract off E_SCF from the energies to get the true eigenvalues
    e0=e0-escf
    
!----------------------------------------------------------------------
! (1) 1-hole configurations -> 1I and 1E configurations
!----------------------------------------------------------------------
    call avec_1h(cfg,Avec,averageii,vec0,csfdim,confdim,refdim,nroots)
    
!----------------------------------------------------------------------
! (2)  2-hole configurations -> 2I, 2E and 1I1E configurations
!----------------------------------------------------------------------
    call avec_2h(cfg,Avec,averageii,vec0,csfdim,confdim,refdim,nroots)

!----------------------------------------------------------------------
! QDPT2 energies and mixing coefficients
!----------------------------------------------------------------------
    if (multistate) call qdpt2(cfg,Avec,hdiag,e0,EQD,mix,csfdim,&
         nroots,refdim,shift)
    
!----------------------------------------------------------------------
! Divide by (H_nn - E^0_I)
!----------------------------------------------------------------------
    call apply_denominator(cfg,Avec,E2,hdiag,e0,csfdim,nroots,refdim,&
         shift,idsp)

!----------------------------------------------------------------------
! Write the damped strong perturber array to disk
!----------------------------------------------------------------------
    ! Register the scratch file
    write(amult,'(i0)') imult
    write(airrep,'(i0)') irrep
    call scratch_name('dsp'//'.mult'//trim(amult)//&
         '.sym'//trim(airrep),dspfile)
    call register_scratch_file(dspscr,dspfile)

    ! Open scratch file
    iscratch=scrunit(dspscr)
    open(iscratch,file=scrname(dspscr),form='unformatted',&
         status='unknown')

    ! Number of damped strong perturbers
    ndsp=sum(idsp)
    write(iscratch) ndsp

    ! Damped strong perturber flags
    write(iscratch) idsp
    
    ! Close scratch file
    close(iscratch)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(e0)
    deallocate(iroots)
    deallocate(vec0)
    deallocate(harr2)

    return
    
  end subroutine enpt2
    
!######################################################################
! avec_1h: Calculation of the elements of the A-vector corresponding
!          to the CSFs generated by the 1-hole configurations
!######################################################################
  subroutine avec_1h(cfg,Avec,averageii,vec0,csfdim,confdim,refdim,&
       nroots)

    use constants
    use bitglobal
    use conftype
    use mrci_integrals
    use mrciutils
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Dimensions
    integer(is), intent(in) :: csfdim,confdim,refdim,nroots

    ! A-vector
    real(dp), intent(inout) :: Avec(csfdim,nroots)

    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: averageii(confdim)
    
    ! Reference space eigenvectors
    real(dp), intent(in)    :: vec0(refdim,nroots)

    ! Difference configuration information
    integer(is)             :: ndiff
    integer(is)             :: Dw(nmo,2)
    
    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)

    ! Indices of the singly-occupied ket MOs
    integer(is)             :: socc(nmo)
    integer(is)             :: nsocc
        
    ! Working arrays
    integer(ib)             :: kconf_full(n_int,2)
    integer(ib)             :: ksop_full(n_int,2)
    
    ! Everything else
    integer(is)             :: ibconf1I,ibconf1E,kconf
    integer(is)             :: knopen,nac
    integer(is)             :: n,bnsp,knsp
    integer(is)             :: n_int_I

!----------------------------------------------------------------------
! Elements of the A-vector corresponding to 1I and 1E configurations
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I

    ! Loop over ket reference configurations
    do kconf=1,cfg%n0h

       ! Number of open shells in the ket configuration
       knopen=sop_nopen(cfg%sop0h(:,:,kconf),n_int_I)

       ! Number of ket CSFs
       knsp=ncsfs(knopen)
       
       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_I,:)=cfg%conf0h(:,:,kconf)
       ksop_full(1:n_int_I,:)=cfg%sop0h(:,:,kconf)

       ! Package the ket configuration information
       call package_confinfo_offdiag(ksop_full,kconf_full,socc,nsocc,&
            Dw,ndiff,nbefore)
       
       ! Loop over 1-hole configurations
       ibconf1I=0
       ibconf1E=0
       do n=1,cfg%n1h

          ! Number of creation and annihilation operators linking the
          ! reference and 1-hole configurations
          nac=n_create_annihilate(cfg%conf0h(1:n_int_I,:,kconf), &
               cfg%conf1h(1:n_int_I,:,n),n_int_I)

          ! 1I elements
          if (nac <= 5) then
             if (cfg%n1I > 0) then
                call avec_1I(n,kconf,ibconf1I,ksop_full,ndiff,Dw,&
                     nbefore,socc,nsocc,knopen,knsp,n_int_I,cfg,&
                     averageii,csfdim,confdim,Avec,vec0,nroots,refdim)
             endif
          else
             ibconf1I=ibconf1I+cfg%off1I(n+1)-cfg%off1I(n)
          endif

          ! 1E elements
          if (nac <= 3) then
             if (cfg%n1E > 0) then
                call avec_1E(n,kconf,ibconf1E,n_int_I,&
                     kconf_full,ksop_full,ndiff,Dw,nbefore,socc,nsocc,&
                     knopen,knsp,averageii,confdim,cfg,&
                     csfdim,Avec,vec0,nroots,refdim)
             endif
          else
             ibconf1E=ibconf1E+cfg%off1E(n+1)-cfg%off1E(n)
          endif
          
       enddo
          
    enddo
       
    return
    
  end subroutine avec_1h

!######################################################################
! avec_2h: Calculation of the elements of the A-vector corresponding
!          to the CSFs generated by the 2-hole configurations
!######################################################################
  subroutine avec_2h(cfg,Avec,averageii,vec0,csfdim,confdim,refdim,&
       nroots)

    use constants
    use bitglobal
    use conftype
    use mrci_integrals
    use mrciutils
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Dimensions
    integer(is), intent(in) :: csfdim,confdim,refdim,nroots

    ! A-vector
    real(dp), intent(inout) :: Avec(csfdim,nroots)

    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: averageii(confdim)
    
    ! Reference space eigenvectors
    real(dp), intent(in)    :: vec0(refdim,nroots)

    ! Difference configuration information
    integer(is)             :: ndiff
    integer(is)             :: Dw(nmo,2)
    
    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)

    ! Indices of the singly-occupied ket MOs
    integer(is)             :: socc(nmo)
    integer(is)             :: nsocc
        
    ! Working arrays
    integer(ib)             :: kconf_full(n_int,2)
    integer(ib)             :: ksop_full(n_int,2)
    
    ! Everything else
    integer(is)             :: ibconf2I,ibconf2E,ibconf1I1E,kconf
    integer(is)             :: knopen,nac
    integer(is)             :: n,bnsp,knsp
    integer(is)             :: n_int_I

!----------------------------------------------------------------------
! Elements of the A-vector corresponding to 2I, 2E, and 1I1E
! configurations
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I
    
    ! Loop over ket reference configurations
    do kconf=1,cfg%n0h
       
       ! Number of open shells in the ket configuration
       knopen=sop_nopen(cfg%sop0h(:,:,kconf),n_int_I)
    
       ! Number of ket CSFs
       knsp=ncsfs(knopen)
       
       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_I,:)=cfg%conf0h(:,:,kconf)
       ksop_full(1:n_int_I,:)=cfg%sop0h(:,:,kconf)
    
       ! Package the ket configuration information
       call package_confinfo_offdiag(ksop_full,kconf_full,socc,nsocc,&
            Dw,ndiff,nbefore)
       
       ! Loop over 2-hole configurations
       ibconf2I=0
       ibconf2E=0
       ibconf1I1E=0
       do n=1,cfg%n2h
          
          ! Number of creation and annihilation operators linking the
          ! reference and 2-hole configurations
          nac=n_create_annihilate(cfg%conf0h(1:n_int_I,:,kconf),&
               cfg%conf2h(1:n_int_I,:,n),n_int_I)
          
          ! 2I elements
          if (nac <= 6) then
             if (cfg%n2I > 0) then
                call avec_2I(n,kconf,ibconf2I,n_int_I,&
                     ksop_full,ndiff,Dw,nbefore,socc,nsocc,knopen,&
                     knsp,averageii,confdim,cfg,csfdim,Avec,vec0,&
                     nroots,refdim)
             endif
          else
             ibconf2I=ibconf2I+cfg%off2I(n+1)-cfg%off2I(n)
          endif
             
          ! 2E elements
          if (nac <= 2) then
             if (cfg%n2E > 0) then
                call avec_2E(n,kconf,ibconf2E,n_int_I,&
                     kconf_full,ksop_full,ndiff,Dw,nbefore,socc,&
                     nsocc,knopen,knsp,averageii,confdim,cfg,csfdim,&
                     Avec,vec0,nroots,refdim)
             endif
          else
             ibconf2E=ibconf2E+cfg%off2E(n+1)-cfg%off2E(n)
          endif
          
          ! 1I1E elements
          if (nac <= 4) then
             if (cfg%n1I1E > 0) then
                call avec_1I1E(n,kconf,ibconf1I1E,n_int_I,&
                     kconf_full,ksop_full,ndiff,Dw,nbefore,socc,nsocc,&
                     knopen,knsp,averageii,confdim,cfg,csfdim,&
                     Avec,vec0,nroots,refdim)
             endif
          else
             ibconf1I1E=ibconf1I1E+cfg%off1I1E(n+1)-cfg%off1I1E(n)
          endif
       
       enddo

    enddo
       
    return
    
  end subroutine avec_2h
    
!######################################################################
! avec_1I: calculation of the A-vector elements corresponding to
!          1I configurations
!######################################################################
  subroutine avec_1I(n,kconf,ibconf1I,ksop_full,ndiff,Dw,nbefore,&
       socc,nsocc,knopen,knsp,n_int_I,cfg,averageii,csfdim,&
       confdim,Avec,vec0,nroots,refdim)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hbuild_mrci
    use mrci_integrals
    use dftmrci
    
    implicit none

    ! Bra configuration counter
    integer(is), intent(inout) :: ibconf1I

    ! Index of the 1-hole configuration
    integer(is), intent(in)    :: n

    ! Dimensions
    integer(is), intent(in)    :: n_int_I
    integer(is), intent(in)    :: csfdim,confdim,nroots,refdim
    
    ! Ket reference space configuration
    integer(is), intent(in)    :: kconf
    integer(ib), intent(in)    :: ksop_full(n_int,2)
    integer(is), intent(in)    :: knopen,knsp
    
    ! Difference configuration information
    integer(is),intent(in)     :: ndiff
    integer(is),intent(in)     :: Dw(nmo,2)
    
    ! Number of open shells preceding each MO
    integer(is), intent(in)    :: nbefore(nmo)

    ! Indices of the singly-occupied ket MOs
    integer(is), intent(in)    :: socc(nmo)
    integer(is), intent(in)    :: nsocc

    ! MRCI configuration derived type
    type(mrcfg), intent(in)    :: cfg

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)       :: averageii(confdim)

    ! A-vector
    real(dp), intent(inout)    :: Avec(csfdim,nroots)

    ! Reference space eigenvectors
    real(dp), intent(in)       :: vec0(refdim,nroots)
    
    ! Working arrays
    integer(ib)                :: bconf_int(n_int_I,2)
    integer(ib)                :: bsop_int(n_int_I,2)
    integer(is), parameter     :: maxexci=2
    integer(is)                :: hlist(maxexci),plist(maxexci)

    ! Everything else
    integer(is)                :: ioff,nexci,bnopen,bnsp

!----------------------------------------------------------------------
! Compute the 1I A-vector elements
!----------------------------------------------------------------------
    ! Loop over 1I configurations generated by the current
    ! 1-hole configuration
    do ioff=cfg%off1I(n),cfg%off1I(n+1)-1
    
       ! Increment the 1I configuration counter
       ibconf1I=ibconf1I+1

       ! Bra 1I configuration in the internal MO space
       bconf_int=0_ib
       bconf_int=cfg%conf1I(1:n_int_I,:,ibconf1I)
       bsop_int=cfg%sop1I(1:n_int_I,:,ibconf1I)
       
       ! Compute the excitation degree between the two
       ! configurations
       nexci=exc_degree_conf(cfg%conf0h(:,:,kconf),bconf_int,&
            n_int_I)
    
       ! Cycle if the excitation degree is greater than 2
       if (nexci > 2) cycle

       ! Number of open shells in the bra configuration
       bnopen=sop_nopen(bsop_int,n_int_I)
    
       ! Number of bra CSFs
       bnsp=ncsfs(bnopen)
    
       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(cfg%conf0h(:,:,kconf),bconf_int,&
            n_int_I,hlist(1:nexci),plist(1:nexci),nexci)
       
       ! Compute the matrix elements between the CSFs generated
       ! by the bra and ket configurations
       call hij_mrci(harr2,harr2dim,nexci,&
            ibconf1I,kconf,&
            cfg%sop1I(:,:,ibconf1I),ksop_full,&
            bnsp,knsp,bnopen,knopen,hlist,plist,cfg%m2c,&
            socc,nsocc,nbefore,Dw,ndiff,&
            cfg%csfs1I,cfg%csfs0h,cfg%n1I+1,cfg%n0h+1,&
            averageii(ibconf1I+cfg%n0h),&
            averageii(kconf))

       ! Contraction of the Hamiltonian matrix elements with
       ! the reference space eigenvectors
       call contract_hmat_vec0(ibconf1I,kconf,&
            cfg%csfs1I,cfg%csfs0h,cfg%n1I+1,cfg%n0h+1,&
            bnsp,knsp,Avec,vec0,csfdim,nroots,refdim)
       
    enddo
       
    return
    
  end subroutine avec_1I

!######################################################################
! avec_1E: calculation of the A-vector elements corresponding to
!          1E configurations
!######################################################################
  subroutine avec_1E(n,kconf,ibconf1E,n_int_I,kconf_full,&
       ksop_full,ndiff,Dw,nbefore,socc,nsocc,knopen,knsp,averageii,&
       confdim,cfg,csfdim,Avec,vec0,nroots,refdim)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hbuild_mrci
    use mrci_integrals
    use dftmrci
    
    implicit none

    ! Bra configuration counter
    integer(is), intent(inout) :: ibconf1E

    ! Index of the 1-hole configuration
    integer(is), intent(in)    :: n

    ! Dimensions
    integer(is), intent(in)    :: n_int_I
    integer(is), intent(in)    :: csfdim,confdim,nroots,refdim
    
    ! Ket reference space configuration
    integer(is), intent(in)    :: kconf
    integer(ib), intent(in)    :: kconf_full(n_int,2)
    integer(ib), intent(in)    :: ksop_full(n_int,2)
    integer(is), intent(in)    :: knopen,knsp
    
    ! Difference configuration information
    integer(is),intent(in)     :: ndiff
    integer(is),intent(in)     :: Dw(nmo,2)
    
    ! Number of open shells preceding each MO
    integer(is), intent(in)    :: nbefore(nmo)

    ! Indices of the singly-occupied ket MOs
    integer(is), intent(in)    :: socc(nmo)
    integer(is), intent(in)    :: nsocc

    ! MRCI configuration derived type
    type(mrcfg), intent(in)    :: cfg

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)       :: averageii(confdim)

    ! A-vector
    real(dp), intent(inout)    :: Avec(csfdim,nroots)

    ! Reference space eigenvectors
    real(dp), intent(in)       :: vec0(refdim,nroots)
    
    ! Working arrays
    integer(ib)                :: bconf_full(n_int,2)
    integer(ib)                :: bsop_full(n_int,2)
    integer(is), parameter     :: maxexci=2
    integer(is)                :: hlist(maxexci),plist(maxexci)

    ! Everything else
    integer(is)                :: ioff,nexci,bnopen,bnsp

!----------------------------------------------------------------------
! Compute the 1E A-vector elements
!----------------------------------------------------------------------
    ! Loop over 1E configurations generated by the current
    ! 1-hole configuration
    do ioff=cfg%off1E(n),cfg%off1E(n+1)-1
    
       ! Increment the 1E configuration counter
       ibconf1E=ibconf1E+1
    
       ! Bra 1E configuration in the full MO space
       bconf_full=0_ib
       bsop_full=0_ib
       bconf_full=cfg%conf1E(:,:,ibconf1E)
       bsop_full=cfg%sop1E(:,:,ibconf1E)
    
       ! Compute the excitation degree between the two
       ! configurations
       nexci=exc_degree_conf(kconf_full,bconf_full,n_int)
    
       ! Cycle if the excitation degree is greater than 2
       if (nexci > 2) cycle

       ! Number of open shells in the bra configuration
       bnopen=sop_nopen(bsop_full,n_int)
    
       ! Number of bra CSFs
       bnsp=ncsfs(bnopen)
    
       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(kconf_full,bconf_full,&
            n_int,hlist(1:nexci),plist(1:nexci),nexci)

       ! Compute the matrix elements between the CSFs generated
       ! by the bra and ket configurations
       call hij_mrci(harr2,harr2dim,nexci,&
            ibconf1E,kconf,&
            cfg%sop1E(:,:,ibconf1E),ksop_full,&
            bnsp,knsp,bnopen,knopen,hlist,plist,cfg%m2c,&
            socc,nsocc,nbefore,Dw,ndiff,&
            cfg%csfs1E,cfg%csfs0h,cfg%n1E+1,cfg%n0h+1,&
            averageii(ibconf1E+cfg%n0h+cfg%n1I+cfg%n2I),&
            averageii(kconf))

       call contract_hmat_vec0(ibconf1E,kconf,&
            cfg%csfs1E,cfg%csfs0h,cfg%n1E+1,cfg%n0h+1,&
            bnsp,knsp,Avec,vec0,csfdim,nroots,refdim)
       
    enddo

    return

  end subroutine avec_1E

!######################################################################
! avec_2I: calculation of the A-vector elements corresponding to 2I
!          configurations
!######################################################################
  subroutine avec_2I(n,kconf,ibconf2I,n_int_I,ksop_full,ndiff,Dw,&
       nbefore,socc,nsocc,knopen,knsp,averageii,confdim,cfg,&
       csfdim,Avec,vec0,nroots,refdim)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hbuild_mrci
    use mrci_integrals
    use dftmrci
    
    implicit none
    
    ! Bra configuration counter
    integer(is), intent(inout) :: ibconf2I

    ! Index of the 2-hole configuration
    integer(is), intent(in)    :: n

    ! Dimensions
    integer(is), intent(in)    :: n_int_I
    integer(is), intent(in)    :: csfdim,confdim,nroots,refdim

    ! Ket reference space configuration
    integer(is), intent(in)    :: kconf
    integer(ib), intent(in)    :: ksop_full(n_int,2)
    integer(is), intent(in)    :: knopen,knsp

    ! Difference configuration information
    integer(is),intent(in)     :: ndiff
    integer(is),intent(in)     :: Dw(nmo,2)
    
    ! Number of open shells preceding each MO
    integer(is), intent(in)    :: nbefore(nmo)

    ! Indices of the singly-occupied ket MOs
    integer(is), intent(in)    :: socc(nmo)
    integer(is), intent(in)    :: nsocc

    ! Spin-coupling averaged on-diagonal Hamiltonian
    ! matrix element values
    real(dp), intent(in)       :: averageii(confdim)

    ! MRCI configuration derived type
    type(mrcfg), intent(in)    :: cfg

    ! A-vector
    real(dp), intent(inout)    :: Avec(csfdim,nroots)

    ! Reference space eigenvectors
    real(dp), intent(in)       :: vec0(refdim,nroots)
    
    ! Working arrays
    integer(ib)                :: bconf_int(n_int_I,2)
    integer(ib)                :: bsop_int(n_int_I,2)
    integer(is), parameter     :: maxexci=2
    integer(is)                :: hlist(maxexci),plist(maxexci)

    ! Everything else
    integer(is)                :: ioff,nexci,bnopen,bnsp

!----------------------------------------------------------------------
! Compute the 2I elements
!----------------------------------------------------------------------
    ! Loop over 2I configurations generated by the current
    ! 2-hole configuration
    do ioff=cfg%off2I(n),cfg%off2I(n+1)-1

       ! Increment the 2I configuration counter
       ibconf2I=ibconf2I+1
    
       ! Bra 2I configuration in the internal MO space
       bconf_int=0_ib
       bconf_int=0_ib
       bconf_int=cfg%conf2I(1:n_int_I,:,ibconf2I)
       bsop_int=cfg%sop2I(1:n_int_I,:,ibconf2I)

       ! Compute the excitation degree between the two
       ! configurations
       nexci=exc_degree_conf(cfg%conf0h(:,:,kconf),bconf_int,n_int_I)
    
       ! Cycle if the excitation degree is greater than 2
       if (nexci > 2) cycle

       ! Number of open shells in the bra configuration
       bnopen=sop_nopen(bsop_int,n_int_I)
    
       ! Number of bra CSFs
       bnsp=ncsfs(bnopen)
    
       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(cfg%conf0h(:,:,kconf),bconf_int,&
            n_int_I,hlist(1:nexci),plist(1:nexci),nexci)

       ! Compute the matrix elements between the CSFs generated
       ! by the bra and ket configurations
       call hij_mrci(harr2,harr2dim,nexci,&
            ibconf2I,kconf,&
            cfg%sop2I(:,:,ibconf2I),ksop_full,&
            bnsp,knsp,bnopen,knopen,hlist,plist,cfg%m2c,&
            socc,nsocc,nbefore,Dw,ndiff,&
            cfg%csfs2I,cfg%csfs0h,cfg%n2I+1,cfg%n0h+1,&
            averageii(ibconf2I+cfg%n0h+cfg%n1I),&
            averageii(kconf))

       ! Contraction of the Hamiltonian matrix elements with
       ! the reference space eigenvectors
       call contract_hmat_vec0(ibconf2I,kconf,&
            cfg%csfs2I,cfg%csfs0h,cfg%n2I+1,cfg%n0h+1,&
            bnsp,knsp,Avec,vec0,csfdim,nroots,refdim)
       
    enddo
       
    return
    
  end subroutine avec_2I

!######################################################################
! avec_2E: calculation of the A-vector elements corresponding to 2E
!          configurations
!######################################################################
  subroutine avec_2E(n,kconf,ibconf2E,n_int_I,kconf_full,ksop_full,&
       ndiff,Dw,nbefore,socc,nsocc,knopen,knsp,averageii,confdim,cfg,&
       csfdim,Avec,vec0,nroots,refdim)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hbuild_mrci
    use mrci_integrals
    use dftmrci
    
    implicit none
    
    ! Bra configuration counter
    integer(is), intent(inout) :: ibconf2E

    ! Index of the 2-hole configuration
    integer(is), intent(in)    :: n

    ! Dimensions
    integer(is), intent(in)    :: n_int_I
    integer(is), intent(in)    :: csfdim,confdim,nroots,refdim

    ! Ket reference space configuration
    integer(is), intent(in)    :: kconf
    integer(ib), intent(in)    :: kconf_full(n_int,2)
    integer(ib), intent(in)    :: ksop_full(n_int,2)
    integer(is), intent(in)    :: knopen,knsp
    
    ! Difference configuration information
    integer(is),intent(in)    :: ndiff
    integer(is),intent(in)    :: Dw(nmo,2)
    
    ! Number of open shells preceding each MO
    integer(is), intent(in)    :: nbefore(nmo)

    ! Indices of the singly-occupied ket MOs
    integer(is), intent(in)    :: socc(nmo)
    integer(is), intent(in)    :: nsocc

    ! Spin-coupling averaged on-diagonal Hamiltonian
    ! matrix element values
    real(dp), intent(in)       :: averageii(confdim)

    ! MRCI configuration derived type
    type(mrcfg), intent(in)    :: cfg

    ! A-vector
    real(dp), intent(inout)    :: Avec(csfdim,nroots)

    ! Reference space eigenvectors
    real(dp), intent(in)       :: vec0(refdim,nroots)
    
    ! Working arrays
    integer(ib)                :: bconf_full(n_int,2)
    integer(ib)                :: bsop_full(n_int,2)
    integer(is), parameter     :: maxexci=2
    integer(is)                :: hlist(maxexci),plist(maxexci)

    ! Everything else
    integer(is)                :: ioff,nexci,bnopen,bnsp

!----------------------------------------------------------------------
! Compute the 2E elements
!----------------------------------------------------------------------
    ! Loop over 2E configurations generated by the current
    ! 2-hole configuration
    do ioff=cfg%off2E(n),cfg%off2E(n+1)-1

       ! Increment the 2E configuration counter
       ibconf2E=ibconf2E+1

       ! Bra 2E configuration in the full MO space
       bconf_full=0_ib
       bsop_full=0_ib
       bconf_full=cfg%conf2E(:,:,ibconf2E)
       bsop_full=cfg%sop2E(:,:,ibconf2E)

       ! Compute the excitation degree between the two
       ! configurations
       nexci=exc_degree_conf(kconf_full,bconf_full,n_int)
       
       ! Cycle if the excitation degree is greater than 2
       if (nexci > 2) cycle
       
       ! Number of open shells in the bra configuration
       bnopen=sop_nopen(bsop_full,n_int)
    
       ! Number of bra CSFs
       bnsp=ncsfs(bnopen)
    
       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(kconf_full,bconf_full,&
            n_int,hlist(1:nexci),plist(1:nexci),nexci)

       ! Compute the matrix elements between the CSFs generated
       ! by the bra and ket configurations
       call hij_mrci(harr2,harr2dim,nexci,&
            ibconf2E,kconf,&
            cfg%sop2E(:,:,ibconf2E),ksop_full,&
            bnsp,knsp,bnopen,knopen,hlist,plist,cfg%m2c,&
            socc,nsocc,nbefore,Dw,ndiff,&
            cfg%csfs2E,cfg%csfs0h,cfg%n2E+1,cfg%n0h+1,&
            averageii(ibconf2E+cfg%n0h+cfg%n1I+cfg%n2I+cfg%n1E),&
            averageii(kconf))

       ! Contraction of the Hamiltonian matrix elements with
       ! the reference space eigenvectors
       call contract_hmat_vec0(ibconf2E,kconf,&
            cfg%csfs2E,cfg%csfs0h,cfg%n2E+1,cfg%n0h+1,&
            bnsp,knsp,Avec,vec0,csfdim,nroots,refdim)
       
    enddo
       
    return
    
  end subroutine avec_2E

!######################################################################
! avec_1I1E: calculation of the A-vector elements corresponding to
!            1I1E configurations
!######################################################################
  subroutine avec_1I1E(n,kconf,ibconf1I1E,n_int_I,kconf_full,&
       ksop_full,ndiff,Dw,nbefore,socc,nsocc,knopen,knsp,averageii,&
       confdim,cfg,csfdim,Avec,vec0,nroots,refdim)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hbuild_mrci
    use mrci_integrals
    use dftmrci
    
    implicit none

    ! Bra configuration counter
    integer(is), intent(inout) :: ibconf1I1E

    ! Index of the 2-hole configuration
    integer(is), intent(in)    :: n

    ! Dimensions
    integer(is), intent(in)    :: n_int_I
    integer(is), intent(in)    :: csfdim,confdim,nroots,refdim

    ! Ket reference space configuration
    integer(is), intent(in)    :: kconf
    integer(ib), intent(in)    :: kconf_full(n_int,2)
    integer(ib), intent(in)    :: ksop_full(n_int,2)
    integer(is), intent(in)    :: knopen,knsp

    ! Difference configuration information
    integer(is),intent(in)    :: ndiff
    integer(is),intent(in)    :: Dw(nmo,2)
    
    ! Number of open shells preceding each MO
    integer(is), intent(in)    :: nbefore(nmo)

    ! Indices of the singly-occupied ket MOs
    integer(is), intent(in)    :: socc(nmo)
    integer(is), intent(in)    :: nsocc

    ! Spin-coupling averaged on-diagonal Hamiltonian
    ! matrix element values
    real(dp), intent(in)       :: averageii(confdim)

    ! MRCI configuration derived type
    type(mrcfg), intent(in)    :: cfg

    ! A-vector
    real(dp), intent(inout)    :: Avec(csfdim,nroots)

    ! Reference space eigenvectors
    real(dp), intent(in)       :: vec0(refdim,nroots)
    
    ! Working arrays
    integer(ib)                :: bconf_full(n_int,2)
    integer(ib)                :: bsop_full(n_int,2)
    integer(is), parameter     :: maxexci=2
    integer(is)                :: hlist(maxexci),plist(maxexci)

    ! Everything else
    integer(is)                :: ioff,nexci,bnopen,bnsp

!----------------------------------------------------------------------
! Compute the 1I1E elements
!----------------------------------------------------------------------
    ! Loop over 1I1E configurations generated by the current
    ! 2-hole configuration
    do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1

       ! Increment the 1I1E configuration counter
       ibconf1I1E=ibconf1I1E+1
    
       ! Bra 1I1E configuration in the full MO space
       bconf_full=0_ib
       bsop_full=0_ib
       bconf_full=cfg%conf1I1E(:,:,ibconf1I1E)
       bsop_full=cfg%sop1I1E(:,:,ibconf1I1E)

       ! Compute the excitation degree between the two
       ! configurations
       nexci=exc_degree_conf(kconf_full,bconf_full,n_int)
       
       ! Cycle if the excitation degree is greater than 2
       if (nexci > 2) cycle       

       ! Number of open shells in the bra configuration
       bnopen=sop_nopen(bsop_full,n_int)
    
       ! Number of bra CSFs
       bnsp=ncsfs(bnopen)
    
       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(kconf_full,bconf_full,&
            n_int,hlist(1:nexci),plist(1:nexci),nexci)

       ! Compute the matrix elements between the CSFs generated
       ! by the bra and ket configurations
       call hij_mrci(harr2,harr2dim,nexci,&
            ibconf1I1E,kconf,&
            cfg%sop1I1E(:,:,ibconf1I1E),ksop_full,&
            bnsp,knsp,bnopen,knopen,hlist,plist,cfg%m2c,&
            socc,nsocc,nbefore,Dw,ndiff,&
            cfg%csfs1I1E,cfg%csfs0h,cfg%n1I1E+1,cfg%n0h+1,&
            averageii(ibconf1I1E+cfg%n0h+cfg%n1I+cfg%n2I+cfg%n1E+cfg%n2E),&
            averageii(kconf))

       ! Contraction of the Hamiltonian matrix elements with
       ! the reference space eigenvectors
       call contract_hmat_vec0(ibconf1I1E,kconf,&
            cfg%csfs1I1E,cfg%csfs0h,cfg%n1I1E+1,cfg%n0h+1,&
            bnsp,knsp,Avec,vec0,csfdim,nroots,refdim)
       
    enddo
       
    return
    
  end subroutine avec_1I1E
  
!######################################################################
! contract_hmat_vec0: contracts a block of the Hamiltonian matrix with
!                     the corresponding part of the reference space
!                     eigenvectors
!######################################################################
  subroutine contract_hmat_vec0(bconf,kconf,bcsfs,kcsfs,bdim,kdim,&
       bnsp,knsp,Avec,vec0,csfdim,nroots,refdim)

    use constants
    use bitglobal

    implicit none

    ! Configuration indices
    integer(is), intent(in) :: bconf,kconf
    
    ! CSF offsets
    integer(is), intent(in) :: kdim,bdim
    integer(is), intent(in) :: bcsfs(bdim),kcsfs(kdim)

    ! Numbers of bra and ket CSFs
    integer(is), intent(in) :: bnsp,knsp

    ! A-vector
    integer(is), intent(in) :: csfdim,nroots
    real(dp), intent(inout) :: Avec(csfdim,nroots)

    ! Reference space eigenvectors
    integer(is), intent(in) :: refdim
    real(dp), intent(in)    :: vec0(refdim,nroots)
    
    ! Everything else
    integer(is)             :: ikcsf,ibcsf
    integer(is)             :: j,counter

    ! Loop over roots
    do j=1,nroots
    
       counter=0
       
       ! Loop over ket CSFs (reference space CSFs)
       do ikcsf=kcsfs(kconf),kcsfs(kconf+1)-1
          
          ! Loop over bra CSFs (MRCI CSFs not in the reference space)
          do ibcsf=bcsfs(bconf),bcsfs(bconf+1)-1
             counter=counter+1
             
             Avec(ibcsf,j)=Avec(ibcsf,j)&
                  +harr2(counter)*vec0(ikcsf,j)
             
          enddo
                    
       enddo
    
    enddo

    return
    
  end subroutine contract_hmat_vec0
    
!######################################################################
! apply_denominator: Application of the (H_nn - E^0_I) to get
!                    the A-vector element values
!######################################################################
  subroutine apply_denominator(cfg,Avec,E2,hdiag,e0,csfdim,nroots,&
       refdim,shift,idsp)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! Dimensions
    integer(is), intent(in)  :: csfdim,nroots,refdim

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: hdiag(csfdim)

    ! ENPT2 wave function and energy corrections
    real(dp), intent(inout)  :: Avec(csfdim,nroots)
    real(dp), intent(out)    :: E2(nroots)
    
    ! Reference space eigenvalues
    real(dp), intent(in)     :: e0(nroots)

    ! ISA shift
    real(dp), intent(in)     :: shift

    ! Damped strong perturbers
    integer(is), intent(out) :: idsp(csfdim)    
    real(dp), parameter      :: cthrsh=0.055d0
    
    ! Everything else
    integer(is)              :: j,icsf
    real(dp)                 :: dj,Aold,ediff
    
!----------------------------------------------------------------------
! ENPT2 energy and wave function corrections
!----------------------------------------------------------------------
    ! Initialisation
    E2=0.0d0
    idsp=0
    
    ! Loop over roots
    do j=1,nroots

       ! Loop over CSFs (excluding the reference space ones)
       do icsf=refdim+1,csfdim

          ! E^(0) - H_ii
          ediff=e0(j)-hdiag(icsf)
          
          ! ISA factor
          dj=shift/ediff

          ! Unshifted A-vector element
          Aold=Avec(icsf,j)/ediff
          
          ! Energy correction
          E2(j)=E2(j)+Avec(icsf,j)**2/(ediff+dj)
          
          ! A-vector element
          Avec(icsf,j)=Avec(icsf,j)/(ediff+dj)

          ! Make sure that all strong perturbers are captured
          if (abs(Aold) >= cthrsh &
               .and. abs(Avec(icsf,j)) < cthrsh) idsp(icsf)=1
          
       enddo
       
    enddo

    return
    
  end subroutine apply_denominator

!######################################################################
! qdpt2: constructs and diagonalises the QDPT2 effective Hamiltonian
!        using the ENPT2 Hamiltonian partitioning  
!######################################################################
  subroutine qdpt2(cfg,Avec,hdiag,e0,EQD,mix,csfdim,nroots,refdim,&
       shift)

    use constants
    use bitglobal
    use conftype
    use utils
    use iomod
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Dimensions
    integer(is), intent(in) :: csfdim,nroots,refdim

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: hdiag(csfdim)

    ! A-vector
    real(dp), intent(inout) :: Avec(csfdim,nroots)
    
    ! Reference space eigenvalues
    real(dp), intent(in)    :: e0(nroots)

    ! QDPT2 energies and mixing coefficients
    real(dp), intent(out)   :: EQD(nroots)
    real(dp), intent(out)   :: mix(nroots,nroots)

    ! ISA shift
    real(dp), intent(in)    :: shift
    
    ! Effective Hamiltonian matrix and eigenpairs
    real(dp), allocatable   :: heff(:,:)

    ! Everything else
    integer(is)             :: i,j,icsf
    real(dp)                :: fac,di,dj
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(heff(nroots,nroots))
    heff=0.0d0

!----------------------------------------------------------------------
! Construct the QDPT2 effective Hamiltonian matrix
!----------------------------------------------------------------------
    ! Loop over pairs of roots
    do i=1,nroots
       do j=i,nroots

          ! E_i^(0) \delta_ij
          if (i == j) heff(i,j)=e0(i)

          ! Loop over FOIS CSFs
          do icsf=refdim+1,csfdim

             ! ISA factors
             di=shift/(e0(i)-hdiag(icsf))
             dj=shift/(e0(j)-hdiag(icsf))

             ! ISA-shifted denominators
             fac=0.5d0/(e0(i)-hdiag(icsf)+di)
             fac=fac+0.5d0/(e0(j)-hdiag(icsf)+dj)

             ! <I|H|i> <i|H|J> / (E_I^0 - H_ii + Delta_Ii) + (I<->J)
             heff(i,j)=heff(i,j)+Avec(icsf,i)*Avec(icsf,j)*fac
             
          enddo
          
          heff(j,i)=heff(i,j)

       enddo
    enddo
    
!----------------------------------------------------------------------
! Diagonalise the effective Hamiltonian
!----------------------------------------------------------------------
    ! Eigenpairs of H_eff
    call diag_matrix_real(heff,EQD,mix,nroots)

    ! Add E_SCF to the eigenvalues
    EQD=EQD+escf
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(heff)
    
    return
    
  end subroutine qdpt2
  
!######################################################################
  
end module epstein_nesbet
