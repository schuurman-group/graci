!***********************************************************************
! Routines common to two or more perturbation theory methods
!***********************************************************************
module pt2_common

  implicit none

contains

!######################################################################
! avec_1h: Calculation of the elements of the A-vector corresponding
!          to the CSFs generated by the 1-hole configurations
!######################################################################
  subroutine avec_1h(cfg,Avec,averageii,vec0,csfdim,confdim,refdim,&
       nroots,harr2,harr2dim)

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

    ! Working array
    integer(is), intent(in) :: harr2dim
    real(dp), intent(inout) :: harr2(harr2dim)
    
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
                     averageii,csfdim,confdim,Avec,vec0,nroots,refdim,&
                     harr2,harr2dim)
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
                     csfdim,Avec,vec0,nroots,refdim,&
                     harr2,harr2dim)
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
       nroots,harr2,harr2dim)

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

    ! Working array
    integer(is), intent(in) :: harr2dim
    real(dp), intent(inout) :: harr2(harr2dim)
    
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
                     nroots,refdim,harr2,harr2dim)
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
                     Avec,vec0,nroots,refdim,harr2,harr2dim)
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
                     Avec,vec0,nroots,refdim,harr2,harr2dim)
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
       confdim,Avec,vec0,nroots,refdim,harr2,harr2dim)

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

    ! Working array
    integer(is), intent(in)    :: harr2dim
    real(dp), intent(inout)    :: harr2(harr2dim)
    
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
            bnsp,knsp,Avec,vec0,csfdim,nroots,refdim,&
            harr2,harr2dim)
       
    enddo
       
    return
    
  end subroutine avec_1I

!######################################################################
! avec_1E: calculation of the A-vector elements corresponding to
!          1E configurations
!######################################################################
  subroutine avec_1E(n,kconf,ibconf1E,n_int_I,kconf_full,&
       ksop_full,ndiff,Dw,nbefore,socc,nsocc,knopen,knsp,averageii,&
       confdim,cfg,csfdim,Avec,vec0,nroots,refdim,harr2,harr2dim)

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

    ! Working array
    integer(is), intent(in)    :: harr2dim
    real(dp), intent(inout)    :: harr2(harr2dim)
    
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
            bnsp,knsp,Avec,vec0,csfdim,nroots,refdim,&
            harr2,harr2dim)
       
    enddo

    return

  end subroutine avec_1E

!######################################################################
! avec_2I: calculation of the A-vector elements corresponding to 2I
!          configurations
!######################################################################
  subroutine avec_2I(n,kconf,ibconf2I,n_int_I,ksop_full,ndiff,Dw,&
       nbefore,socc,nsocc,knopen,knsp,averageii,confdim,cfg,&
       csfdim,Avec,vec0,nroots,refdim,harr2,harr2dim)

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

    ! Working array
    integer(is), intent(in)    :: harr2dim
    real(dp), intent(inout)    :: harr2(harr2dim)
    
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
            bnsp,knsp,Avec,vec0,csfdim,nroots,refdim,&
            harr2,harr2dim)
       
    enddo
       
    return
    
  end subroutine avec_2I

!######################################################################
! avec_2E: calculation of the A-vector elements corresponding to 2E
!          configurations
!######################################################################
  subroutine avec_2E(n,kconf,ibconf2E,n_int_I,kconf_full,ksop_full,&
       ndiff,Dw,nbefore,socc,nsocc,knopen,knsp,averageii,confdim,cfg,&
       csfdim,Avec,vec0,nroots,refdim,harr2,harr2dim)

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

    ! Working array
    integer(is), intent(in)    :: harr2dim
    real(dp), intent(inout)    :: harr2(harr2dim)
    
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
            bnsp,knsp,Avec,vec0,csfdim,nroots,refdim,&
            harr2,harr2dim)
       
    enddo
       
    return
    
  end subroutine avec_2E

!######################################################################
! avec_1I1E: calculation of the A-vector elements corresponding to
!            1I1E configurations
!######################################################################
  subroutine avec_1I1E(n,kconf,ibconf1I1E,n_int_I,kconf_full,&
       ksop_full,ndiff,Dw,nbefore,socc,nsocc,knopen,knsp,averageii,&
       confdim,cfg,csfdim,Avec,vec0,nroots,refdim,harr2,harr2dim)

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

    ! Working array
    integer(is), intent(in)    :: harr2dim
    real(dp), intent(inout)    :: harr2(harr2dim)
    
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
            bnsp,knsp,Avec,vec0,csfdim,nroots,refdim,&
            harr2,harr2dim)
       
    enddo
       
    return
    
  end subroutine avec_1I1E
  
!######################################################################
! contract_hmat_vec0: contracts a block of the Hamiltonian matrix with
!                     the corresponding part of the reference space
!                     eigenvectors
!######################################################################
  subroutine contract_hmat_vec0(bconf,kconf,bcsfs,kcsfs,bdim,kdim,&
       bnsp,knsp,Avec,vec0,csfdim,nroots,refdim,harr2,harr2dim)

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

    ! Working array
    integer(is), intent(in) :: harr2dim
    real(dp), intent(inout) :: harr2(harr2dim)
    
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
  
end module pt2_common
