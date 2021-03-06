!**********************************************************************
! Calculation of the on-diagonal matrix elements of the MRCI
! Hamiltonian matrix
!**********************************************************************
module hii

  implicit none

contains

!######################################################################
! hmat_diagonal: computes the on-diagonal Hamiltonian matrix
!                elements plus their spin-coupling-averaged values
!######################################################################
  subroutine hmat_diagonal(hdiag,csfdim,averageii,confdim,cfg)
    
    use constants
    use bitglobal
    use mrciutils
    use bitutils
    use conftype
    use hparam
    use hbuild_mrci
    use mrci_integrals
    use dftmrci
    
    implicit none

    ! Total number of CSFs and configurations
    integer(is), intent(in) :: csfdim,confdim
    
    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(out)   :: hdiag(csfdim)

    ! On-diagonal Hamiltonian matrix elements
    ! averaged over spin couplings
    real(dp), intent(out)   :: averageii(confdim)

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg
    
    ! Working arrays
    integer(ib)             :: conf_full(n_int,2)
    integer(ib)             :: sop_full(n_int,2)
    real(dp), allocatable   :: harr(:)
    
    ! MO classes
    integer(is)             :: socc(nmo),docc(nmo),unocc(nmo)
    integer(is)             :: nopen,nsocc,ndocc,nunocc

    ! Difference configuration information
    integer(is)             :: ndiff
    integer(is)             :: Dw(nmo,2)

    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)
    
    ! Everything else
    integer(is)             :: iconf,iomega,icsf,nsp
    integer(is)             :: n,ioff,counter,ilast
    integer(is)             :: n_int_I
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(harr(maxval(ncsfs(0:nomax))))
    harr=0.0d0
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I
    hdiag=0.0d0
    averageii=0.0d0
    
!----------------------------------------------------------------------
! (1) Reference space configurations
!----------------------------------------------------------------------
    ! Loop over reference configurations
    do iconf=1,cfg%n0h

       ! Number of open shells
       nopen=sop_nopen(cfg%sop0h(:,:,iconf),n_int_I)

       ! Configuration and SOP in the full MO space
       conf_full=0_ib
       sop_full=0_ib
       conf_full(1:n_int_I,:)=cfg%conf0h(:,:,iconf)
       sop_full(1:n_int_I,:)=cfg%sop0h(:,:,iconf)

       ! Get all the configuration information needed to evaluate
       ! the on-diagonal Hamiltonian matrix element
       call package_confinfo(sop_full,conf_full,unocc,socc,docc,nunocc,&
            nsocc,ndocc,Dw,ndiff,nbefore)

       ! Number of CSFs generated by the configuration
       nsp=ncsfs(nopen)

       ! On-diagonal Hamiltonian matrix elements for all CSFs generated
       ! by the current configuration
       call hii_mrci(harr(1:nsp),nsp,sop_full,nopen,socc,nsocc,&
            docc,ndocc,unocc,nunocc,nbefore,Dw,ndiff,cfg%m2c)

       ! Apply any DFT/MRCI corrections to the on-diagonal
       ! matrix elements
       if (ldftmrci) call hii_dftmrci(harr(1:nsp),nsp,Dw,ndiff,&
            nopen,cfg%m2c,sop_full,socc,nsocc,nbefore)

       ! Fill in the hdiag and averageii arrays
       iomega=0
       do icsf=cfg%csfs0h(iconf),cfg%csfs0h(iconf+1)-1
          iomega=iomega+1
          hdiag(icsf)=harr(iomega)
          averageii(iconf)=averageii(iconf)+harr(iomega)
       enddo
       averageii(iconf)=averageii(iconf)/nsp

    enddo
    
!----------------------------------------------------------------------
! (2) 1I configurations
!----------------------------------------------------------------------
    ilast=cfg%n0h
    
    if (cfg%n1I > 0) then

       ! Initialise the configuration counter
       counter=0
       
       ! Loop over 1-hole configurations
       do n=1,cfg%n1h

          ! Loop over the 1I configurations generated by the
          ! current 1-hole configuration
          do ioff=cfg%off1I(n),cfg%off1I(n+1)-1

             ! Increment the configuration counter
             counter=counter+1

             ! Configuration and SOP
             conf_full=cfg%conf1I(:,:,counter)
             sop_full=cfg%sop1I(:,:,counter)

             ! Number of open shells
             nopen=sop_nopen(sop_full,n_int)
             
             ! Get all the configuration information needed to evaluate
             ! the on-diagonal Hamiltonian matrix element
             call package_confinfo(sop_full,conf_full,unocc,socc,docc,&
                  nunocc,nsocc,ndocc,Dw,ndiff,nbefore)
             
             ! Number of CSFs generated by the configuration
             nsp=ncsfs(nopen)
             
             ! On-diagonal Hamiltonian matrix elements for all CSFs 
             ! generated by the current configuration
             call hii_mrci(harr(1:nsp),nsp,sop_full,nopen,socc,nsocc,&
                  docc,ndocc,unocc,nunocc,nbefore,Dw,ndiff,cfg%m2c)
             
             ! Apply any DFT/MRCI corrections to the on-diagonal
             ! matrix elements
             if (ldftmrci) call hii_dftmrci(harr(1:nsp),nsp,Dw,ndiff,&
                  nopen,cfg%m2c,sop_full,socc,nsocc,nbefore)

             ! Fill in the hdiag and averageii arrays
             iomega=0
             do icsf=cfg%csfs1I(counter),cfg%csfs1I(counter+1)-1
                iomega=iomega+1
                hdiag(icsf)=harr(iomega)
                averageii(counter+ilast)=&
                     averageii(counter+ilast)+harr(iomega)
             enddo
             averageii(counter+ilast)=averageii(counter+ilast)/nsp
             
          enddo
          
       enddo
       
    endif

!----------------------------------------------------------------------
! (3) 2I configurations
!----------------------------------------------------------------------
    ilast=cfg%n0h+cfg%n1I
    
    if (cfg%n2I > 0) then

       ! Initialise the configuration counter
       counter=0
       
       ! Loop over 2-hole configurations
       do n=1,cfg%n2h

          ! Loop over the 2I configurations generated by the
          ! current 2-hole configuration
          do ioff=cfg%off2I(n),cfg%off2I(n+1)-1

             ! Increment the configuration counter
             counter=counter+1

             ! Configuration and SOP
             conf_full=cfg%conf2I(:,:,counter)
             sop_full=cfg%sop2I(:,:,counter)

             ! Number of open shells
             nopen=sop_nopen(sop_full,n_int)

             ! Get all the configuration information needed to evaluate
             ! the on-diagonal Hamiltonian matrix element
             call package_confinfo(sop_full,conf_full,unocc,socc,docc,&
                  nunocc,nsocc,ndocc,Dw,ndiff,nbefore)
             
             ! Number of CSFs generated by the configuration
             nsp=ncsfs(nopen)
             
             ! On-diagonal Hamiltonian matrix elements for all CSFs 
             ! generated by the current configuration
             call hii_mrci(harr(1:nsp),nsp,sop_full,nopen,socc,nsocc,&
                  docc,ndocc,unocc,nunocc,nbefore,Dw,ndiff,cfg%m2c)
             
             ! Apply any DFT/MRCI corrections to the on-diagonal
             ! matrix elements
             if (ldftmrci) call hii_dftmrci(harr(1:nsp),nsp,Dw,ndiff,&
                  nopen,cfg%m2c,sop_full,socc,nsocc,nbefore)

             ! Fill in the hdiag and averageii arrays
             iomega=0
             do icsf=cfg%csfs2I(counter),cfg%csfs2I(counter+1)-1
                iomega=iomega+1
                hdiag(icsf)=harr(iomega)
                averageii(counter+ilast)=&
                     averageii(counter+ilast)+harr(iomega)
             enddo
             averageii(counter+ilast)=averageii(counter+ilast)/nsp
             
          enddo
             
       enddo
          
    endif

!----------------------------------------------------------------------
! (4) 1E configurations
!----------------------------------------------------------------------
    ilast=cfg%n0h+cfg%n1I+cfg%n2I
    
    if (cfg%n1E > 0) then

       ! Initialise the configuration counter
       counter=0
       
       ! Loop over 1-hole configurations
       do n=1,cfg%n1h

          ! Loop over the 1E configurations generated by the
          ! current 1-hole configuration
          do ioff=cfg%off1E(n),cfg%off1E(n+1)-1

             ! Increment the configuration counter
             counter=counter+1

             ! Configuration and SOP
             conf_full=cfg%conf1E(:,:,counter)
             sop_full=cfg%sop1E(:,:,counter)

             ! Number of open shells
             nopen=sop_nopen(sop_full,n_int)
             
             ! Get all the configuration information needed to evaluate
             ! the on-diagonal Hamiltonian matrix element
             call package_confinfo(sop_full,conf_full,unocc,socc,docc,&
                  nunocc,nsocc,ndocc,Dw,ndiff,nbefore)
             
             ! Number of CSFs generated by the configuration
             nsp=ncsfs(nopen)
             
             ! On-diagonal Hamiltonian matrix elements for all CSFs 
             ! generated by the current configuration
             call hii_mrci(harr(1:nsp),nsp,sop_full,nopen,socc,nsocc,&
                  docc,ndocc,unocc,nunocc,nbefore,Dw,ndiff,cfg%m2c)
             
             ! Apply any DFT/MRCI corrections to the on-diagonal
             ! matrix elements
             if (ldftmrci) call hii_dftmrci(harr(1:nsp),nsp,Dw,ndiff,&
                  nopen,cfg%m2c,sop_full,socc,nsocc,nbefore)

             ! Fill in the hdiag and averageii arrays
             iomega=0
             do icsf=cfg%csfs1E(counter),cfg%csfs1E(counter+1)-1
                iomega=iomega+1
                hdiag(icsf)=harr(iomega)
                averageii(counter+ilast)=&
                     averageii(counter+ilast)+harr(iomega)
             enddo
             averageii(counter+ilast)=averageii(counter+ilast)/nsp
             
          enddo
             
       enddo
          
    endif

!----------------------------------------------------------------------
! (5) 2E configurations
!----------------------------------------------------------------------
    ilast=cfg%n0h+cfg%n1I+cfg%n2I+cfg%n1E

    if (cfg%n2E > 0) then

       ! Initialise the configuration counter
       counter=0
       
       ! Loop over 2-hole configurations
       do n=1,cfg%n2h

          ! Loop over the 2E configurations generated by the
          ! current 2-hole configuration
          do ioff=cfg%off2E(n),cfg%off2E(n+1)-1

             ! Increment the configuration counter
             counter=counter+1

             ! Configuration and SOP
             conf_full=cfg%conf2E(:,:,counter)
             sop_full=cfg%sop2E(:,:,counter)

             ! Number of open shells
             nopen=sop_nopen(sop_full,n_int)

             ! Get all the configuration information needed to evaluate
             ! the on-diagonal Hamiltonian matrix element
             call package_confinfo(sop_full,conf_full,unocc,socc,docc,&
                  nunocc,nsocc,ndocc,Dw,ndiff,nbefore)
             
             ! Number of CSFs generated by the configuration
             nsp=ncsfs(nopen)
             
             ! On-diagonal Hamiltonian matrix elements for all CSFs 
             ! generated by the current configuration
             call hii_mrci(harr(1:nsp),nsp,sop_full,nopen,socc,nsocc,&
                  docc,ndocc,unocc,nunocc,nbefore,Dw,ndiff,cfg%m2c)
             
             ! Apply any DFT/MRCI corrections to the on-diagonal
             ! matrix elements
             if (ldftmrci) call hii_dftmrci(harr(1:nsp),nsp,Dw,ndiff,&
                  nopen,cfg%m2c,sop_full,socc,nsocc,nbefore)

             ! Fill in the hdiag and averageii arrays
             iomega=0
             do icsf=cfg%csfs2E(counter),cfg%csfs2E(counter+1)-1
                iomega=iomega+1
                hdiag(icsf)=harr(iomega)
                averageii(counter+ilast)=&
                     averageii(counter+ilast)+harr(iomega)
             enddo
             averageii(counter+ilast)=averageii(counter+ilast)/nsp
             
          enddo

       enddo
             
    endif

!----------------------------------------------------------------------
! (6) 1I1E configurations
!----------------------------------------------------------------------
    ilast=cfg%n0h+cfg%n1I+cfg%n2I+cfg%n1E+cfg%n2E
    
    if (cfg%n1I1E > 0) then

       ! Initialise the configuration counter
       counter=0
       
       ! Loop over 2-hole configurations
       do n=1,cfg%n2h

          ! Loop over the 1I1E configurations generated by
          ! the current 2-hole configuration
          do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1

             ! Increment the configuration counter
             counter=counter+1

             ! Configuration and SOP
             conf_full=cfg%conf1I1E(:,:,counter)
             sop_full=cfg%sop1I1E(:,:,counter)

             ! Number of open shells
             nopen=sop_nopen(sop_full,n_int)

             ! Get all the configuration information needed to evaluate
             ! the on-diagonal Hamiltonian matrix element
             call package_confinfo(sop_full,conf_full,unocc,socc,docc,&
                  nunocc,nsocc,ndocc,Dw,ndiff,nbefore)
             
             ! Number of CSFs generated by the configuration
             nsp=ncsfs(nopen)
             
             ! On-diagonal Hamiltonian matrix elements for all CSFs 
             ! generated by the current configuration
             call hii_mrci(harr(1:nsp),nsp,sop_full,nopen,socc,nsocc,&
                  docc,ndocc,unocc,nunocc,nbefore,Dw,ndiff,cfg%m2c)
             
             ! Apply any DFT/MRCI corrections to the on-diagonal
             ! matrix elements
             if (ldftmrci) call hii_dftmrci(harr(1:nsp),nsp,Dw,ndiff,&
                  nopen,cfg%m2c,sop_full,socc,nsocc,nbefore)

             ! Fill in the hdiag and averageii arrays
             iomega=0
             do icsf=cfg%csfs1I1E(counter),cfg%csfs1I1E(counter+1)-1
                iomega=iomega+1
                hdiag(icsf)=harr(iomega)
                averageii(counter+ilast)=&
                     averageii(counter+ilast)+harr(iomega)
             enddo
             averageii(counter+ilast)=averageii(counter+ilast)/nsp
             
          enddo

       enddo
          
    endif

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(harr)
    
    return
    
  end subroutine hmat_diagonal
  
!######################################################################
  
end module hii
