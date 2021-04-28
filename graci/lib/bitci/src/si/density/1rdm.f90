!**********************************************************************
! Routines for the calculation of 1-RDMs for MRCI wavefunctions
!**********************************************************************
module rdm

  use constants
  
  implicit none

  ! Spin-coupling coefficients
  real(dp), allocatable, private :: spincp(:,:)
  
contains

!######################################################################
! 1rdm_mrci: Master routine for the calculation of MRCI 1-RDMs for
!            a specified set of states of a single irrep
!######################################################################
  subroutine rdm_mrci(cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use iomod
    use timing
    
    implicit none

    ! No. CSFs
    integer(is), intent(in) :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in) :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(out)   :: dmat(nmo,nmo,nroots)
   
    ! Timing variables
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                                twall_end 
    
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)  
    
!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    dmat=0.0d0

    allocate(spincp(ncsfs(nomax),ncsfs(nomax)))
    spincp=0.0d0

!----------------------------------------------------------------------
! (1) On-diagonal elements
!----------------------------------------------------------------------
    call rdm_diag(cfg,csfdim,nroots,vec,dmat)
    
!----------------------------------------------------------------------
! (2) Ref-Ref contributions to the 1-RDMs
!----------------------------------------------------------------------
    call rdm_0h_0h(cfg,csfdim,nroots,vec,dmat)

!----------------------------------------------------------------------
! (3) Ref - 1H contributions to the 1-RDMs
!----------------------------------------------------------------------
    call rdm_0h_1h(cfg,csfdim,nroots,vec,dmat)

!----------------------------------------------------------------------
! (4) Ref - 2H contributions to the 1-RDMs
!----------------------------------------------------------------------
    call rdm_0h_2h(cfg,csfdim,nroots,vec,dmat)

!----------------------------------------------------------------------
! (5) 1H - 1H contributions to the 1-RDMs
!----------------------------------------------------------------------
    call rdm_1h_1h(cfg,csfdim,nroots,vec,dmat)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(spincp)
    
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'rdm_mrci')
    
    return
    
  end subroutine rdm_mrci

!######################################################################
! rdm_mrci_ondiag: Calculation of the on-diagonal elements of the MRCI
!                  1-RDMs
!######################################################################
  subroutine rdm_diag(cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use iomod
    
    implicit none

    ! No. CSFs
    integer(is), intent(in) :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in) :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)

    ! Working arrays
    integer(ib)             :: conf_full(n_int,2)
    integer(ib)             :: sop_full(n_int,2)
        
    ! MO classes
    integer(is)             :: socc(nmo),docc(nmo),unocc(nmo)
    integer(is)             :: nsocc,ndocc,nunocc
    
    ! Everything else
    integer(is)             :: iconf,iomega,icsf
    integer(is)             :: i,n,ista,imo,ioff
    integer(is)             :: n_int_I
    real(dp)                :: c2,trace

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ! N_int^I
    n_int_I=cfg%n_int_I
    
!----------------------------------------------------------------------
! Reference space configurations
!----------------------------------------------------------------------
    ! Loop over reference configurationd
    do iconf=1,cfg%n0h

       ! Get the lists of singly-occupied and doubly-occupied MOs
       call sop_socc_list(cfg%sop0h(:,:,iconf),n_int_I,socc,nmo,nsocc)
       call sop_docc_list(cfg%sop0h(:,:,iconf),n_int_I,docc,nmo,ndocc)

       ! Loop over roots
       do ista=1,nroots

          ! Loop over CSFs generated by this configuration
          do icsf=cfg%csfs0h(iconf),cfg%csfs0h(iconf+1)-1

             ! Coefficient squared
             c2=vec(icsf,ista)**2

             ! Singly-occupied MO contributions to the 1-RDM
             do i=1,nsocc
                imo=socc(i)
                dmat(imo,imo,ista)=dmat(imo,imo,ista)+c2
             enddo

             ! Doubly-occupied MO contributions to the 1-RDM
             do i=1,ndocc
                imo=docc(i)
                dmat(imo,imo,ista)=dmat(imo,imo,ista)+2.0d0*c2
             enddo
             
          enddo
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
    if (cfg%n1I > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h

          ! Loop over the 1I configurations generated by the
          ! current 1-hole configuration
          do ioff=cfg%off1I(n),cfg%off1I(n+1)-1

             ! Get the lists of singly-occupied and doubly-occupied MOs
             call sop_socc_list(cfg%sop1I(:,:,ioff),n_int,socc,&
                  nmo,nsocc)
             call sop_docc_list(cfg%sop1I(:,:,ioff),n_int,docc,&
                  nmo,ndocc)

             ! Loop over roots
             do ista=1,nroots
                
                ! Loop over CSFs generated by this configuration
                do icsf=cfg%csfs1I(ioff),cfg%csfs1I(ioff+1)-1
                   
                   ! Coefficient squared
                   c2=vec(icsf,ista)**2
                   
                   ! Singly-occupied MO contributions to the 1-RDM
                   do i=1,nsocc
                      imo=socc(i)
                      dmat(imo,imo,ista)=dmat(imo,imo,ista)+c2
                   enddo
                   
                   ! Doubly-occupied MO contributions to the 1-RDM
                   do i=1,ndocc
                      imo=docc(i)
                      dmat(imo,imo,ista)=dmat(imo,imo,ista)+2.0d0*c2
                   enddo
                   
                enddo
          
             enddo
             
          enddo
             
       enddo

    endif

!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
    if (cfg%n2I > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h

          ! Loop over the 2I configurations generated by the
          ! current 2-hole configuration
          do ioff=cfg%off2I(n),cfg%off2I(n+1)-1

             ! Get the lists of singly-occupied and doubly-occupied MOs
             call sop_socc_list(cfg%sop2I(:,:,ioff),n_int,socc,&
                  nmo,nsocc)
             call sop_docc_list(cfg%sop2I(:,:,ioff),n_int,docc,&
                  nmo,ndocc)

             ! Loop over roots
             do ista=1,nroots
                
                ! Loop over CSFs generated by this configuration
                do icsf=cfg%csfs2I(ioff),cfg%csfs2I(ioff+1)-1
                   
                   ! Coefficient squared
                   c2=vec(icsf,ista)**2
                   
                   ! Singly-occupied MO contributions to the 1-RDM
                   do i=1,nsocc
                      imo=socc(i)
                      dmat(imo,imo,ista)=dmat(imo,imo,ista)+c2
                   enddo
                   
                   ! Doubly-occupied MO contributions to the 1-RDM
                   do i=1,ndocc
                      imo=docc(i)
                      dmat(imo,imo,ista)=dmat(imo,imo,ista)+2.0d0*c2
                   enddo
                   
                enddo
          
             enddo
             
          enddo
             
       enddo

    endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
    if (cfg%n1E > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h

          ! Loop over the 1E configurations generated by the
          ! current 1-hole configuration
          do ioff=cfg%off1E(n),cfg%off1E(n+1)-1

             ! Get the lists of singly-occupied and doubly-occupied MOs
             call sop_socc_list(cfg%sop1E(:,:,ioff),n_int,socc,&
                  nmo,nsocc)
             call sop_docc_list(cfg%sop1E(:,:,ioff),n_int,docc,&
                  nmo,ndocc)

             ! Loop over roots
             do ista=1,nroots
                
                ! Loop over CSFs generated by this configuration
                do icsf=cfg%csfs1E(ioff),cfg%csfs1E(ioff+1)-1
                   
                   ! Coefficient squared
                   c2=vec(icsf,ista)**2
                   
                   ! Singly-occupied MO contributions to the 1-RDM
                   do i=1,nsocc
                      imo=socc(i)
                      dmat(imo,imo,ista)=dmat(imo,imo,ista)+c2
                   enddo
                   
                   ! Doubly-occupied MO contributions to the 1-RDM
                   do i=1,ndocc
                      imo=docc(i)
                      dmat(imo,imo,ista)=dmat(imo,imo,ista)+2.0d0*c2
                   enddo
                   
                enddo
          
             enddo
             
          enddo
             
       enddo

    endif

!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
    if (cfg%n2E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h

          ! Loop over the 2E configurations generated by the
          ! current 2-hole configuration
          do ioff=cfg%off2E(n),cfg%off2E(n+1)-1

             ! Get the lists of singly-occupied and doubly-occupied MOs
             call sop_socc_list(cfg%sop2E(:,:,ioff),n_int,socc,&
                  nmo,nsocc)
             call sop_docc_list(cfg%sop2E(:,:,ioff),n_int,docc,&
                  nmo,ndocc)

             ! Loop over roots
             do ista=1,nroots
                
                ! Loop over CSFs generated by this configuration
                do icsf=cfg%csfs2E(ioff),cfg%csfs2E(ioff+1)-1
                   
                   ! Coefficient squared
                   c2=vec(icsf,ista)**2
                   
                   ! Singly-occupied MO contributions to the 1-RDM
                   do i=1,nsocc
                      imo=socc(i)
                      dmat(imo,imo,ista)=dmat(imo,imo,ista)+c2
                   enddo
                   
                   ! Doubly-occupied MO contributions to the 1-RDM
                   do i=1,ndocc
                      imo=docc(i)
                      dmat(imo,imo,ista)=dmat(imo,imo,ista)+2.0d0*c2
                   enddo
                   
                enddo
          
             enddo
             
          enddo
             
       enddo

    endif

!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
    if (cfg%n1I1E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h

          ! Loop over the 1I1E configurations generated by the
          ! current 2-hole configuration
          do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1

             ! Get the lists of singly-occupied and doubly-occupied MOs
             call sop_socc_list(cfg%sop1I1E(:,:,ioff),n_int,socc,&
                  nmo,nsocc)
             call sop_docc_list(cfg%sop1I1E(:,:,ioff),n_int,docc,&
                  nmo,ndocc)

             ! Loop over roots
             do ista=1,nroots
                
                ! Loop over CSFs generated by this configuration
                do icsf=cfg%csfs1I1E(ioff),cfg%csfs1I1E(ioff+1)-1
                   
                   ! Coefficient squared
                   c2=vec(icsf,ista)**2
                   
                   ! Singly-occupied MO contributions to the 1-RDM
                   do i=1,nsocc
                      imo=socc(i)
                      dmat(imo,imo,ista)=dmat(imo,imo,ista)+c2
                   enddo
                   
                   ! Doubly-occupied MO contributions to the 1-RDM
                   do i=1,ndocc
                      imo=docc(i)
                      dmat(imo,imo,ista)=dmat(imo,imo,ista)+2.0d0*c2
                   enddo
                   
                enddo
          
             enddo
             
          enddo
             
       enddo

    endif

!----------------------------------------------------------------------
! Sanity check on the traces of the 1-RDMs
!----------------------------------------------------------------------
    ! Loop over roots
    do ista=1,nroots
       trace=0.0d0
       ! Loop over MOs
       do imo=1,nmo
          trace=trace+dmat(imo,imo,ista)
       enddo
       ! Exit here if the trace of the 1-RDM does not equal the
       ! no. electrons
       if (abs(trace-nelB) > 1e-10_dp) then
          errmsg='Incorrect Tr(rho) in rdm_mrci_diag'
          call error_control
       endif
    enddo
    
    return
    
  end subroutine rdm_diag
    
!######################################################################
! rdm_mrci_0h_0h: Calculation of the Ref-Ref contributions to the MRCI
!                 1-RDMs
!######################################################################  
  subroutine rdm_0h_0h(cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! No. CSFs
    integer(is), intent(in) :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in) :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)

    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)
    
    ! Creation/annihilation operator indices
    integer(is), parameter  :: maxexci=1
    integer(is)             :: hlist(maxexci),plist(maxexci)

    ! Working arrays
    integer(ib)             :: kconf_full(n_int,2)
    integer(ib)             :: ksop_full(n_int,2)
        
    ! Everything else
    integer(is)             :: kconf,bconf,nexci,n_int_I
    integer(is)             :: knsp,bnsp,knopen,bnopen
    integer(is)             :: i,a,ista,ikcsf,ibcsf,komega,bomega
    real(dp)                :: kcoe,bcoe
    real(dp)                :: prod
    
!----------------------------------------------------------------------
! Contributions from bra and ket reference space CSFs
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I

    ! Loop over ket configurations
    do kconf=1,cfg%n0h-1

       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_I,:)=cfg%conf0h(:,:,kconf)
       ksop_full(1:n_int_I,:)=cfg%sop0h(:,:,kconf)

       ! Number of open shells in the ket configuration
       knopen=sop_nopen(cfg%sop0h(:,:,kconf),n_int_I)

       ! Number of ket CSFs
       knsp=ncsfs(knopen)
       
       ! Get the number of open shells preceding each ket conf MO
       call nobefore(ksop_full,nbefore)
       
       ! Loop over bra configurations
       do bconf=kconf+1,cfg%n0h

          ! Compute the excitation degree between the two
          ! configurations
          nexci=exc_degree_conf(cfg%conf0h(:,:,kconf),&
               cfg%conf0h(:,:,bconf),n_int_I)

          ! Cycle if the excitation degree is not equal to 1
          if (nexci /= 1) cycle

          ! Number of open shells in the bra configuration
          bnopen=sop_nopen(cfg%sop0h(:,:,bconf),n_int_I)

          ! Number of bra CSFs
          bnsp=ncsfs(bnopen)
          
          ! Get the indices of the MOs involved in the excitation
          hlist=0
          plist=0
          call get_exci_indices(cfg%conf0h(:,:,kconf),&
               cfg%conf0h(:,:,bconf),n_int_I,hlist(1:nexci),&
               plist(1:nexci),nexci)

          ! Get the spin-coupling coefficients
          spincp(1:knsp,1:bnsp)=spincp_coeff(knsp,bnsp,ksop_full,&
               plist(1),hlist(1),knopen,nbefore)

          ! Idices of the 1-RDM elements
          i=hlist(1)
          a=plist(1)

          ! Loop over roots
          do ista=1,nroots

             ! Loop over bra CSFs
             bomega=0
             do ibcsf=cfg%csfs0h(bconf),cfg%csfs0h(bconf+1)-1
                bomega=bomega+1
                bcoe=vec(ibcsf,ista)

                ! Loop over ket CSFs
                komega=0
                do ikcsf=cfg%csfs0h(kconf),cfg%csfs0h(kconf+1)-1
                   komega=komega+1
                   kcoe=vec(ikcsf,ista)

                   ! Contribution to the 1-RDM
                   prod=kcoe*bcoe*spincp(komega,bomega)
                   dmat(i,a,ista)=dmat(i,a,ista)+prod
                   dmat(a,i,ista)=dmat(a,i,ista)+prod
                   
                enddo
                   
             enddo
                
          enddo
          
       enddo
       
    enddo

    return
    
  end subroutine rdm_0h_0h

!######################################################################
! rdm_mrci_0h_1h: Calculation of the Ref-1I and Ref-1E contributions
!                 to the MRCI 1-RDMs
!######################################################################  
  subroutine rdm_0h_1h(cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! No. CSFs
    integer(is), intent(in) :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in) :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)

    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)

    ! Working arrays
    integer(ib)             :: kconf_full(n_int,2)
    integer(ib)             :: ksop_full(n_int,2)
    
    ! Everything else
    integer(is)             :: kconf,n,nac,nexci,n_int_I
    integer(is)             :: knsp,knopen

!----------------------------------------------------------------------
! Contributions ket reference and bra 1I & 1E CSFs
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
       
       ! Get the number of open shells preceding each ket conf MO
       call nobefore(ksop_full,nbefore)

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h

          ! Number of creation and annihilation operators linking the
          ! reference and 1-hole configurations
          nac=n_create_annihilate(cfg%conf0h(1:n_int_I,:,kconf), &
               cfg%conf1h(1:n_int_I,:,n),n_int_I)

          ! Ref - 1I contributions
          if (nac <= 3 .and. cfg%n1I >0 &
               .and. cfg%off1I(n) /= cfg%off1I(n+1)) then
             call rdm_0h_1I(n,kconf,ksop_full,nbefore,knopen,knsp,&
                  n_int_I,cfg,csfdim,nroots,vec,dmat)
          endif

          ! Ref - 1E contributions
          if (nac <= 1 .and. cfg%n1E >0 &
               .and. cfg%off1E(n) /= cfg%off1E(n+1)) then
             call rdm_0h_1E(n,kconf,ksop_full,nbefore,knopen,knsp,&
                  n_int_I,cfg,csfdim,nroots,vec,dmat)
          endif
             
       enddo
       
    enddo
    
    return
    
  end subroutine rdm_0h_1h

!######################################################################
! rdm_mrci_0h_2h: Calculation of the Ref-2I, Ref-2E and Ref-1I1E
!                 contributions to the MRCI 1-RDMs
!######################################################################  
  subroutine rdm_0h_2h(cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! No. CSFs
    integer(is), intent(in) :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in) :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)

    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)

    ! Working arrays
    integer(ib)             :: kconf_full(n_int,2)
    integer(ib)             :: ksop_full(n_int,2)
    
    ! Everything else
    integer(is)             :: kconf,n,nac,nexci,n_int_I
    integer(is)             :: knsp,knopen

!----------------------------------------------------------------------
! Contributions ket reference and bra 2I, 2E & 1I1E CSFs
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
       
       ! Get the number of open shells preceding each ket conf MO
       call nobefore(ksop_full,nbefore)

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h

          ! Number of creation and annihilation operators linking the
          ! reference and 2-hole configurations
          nac=n_create_annihilate(cfg%conf0h(1:n_int_I,:,kconf), &
               cfg%conf2h(1:n_int_I,:,n),n_int_I)

          ! Ref - 2I contributions
          if (nac <= 4 .and. cfg%n2I >0 &
               .and. cfg%off2I(n) /= cfg%off2I(n+1)) then
             call rdm_0h_2I(n,kconf,ksop_full,nbefore,knopen,knsp,&
                  n_int_I,cfg,csfdim,nroots,vec,dmat)
          endif

          ! Ref - 2E contributions
          if (nac == 0 .and. cfg%n2E >0 &
               .and. cfg%off2E(n) /= cfg%off2E(n+1)) then
             call rdm_0h_2E(n,kconf,ksop_full,nbefore,knopen,knsp,&
                  n_int_I,cfg,csfdim,nroots,vec,dmat)
          endif

          ! Ref - 1I1E contributions
          if (nac <= 2 .and. cfg%n1I1E >0 &
               .and. cfg%off1I1E(n) /= cfg%off1I1E(n+1)) then
             call rdm_0h_1I1E(n,kconf,ksop_full,nbefore,knopen,knsp,&
                  n_int_I,cfg,csfdim,nroots,vec,dmat)
          endif
          
       enddo
       
    enddo
    
    return
    
  end subroutine rdm_0h_2h

!######################################################################
! rdm_1h_1h: Calculation of the 1I-1I, 1I-1E and 1E-1E contributions
!            to the MRCI 1-RDMs
!######################################################################
  subroutine rdm_1h_1h(cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! No. CSFs
    integer(is), intent(in) :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in) :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)

    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)

    ! Working arrays
    integer(ib), allocatable   :: kconf_int(:,:)
    integer(ib), allocatable   :: ksop_int(:,:)
    integer(ib)                :: kconf_full(n_int,2)
    integer(ib)                :: bconf1h_full(n_int,2)
    integer(ib)                :: ksop_full(n_int,2)
    integer(is), parameter     :: maxexci=1
    integer(is)                :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)                :: n_int_I,knopen,knsp,kn,bn,nac,nac1
    integer(is)                :: ioff
    
!----------------------------------------------------------------------
! Contributions ket reference and bra 2I, 2E & 1I1E CSFs
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I

    ! Loop over ket 1-hole configurations
    do kn=1,cfg%n1h
    
       ! Loop over bra 1-hole configurations
       do bn=kn,cfg%n1h

          ! Number of creation and annihilation operators linking the
          ! bra and ket 1-hole configurations
          nac1=n_create_annihilate(cfg%conf1h(1:n_int_I,:,kn),&
               cfg%conf1h(1:n_int_I,:,bn),n_int_I)

          ! Cycle if the full configurations generated by the 1-hole
          ! configurations cannot be interacting wrt the singlet
          ! excitation operators E_a^i
          if (nac1 > 4) cycle

          !
          ! Ket 1I, bra 1I and 1E contributions
          !
          if (cfg%n1I > 0) then

             ! Loop over ket 1I configurations
             do ioff=cfg%off1I(kn),cfg%off1I(kn+1)-1

                ! Ket 1I configuration
                kconf_int=cfg%conf1I(1:n_int_I,:,ioff)
                ksop_int=cfg%sop1I(1:n_int_I,:,ioff)
                kconf_full=0_ib
                ksop_full=0_ib
                kconf_full(1:n_int_I,:)=kconf_int
                ksop_full(1:n_int_I,:)=ksop_int

                ! Number of creation and annihilation operators linking
                ! the ket 1I and bra 1-hole configurations
                nac=n_create_annihilate(kconf_int,&
                     cfg%conf1h(1:n_int_I,:,bn),n_int_I)

                ! Cycle if the the bra 1-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 1I configuration wrt the singlet excitation
                ! operators E_a^i
                if (nac1 > 3) cycle

                ! Number of open shells in the ket 1I configuration
                knopen=sop_nopen(ksop_int(1:n_int_I,:),n_int_I)
                
                ! Number of ket 1I CSFs
                knsp=ncsfs(knopen)

                ! Get the number of open shells preceding each ket
                ! conf MO
                call nobefore(ksop_full,nbefore)

                ! 1I - 1I contributions
                if (nac <= 3 &
                     .and. cfg%off1I(bn) /= cfg%off1I(bn+1)) then
                   call rdm_1I_1I(bn,ioff,kconf_int,ksop_int,n_int_I,&
                        nbefore,knopen,knsp,cfg,csfdim,nroots,vec,dmat)
                endif
                
                ! 1I - 1E contributions
                if (nac <= 1 &
                     .and. cfg%off1E(bn) /= cfg%off1E(bn+1)) then

                endif
                   
             enddo
                
          endif

          !
          ! Ket 1E, bra 1E matrix elements
          !
          ! Cycle if the bra 1-hole configuration doesn't generate
          ! any 1E configurations
          if (cfg%off1E(bn) == cfg%off1E(bn+1)) cycle
          
          ! Cycle if the the bra and ket 1-hole configurations
          ! cannot generate interacting 1E configurations wrt the
          ! singlet excitation operators E_a^i
          if (nac1 > 2) cycle
          
          ! Loop over ket 1E configurations
          do ioff=cfg%off1E(kn),cfg%off1E(kn+1)-1
             
             ! Ket 1E configuration
             kconf_full=cfg%conf1E(:,:,ioff)
             ksop_full=cfg%sop1E(:,:,ioff)
             
             ! Number of open shells in the ket 1E configuration
             knopen=sop_nopen(ksop_full,n_int)
                
             ! Number of ket 1E CSFs
             knsp=ncsfs(knopen)
             
             ! 1E - 1E matrix contributions
             
             
          enddo
             
       enddo

    enddo
          
    return
    
  end subroutine rdm_1h_1h
    
!######################################################################
! rdm_0h_1I: Calculation of the Ref-1I contributions to the 1-RDMs
!###################################################################### 
  subroutine rdm_0h_1I(n,kconf,ksop_full,nbefore,knopen,knsp,&
       n_int_I,cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! Index of the 1-hole configuration
    integer(is), intent(in) :: n

    ! Dimensions
    integer(is), intent(in) :: n_int_I,csfdim,nroots

    ! Ket reference space configuration
    integer(is), intent(in) :: kconf
    integer(ib), intent(in) :: ksop_full(n_int,2)
    integer(is), intent(in) :: knopen,knsp
    integer(is), intent(in) :: nbefore(nmo)

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)
    
    ! Working arrays
    integer(ib)             :: bconf_int(n_int_I,2)
    integer(ib)             :: bsop_int(n_int_I,2)
    integer(is), parameter  :: maxexci=1
    integer(is)             :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: bnsp,bnopen
    integer(is)             :: i,a,ista,ioff,ikcsf,ibcsf,komega,bomega
    integer(is)             :: nexci
    real(dp)                :: kcoe,bcoe
    real(dp)                :: prod

!----------------------------------------------------------------------
! Compute the R-1I contributions to the 1-RDMs
!----------------------------------------------------------------------
    ! Loop over 1I configurations generated by the current
    ! 1-hole configuration
    do ioff=cfg%off1I(n),cfg%off1I(n+1)-1

       ! Bra 1I configuration in the internal MO space
       bconf_int=0_ib
       bsop_int=0_ib
       bconf_int=cfg%conf1I(1:n_int_I,:,ioff)
       bsop_int=cfg%sop1I(1:n_int_I,:,ioff)
       
       ! Compute the excitation degree between the two
       ! configurations
       nexci=exc_degree_conf(cfg%conf0h(:,:,kconf),bconf_int,n_int_I)

       ! Cycle if the excitation degree is not equal to 1
       if (nexci /= 1) cycle

       ! Number of open shells in the bra configuration
       bnopen=sop_nopen(bsop_int,n_int_I)
    
       ! Number of bra CSFs
       bnsp=ncsfs(bnopen)
    
       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(cfg%conf0h(:,:,kconf),bconf_int,&
            n_int_I,hlist(1:nexci),plist(1:nexci),nexci)

       ! Get the spin-coupling coefficients
       spincp(1:knsp,1:bnsp)=spincp_coeff(knsp,bnsp,ksop_full,&
            plist(1),hlist(1),knopen,nbefore)

       ! Idices of the 1-RDM elements
       i=hlist(1)
       a=plist(1)
       
       ! Loop over roots
       do ista=1,nroots
          
          ! Loop over bra CSFs
          bomega=0
          do ibcsf=cfg%csfs1I(ioff),cfg%csfs1I(ioff+1)-1
             bomega=bomega+1
             bcoe=vec(ibcsf,ista)
             
             ! Loop over ket CSFs
             komega=0
             do ikcsf=cfg%csfs0h(kconf),cfg%csfs0h(kconf+1)-1
                komega=komega+1
                kcoe=vec(ikcsf,ista)
                
                ! Contribution to the 1-RDM
                prod=kcoe*bcoe*spincp(komega,bomega)
                dmat(i,a,ista)=dmat(i,a,ista)+prod
                dmat(a,i,ista)=dmat(a,i,ista)+prod
                
             enddo
             
          enddo
          
       enddo

    enddo
       
    return
    
  end subroutine rdm_0h_1I

!######################################################################
! rdm_0h_1E: Calculation of the Ref-1E contributions to the 1-RDMs
!###################################################################### 
  subroutine rdm_0h_1E(n,kconf,ksop_full,nbefore,knopen,knsp,&
       n_int_I,cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! Index of the 1-hole configuration
    integer(is), intent(in) :: n

    ! Dimensions
    integer(is), intent(in) :: n_int_I,csfdim,nroots

    ! Ket reference space configuration
    integer(is), intent(in) :: kconf
    integer(ib), intent(in) :: ksop_full(n_int,2)
    integer(is), intent(in) :: knopen,knsp
    integer(is), intent(in) :: nbefore(nmo)

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)
    
    ! Working arrays
    integer(ib)             :: bconf_int(n_int_I,2)
    integer(ib)             :: bsop_int(n_int_I,2)
    integer(is), parameter  :: maxexci=1
    integer(is)             :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: bnsp,bnopen
    integer(is)             :: i,a,ista,ioff,ikcsf,ibcsf,komega,bomega
    integer(is)             :: nexci
    real(dp)                :: kcoe,bcoe
    real(dp)                :: prod

!----------------------------------------------------------------------
! Compute the R-1E contributions to the 1-RDMs
!----------------------------------------------------------------------
    ! Loop over 1E configurations generated by the current
    ! 1-hole configuration
    do ioff=cfg%off1E(n),cfg%off1E(n+1)-1

       ! Bra 1E configuration in the internal MO space
       bconf_int=0_ib
       bsop_int=0_ib
       bconf_int=cfg%conf1E(1:n_int_I,:,ioff)
       bsop_int=cfg%sop1E(1:n_int_I,:,ioff)
       
       ! Compute the excitation degree between the two
       ! configurations
       nexci=exc_degree_conf(cfg%conf0h(:,:,kconf),bconf_int,n_int_I)

       ! Cycle if the excitation degree is not equal to 1
       if (nexci /= 1) cycle

       ! Number of open shells in the bra configuration
       bnopen=sop_nopen(bsop_int,n_int_I)
    
       ! Number of bra CSFs
       bnsp=ncsfs(bnopen)
    
       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(cfg%conf0h(:,:,kconf),bconf_int,&
            n_int_I,hlist(1:nexci),plist(1:nexci),nexci)

       ! Get the spin-coupling coefficients
       spincp(1:knsp,1:bnsp)=spincp_coeff(knsp,bnsp,ksop_full,&
            plist(1),hlist(1),knopen,nbefore)

       ! Idices of the 1-RDM elements
       i=hlist(1)
       a=plist(1)
       
       ! Loop over roots
       do ista=1,nroots
          
          ! Loop over bra CSFs
          bomega=0
          do ibcsf=cfg%csfs1E(ioff),cfg%csfs1E(ioff+1)-1
             bomega=bomega+1
             bcoe=vec(ibcsf,ista)
             
             ! Loop over ket CSFs
             komega=0
             do ikcsf=cfg%csfs0h(kconf),cfg%csfs0h(kconf+1)-1
                komega=komega+1
                kcoe=vec(ikcsf,ista)
                
                ! Contribution to the 1-RDM
                prod=kcoe*bcoe*spincp(komega,bomega)
                dmat(i,a,ista)=dmat(i,a,ista)+prod
                dmat(a,i,ista)=dmat(a,i,ista)+prod
                
             enddo
             
          enddo
          
       enddo
       
    enddo
       
    return
    
  end subroutine rdm_0h_1E

!######################################################################
! rdm_0h_2I: Calculation of the Ref-2I contributions to the 1-RDMs
!###################################################################### 
  subroutine rdm_0h_2I(n,kconf,ksop_full,nbefore,knopen,knsp,&
       n_int_I,cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! Index of the 2-hole configuration
    integer(is), intent(in) :: n

    ! Dimensions
    integer(is), intent(in) :: n_int_I,csfdim,nroots

    ! Ket reference space configuration
    integer(is), intent(in) :: kconf
    integer(ib), intent(in) :: ksop_full(n_int,2)
    integer(is), intent(in) :: knopen,knsp
    integer(is), intent(in) :: nbefore(nmo)

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)
    
    ! Working arrays
    integer(ib)             :: bconf_int(n_int_I,2)
    integer(ib)             :: bsop_int(n_int_I,2)
    integer(is), parameter  :: maxexci=1
    integer(is)             :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: bnsp,bnopen
    integer(is)             :: i,a,ista,ioff,ikcsf,ibcsf,komega,bomega
    integer(is)             :: nexci
    real(dp)                :: kcoe,bcoe
    real(dp)                :: prod

!----------------------------------------------------------------------
! Compute the R-2I contributions to the 1-RDMs
!----------------------------------------------------------------------
    ! Loop over 2I configurations generated by the current
    ! 2-hole configuration
    do ioff=cfg%off2I(n),cfg%off2I(n+1)-1

       ! Bra 2I configuration in the internal MO space
       bconf_int=0_ib
       bsop_int=0_ib
       bconf_int=cfg%conf2I(1:n_int_I,:,ioff)
       bsop_int=cfg%sop2I(1:n_int_I,:,ioff)
       
       ! Compute the excitation degree between the two
       ! configurations
       nexci=exc_degree_conf(cfg%conf0h(:,:,kconf),bconf_int,n_int_I)

       ! Cycle if the excitation degree is not equal to 1
       if (nexci /= 1) cycle

       ! Number of open shells in the bra configuration
       bnopen=sop_nopen(bsop_int,n_int_I)
    
       ! Number of bra CSFs
       bnsp=ncsfs(bnopen)
    
       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(cfg%conf0h(:,:,kconf),bconf_int,&
            n_int_I,hlist(1:nexci),plist(1:nexci),nexci)

       ! Get the spin-coupling coefficients
       spincp(1:knsp,1:bnsp)=spincp_coeff(knsp,bnsp,ksop_full,&
            plist(1),hlist(1),knopen,nbefore)

       ! Idices of the 1-RDM elements
       i=hlist(1)
       a=plist(1)
       
       ! Loop over roots
       do ista=1,nroots
          
          ! Loop over bra CSFs
          bomega=0
          do ibcsf=cfg%csfs2I(ioff),cfg%csfs2I(ioff+1)-1
             bomega=bomega+1
             bcoe=vec(ibcsf,ista)
             
             ! Loop over ket CSFs
             komega=0
             do ikcsf=cfg%csfs0h(kconf),cfg%csfs0h(kconf+1)-1
                komega=komega+1
                kcoe=vec(ikcsf,ista)
                
                ! Contribution to the 1-RDM
                prod=kcoe*bcoe*spincp(komega,bomega)
                dmat(i,a,ista)=dmat(i,a,ista)+prod
                dmat(a,i,ista)=dmat(a,i,ista)+prod
                
             enddo
             
          enddo
          
       enddo

    enddo
       
    return
    
  end subroutine rdm_0h_2I

!######################################################################
! rdm_0h_2E: Calculation of the Ref-2E contributions to the 1-RDMs
!###################################################################### 
  subroutine rdm_0h_2E(n,kconf,ksop_full,nbefore,knopen,knsp,&
       n_int_I,cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! Index of the 2-hole configuration
    integer(is), intent(in) :: n

    ! Dimensions
    integer(is), intent(in) :: n_int_I,csfdim,nroots

    ! Ket reference space configuration
    integer(is), intent(in) :: kconf
    integer(ib), intent(in) :: ksop_full(n_int,2)
    integer(is), intent(in) :: knopen,knsp
    integer(is), intent(in) :: nbefore(nmo)

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)
    
    ! Working arrays
    integer(ib)             :: bconf_int(n_int_I,2)
    integer(ib)             :: bsop_int(n_int_I,2)
    integer(is), parameter  :: maxexci=1
    integer(is)             :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: bnsp,bnopen
    integer(is)             :: i,a,ista,ioff,ikcsf,ibcsf,komega,bomega
    integer(is)             :: nexci
    real(dp)                :: kcoe,bcoe
    real(dp)                :: prod

!----------------------------------------------------------------------
! Compute the R-2E contributions to the 1-RDMs
!----------------------------------------------------------------------
    ! Loop over 2E configurations generated by the current
    ! 2-hole configuration
    do ioff=cfg%off2E(n),cfg%off2E(n+1)-1

       ! Bra 2E configuration in the internal MO space
       bconf_int=0_ib
       bsop_int=0_ib
       bconf_int=cfg%conf2E(1:n_int_I,:,ioff)
       bsop_int=cfg%sop2E(1:n_int_I,:,ioff)
       
       ! Compute the excitation degree between the two
       ! configurations
       nexci=exc_degree_conf(cfg%conf0h(:,:,kconf),bconf_int,n_int_I)

       ! Cycle if the excitation degree is not equal to 1
       if (nexci /= 1) cycle

       ! Number of open shells in the bra configuration
       bnopen=sop_nopen(bsop_int,n_int_I)
    
       ! Number of bra CSFs
       bnsp=ncsfs(bnopen)
    
       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(cfg%conf0h(:,:,kconf),bconf_int,&
            n_int_I,hlist(1:nexci),plist(1:nexci),nexci)

       ! Get the spin-coupling coefficients
       spincp(1:knsp,1:bnsp)=spincp_coeff(knsp,bnsp,ksop_full,&
            plist(1),hlist(1),knopen,nbefore)

       ! Idices of the 1-RDM elements
       i=hlist(1)
       a=plist(1)
       
       ! Loop over roots
       do ista=1,nroots
          
          ! Loop over bra CSFs
          bomega=0
          do ibcsf=cfg%csfs2E(ioff),cfg%csfs2E(ioff+1)-1
             bomega=bomega+1
             bcoe=vec(ibcsf,ista)
             
             ! Loop over ket CSFs
             komega=0
             do ikcsf=cfg%csfs0h(kconf),cfg%csfs0h(kconf+1)-1
                komega=komega+1
                kcoe=vec(ikcsf,ista)
                
                ! Contribution to the 1-RDM
                prod=kcoe*bcoe*spincp(komega,bomega)
                dmat(i,a,ista)=dmat(i,a,ista)+prod
                dmat(a,i,ista)=dmat(a,i,ista)+prod
                
             enddo
             
          enddo
          
       enddo

    enddo
       
    return
    
  end subroutine rdm_0h_2E

!######################################################################
! rdm_0h_1I1E: Calculation of the Ref-1I1E contributions to the 1-RDMs
!###################################################################### 
  subroutine rdm_0h_1I1E(n,kconf,ksop_full,nbefore,knopen,knsp,&
       n_int_I,cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! Index of the 2-hole configuration
    integer(is), intent(in) :: n

    ! Dimensions
    integer(is), intent(in) :: n_int_I,csfdim,nroots

    ! Ket reference space configuration
    integer(is), intent(in) :: kconf
    integer(ib), intent(in) :: ksop_full(n_int,2)
    integer(is), intent(in) :: knopen,knsp
    integer(is), intent(in) :: nbefore(nmo)

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)
    
    ! Working arrays
    integer(ib)             :: bconf_int(n_int_I,2)
    integer(ib)             :: bsop_int(n_int_I,2)
    integer(is), parameter  :: maxexci=1
    integer(is)             :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: bnsp,bnopen
    integer(is)             :: i,a,ista,ioff,ikcsf,ibcsf,komega,bomega
    integer(is)             :: nexci
    real(dp)                :: kcoe,bcoe
    real(dp)                :: prod

!----------------------------------------------------------------------
! Compute the R-1I1E contributions to the 1-RDMs
!----------------------------------------------------------------------
    ! Loop over 1I1E configurations generated by the current
    ! 2-hole configuration
    do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1

       ! Bra 1I1E configuration in the internal MO space
       bconf_int=0_ib
       bsop_int=0_ib
       bconf_int=cfg%conf1I1E(1:n_int_I,:,ioff)
       bsop_int=cfg%sop1I1E(1:n_int_I,:,ioff)
       
       ! Compute the excitation degree between the two
       ! configurations
       nexci=exc_degree_conf(cfg%conf0h(:,:,kconf),bconf_int,n_int_I)

       ! Cycle if the excitation degree is not equal to 1
       if (nexci /= 1) cycle

       ! Number of open shells in the bra configuration
       bnopen=sop_nopen(bsop_int,n_int_I)
    
       ! Number of bra CSFs
       bnsp=ncsfs(bnopen)
    
       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(cfg%conf0h(:,:,kconf),bconf_int,&
            n_int_I,hlist(1:nexci),plist(1:nexci),nexci)

       ! Get the spin-coupling coefficients
       spincp(1:knsp,1:bnsp)=spincp_coeff(knsp,bnsp,ksop_full,&
            plist(1),hlist(1),knopen,nbefore)

       ! Idices of the 1-RDM elements
       i=hlist(1)
       a=plist(1)
       
       ! Loop over roots
       do ista=1,nroots
          
          ! Loop over bra CSFs
          bomega=0
          do ibcsf=cfg%csfs1I1E(ioff),cfg%csfs1I1E(ioff+1)-1
             bomega=bomega+1
             bcoe=vec(ibcsf,ista)
             
             ! Loop over ket CSFs
             komega=0
             do ikcsf=cfg%csfs0h(kconf),cfg%csfs0h(kconf+1)-1
                komega=komega+1
                kcoe=vec(ikcsf,ista)
                
                ! Contribution to the 1-RDM
                prod=kcoe*bcoe*spincp(komega,bomega)
                dmat(i,a,ista)=dmat(i,a,ista)+prod
                dmat(a,i,ista)=dmat(a,i,ista)+prod
                
             enddo
             
          enddo
          
       enddo

    enddo
       
    return
    
  end subroutine rdm_0h_1I1E

!######################################################################
! rdm_1I_1I: Calculation of the 1I-1I contributions to the 1-RDMs
!###################################################################### 
  subroutine rdm_1I_1I(bn,ikconf1I,kconf_int,ksop_int,n_int_I,&
       nbefore,knopen,knsp,cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! Index of the bra 2-hole configuration
    integer(is), intent(in) :: bn

    ! Dimensions
    integer(is), intent(in) :: n_int_I,csfdim,nroots

    ! Ket reference space configuration
    integer(is), intent(in) :: ikconf1I
    integer(ib), intent(in) :: kconf_int(n_int_I,2)
    integer(ib), intent(in) :: ksop_int(n_int_I,2)
    integer(is), intent(in) :: knopen,knsp
    integer(is), intent(in) :: nbefore(nmo)

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)
    
    ! Working arrays
    integer(ib)             :: bconf_int(n_int_I,2)
    integer(ib)             :: bsop_int(n_int_I,2)
    integer(is), parameter  :: maxexci=1
    integer(is)             :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: bnsp,bnopen
    integer(is)             :: i,a,ista,ioff,ikcsf,ibcsf,komega,bomega
    integer(is)             :: nexci
    real(dp)                :: kcoe,bcoe
    real(dp)                :: prod

    return
    
  end subroutine rdm_1I_1I
    
!######################################################################
! spincp_coeff: Given a SOP and pair of creation/annihilation operator
!               indices, returns the complete set of spin-coupling
!               coefficients
!######################################################################
  function spincp_coeff(knsp,bnsp,sop,ac,ia,nopen,nbefore)

    use constants
    use bitglobal
    use pattern_indices
    use bitstrings
    use iomod
    
    implicit none

    ! Numbers of ket and bra CSFs
    integer(is), intent(in) :: knsp,bnsp

    ! SOP
    integer(ib), intent(in) :: sop(n_int,2)
    
    ! Creation/annihilation operator indices
    integer(is), intent(in) :: ac,ia

    ! Number of open shells
    integer(is), intent(in) :: nopen
    
    ! Number of open shells preceding each MO
    integer(is), intent(in) :: nbefore(nmo)

    ! Spin-coupling sub-case bit string encodings
    integer(is)             :: pattern
    integer(ib)             :: icase

    ! Function result
    real(dp)                :: spincp_coeff(knsp,bnsp)
    
    ! Everything else
    integer(is)             :: nc,na

!----------------------------------------------------------------------
! Get the pattern index and sub-case bit string for the spin-coupling
! coefficients
!----------------------------------------------------------------------
    ! No. open shells before the created electron
    nc=nbefore(ac)
    
    ! No. open shells before the annihilated electron
    na=nbefore(ia)

    ! Spin-coupling sub-case bit string
    icase=get_icase(sop,ac,ia)

    ! Pattern index
    pattern=pattern_index(sop,ac,ia,nc,na,nopen,icase)
    
!----------------------------------------------------------------------
! Fill in the array of spin-coupling coefficients
!----------------------------------------------------------------------
    select case(icase)
    case(i1a)
       spincp_coeff(1:knsp,1:bnsp)=spincp1(1:knsp,1:bnsp,pattern)
    case(i1b)
       spincp_coeff(1:knsp,1:bnsp)=-spincp1(1:knsp,1:bnsp,pattern)
    case(i2a)
       spincp_coeff(1:knsp,1:bnsp)=spincp2(1:knsp,1:bnsp,pattern)
    case(i2b)
       spincp_coeff(1:knsp,1:bnsp)=transpose(spincp2(1:bnsp,1:knsp,pattern))
    case default
       errmsg='Unrecognised icase value in spincp_coeff'
       call error_control
    end select
    
    return
    
  end function spincp_coeff
  
!######################################################################
  
end module rdm
