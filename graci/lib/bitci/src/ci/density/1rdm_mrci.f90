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
!            a specified set of states of a single irrep.
!######################################################################
!            Note that the 1-RDMS have the 'Canonical' MO indexing
!######################################################################
  subroutine rdm_mrci(cfg,csfdim,nroots,vec,rho)

    use constants
    use bitglobal
    use conftype
    use iomod
    use timing
    use utils
    
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
    real(dp), intent(out)   :: rho(nmo,nmo,nroots)
   
    ! Timing variables
    real(dp)                :: tcpu_start,tcpu_end,twall_start,&
                               twall_end
    
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)
    
!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    rho=0.0d0

    allocate(spincp(ncsfs(nomax),ncsfs(nomax)))
    spincp=0.0d0

!----------------------------------------------------------------------
! (1) On-diagonal elements
!----------------------------------------------------------------------
    call rdm_diag(cfg,csfdim,nroots,vec,rho)

!----------------------------------------------------------------------
! (2) Ref - Ref contributions to the 1-RDMs
!----------------------------------------------------------------------
    call rdm_0h_0h(cfg,csfdim,nroots,vec,rho)

!----------------------------------------------------------------------
! (3) Ref - 1H contributions to the 1-RDMs
!----------------------------------------------------------------------
    call rdm_0h_1h(cfg,csfdim,nroots,vec,rho)

!----------------------------------------------------------------------
! (4) Ref - 2H contributions to the 1-RDMs
!----------------------------------------------------------------------
    call rdm_0h_2h(cfg,csfdim,nroots,vec,rho)

!----------------------------------------------------------------------
! (5) 1H - 1H contributions to the 1-RDMs
!----------------------------------------------------------------------
    call rdm_1h_1h(cfg,csfdim,nroots,vec,rho)

!----------------------------------------------------------------------
! (6) 2H - 1H contributions to the 1-RDMs
!----------------------------------------------------------------------
    call rdm_2h_1h(cfg,csfdim,nroots,vec,rho)
    
!----------------------------------------------------------------------
! (6) 2H - 2H contributions to the 1-RDMs
!----------------------------------------------------------------------
    call rdm_2h_2h(cfg,csfdim,nroots,vec,rho)

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
  subroutine rdm_diag(cfg,csfdim,nroots,vec,rho)

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
    real(dp), intent(inout) :: rho(nmo,nmo,nroots)

    ! Working arrays
    integer(ib)             :: conf_full(n_int,2)
    integer(ib)             :: sop_full(n_int,2)
        
    ! MO classes
    integer(is)             :: socc(nmo),docc(nmo),unocc(nmo)
    integer(is)             :: nsocc,ndocc,nunocc
    
    ! Everything else
    integer(is)             :: iconf,icsf
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
    ! Loop over reference configurations
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
                imo=cfg%m2c(socc(i))
                rho(imo,imo,ista)=rho(imo,imo,ista)+c2
             enddo

             ! Doubly-occupied MO contributions to the 1-RDM
             do i=1,ndocc
                imo=cfg%m2c(docc(i))
                rho(imo,imo,ista)=rho(imo,imo,ista)+2.0d0*c2
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
                      imo=cfg%m2c(socc(i))
                      rho(imo,imo,ista)=rho(imo,imo,ista)+c2
                   enddo
                   
                   ! Doubly-occupied MO contributions to the 1-RDM
                   do i=1,ndocc
                      imo=cfg%m2c(docc(i))
                      rho(imo,imo,ista)=rho(imo,imo,ista)+2.0d0*c2
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
                      imo=cfg%m2c(socc(i))
                      rho(imo,imo,ista)=rho(imo,imo,ista)+c2
                   enddo
                   
                   ! Doubly-occupied MO contributions to the 1-RDM
                   do i=1,ndocc
                      imo=cfg%m2c(docc(i))
                      rho(imo,imo,ista)=rho(imo,imo,ista)+2.0d0*c2
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
                      imo=cfg%m2c(socc(i))
                      rho(imo,imo,ista)=rho(imo,imo,ista)+c2
                   enddo
                   
                   ! Doubly-occupied MO contributions to the 1-RDM
                   do i=1,ndocc
                      imo=cfg%m2c(docc(i))
                      rho(imo,imo,ista)=rho(imo,imo,ista)+2.0d0*c2
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
                      imo=cfg%m2c(socc(i))
                      rho(imo,imo,ista)=rho(imo,imo,ista)+c2
                   enddo
                   
                   ! Doubly-occupied MO contributions to the 1-RDM
                   do i=1,ndocc
                      imo=cfg%m2c(docc(i))
                      rho(imo,imo,ista)=rho(imo,imo,ista)+2.0d0*c2
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
                      imo=cfg%m2c(socc(i))
                      rho(imo,imo,ista)=rho(imo,imo,ista)+c2
                   enddo
                   
                   ! Doubly-occupied MO contributions to the 1-RDM
                   do i=1,ndocc
                      imo=cfg%m2c(docc(i))
                      rho(imo,imo,ista)=rho(imo,imo,ista)+2.0d0*c2
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
          trace=trace+rho(imo,imo,ista)
       enddo
       ! Exit here if the trace of the 1-RDM does not equal the
       ! no. electrons
       if (abs(trace-nel) > 1e-10_dp) then
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
  subroutine rdm_0h_0h(cfg,csfdim,nroots,vec,rho)

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
    real(dp), intent(inout) :: rho(nmo,nmo,nroots)

    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)
    
    ! Creation/annihilation operator indices
    integer(is), parameter  :: maxexci=1
    integer(is)             :: hlist(maxexci),plist(maxexci)

    ! Working arrays
    integer(ib)             :: kconf_full(n_int,2)
    integer(ib)             :: ksop_full(n_int,2)
        
    ! Everything else
    integer(is)             :: ikconf,ibconf,nexci,n_int_I
    integer(is)             :: knsp,bnsp,knopen,bnopen
    integer(is)             :: i,a,ista,ikcsf,ibcsf,komega,bomega
    real(dp)                :: kcoe,bcoe
    real(dp)                :: prod
    
!----------------------------------------------------------------------
! Contributions from bra and ket reference space CSFs
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I

    ! Loop over ket configurations
    do ikconf=1,cfg%n0h-1

       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_I,:)=cfg%conf0h(:,:,ikconf)
       ksop_full(1:n_int_I,:)=cfg%sop0h(:,:,ikconf)

       ! Number of open shells in the ket configuration
       knopen=sop_nopen(cfg%sop0h(:,:,ikconf),n_int_I)

       ! Number of ket CSFs
       knsp=ncsfs(knopen)
       
       ! Get the number of open shells preceding each ket conf MO
       call nobefore(ksop_full,nbefore)
       
       ! Loop over bra configurations
       do ibconf=ikconf+1,cfg%n0h

          ! Compute the excitation degree between the two
          ! configurations
          nexci=exc_degree_conf(cfg%conf0h(:,:,ikconf),&
               cfg%conf0h(:,:,ibconf),n_int_I)

          ! Cycle if the excitation degree is not equal to 1
          if (nexci /= 1) cycle

          ! Number of open shells in the bra configuration
          bnopen=sop_nopen(cfg%sop0h(:,:,ibconf),n_int_I)

          ! Number of bra CSFs
          bnsp=ncsfs(bnopen)
          
          ! Get the indices of the MOs involved in the excitation
          hlist=0
          plist=0
          call get_exci_indices(cfg%conf0h(:,:,ikconf),&
               cfg%conf0h(:,:,ibconf),n_int_I,hlist(1:nexci),&
               plist(1:nexci),nexci)

          ! Get the spin-coupling coefficients
          spincp(1:knsp,1:bnsp)=spincp_coeff(knsp,bnsp,ksop_full,&
               plist(1),hlist(1),knopen,nbefore)

          ! Idices of the 1-RDM elements
          i=cfg%m2c(hlist(1))
          a=cfg%m2c(plist(1))

          ! Loop over roots
          do ista=1,nroots

             ! Loop over bra CSFs
             bomega=0
             do ibcsf=cfg%csfs0h(ibconf),cfg%csfs0h(ibconf+1)-1
                bomega=bomega+1
                bcoe=vec(ibcsf,ista)

                ! Loop over ket CSFs
                komega=0
                do ikcsf=cfg%csfs0h(ikconf),cfg%csfs0h(ikconf+1)-1
                   komega=komega+1
                   kcoe=vec(ikcsf,ista)

                   ! Contribution to the 1-RDM
                   prod=kcoe*bcoe*spincp(komega,bomega)
                   rho(i,a,ista)=rho(i,a,ista)+prod
                   rho(a,i,ista)=rho(a,i,ista)+prod
                   
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
  subroutine rdm_0h_1h(cfg,csfdim,nroots,vec,rho)

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
    real(dp), intent(inout) :: rho(nmo,nmo,nroots)

    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)

    ! Working arrays
    integer(ib)             :: kconf_full(n_int,2)
    integer(ib)             :: ksop_full(n_int,2)
    
    ! Everything else
    integer(is)             :: ikconf,n,nac,nexci,n_int_I
    integer(is)             :: knsp,knopen

!----------------------------------------------------------------------
! Contributions ket reference and bra 1I & 1E CSFs
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I
    
    ! Loop over ket reference configurations
    do ikconf=1,cfg%n0h

       ! Number of open shells in the ket configuration
       knopen=sop_nopen(cfg%sop0h(:,:,ikconf),n_int_I)
       
       ! Number of ket CSFs
       knsp=ncsfs(knopen)

       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_I,:)=cfg%conf0h(:,:,ikconf)
       ksop_full(1:n_int_I,:)=cfg%sop0h(:,:,ikconf)
       
       ! Get the number of open shells preceding each ket conf MO
       call nobefore(ksop_full,nbefore)

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h

          ! Number of creation and annihilation operators linking the
          ! reference and 1-hole configurations
          nac=n_create_annihilate(cfg%conf0h(1:n_int_I,:,ikconf), &
               cfg%conf1h(1:n_int_I,:,n),n_int_I)

          ! Ref - 1I contributions
          if (nac <= 3 .and. cfg%n1I > 0 &
               .and. cfg%off1I(n) /= cfg%off1I(n+1)) then
             call rdm_batch(&
                  n,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                  cfg%n1I,cfg%n0h,&       ! no. bra and ket confs
                  cfg%conf1I,cfg%sop1I,&  ! bra confs and SOPs
                  cfg%n1h,cfg%off1I,&     ! no. bra hole confs and offsets
                  cfg%csfs1I,cfg%csfs0h,& ! bra and ket CSF offsets
                  csfdim,nroots,vec,rho,cfg%m2c,.false.)
          endif
          
          ! Ref - 1E contributions
          if (nac <= 1 .and. cfg%n1E > 0 &
               .and. cfg%off1E(n) /= cfg%off1E(n+1)) then
             call rdm_batch(&
                  n,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                  cfg%n1E,cfg%n0h,&       ! no. bra and ket confs
                  cfg%conf1E,cfg%sop1E,&  ! bra confs and SOPs
                  cfg%n1h,cfg%off1E,&     ! no. bra hole confs and offsets
                  cfg%csfs1E,cfg%csfs0h,& ! bra and ket CSF offsets
                  csfdim,nroots,vec,rho,cfg%m2c,.false.)
          endif
             
       enddo
       
    enddo
    
    return
    
  end subroutine rdm_0h_1h

!######################################################################
! rdm_mrci_0h_2h: Calculation of the Ref-2I, Ref-2E and Ref-1I1E
!                 contributions to the MRCI 1-RDMs
!######################################################################  
  subroutine rdm_0h_2h(cfg,csfdim,nroots,vec,rho)

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
    real(dp), intent(inout) :: rho(nmo,nmo,nroots)

    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)

    ! Working arrays
    integer(ib)             :: kconf_full(n_int,2)
    integer(ib)             :: ksop_full(n_int,2)
    
    ! Everything else
    integer(is)             :: ikconf,n,nac,nexci,n_int_I
    integer(is)             :: knsp,knopen

!----------------------------------------------------------------------
! Contributions ket reference and bra 2I, 2E & 1I1E CSFs
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I
    
    ! Loop over ket reference configurations
    do ikconf=1,cfg%n0h

       ! Number of open shells in the ket configuration
       knopen=sop_nopen(cfg%sop0h(:,:,ikconf),n_int_I)
       
       ! Number of ket CSFs
       knsp=ncsfs(knopen)

       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_I,:)=cfg%conf0h(:,:,ikconf)
       ksop_full(1:n_int_I,:)=cfg%sop0h(:,:,ikconf)
       
       ! Get the number of open shells preceding each ket conf MO
       call nobefore(ksop_full,nbefore)

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h

          ! Number of creation and annihilation operators linking the
          ! reference and 2-hole configurations
          nac=n_create_annihilate(cfg%conf0h(1:n_int_I,:,ikconf), &
               cfg%conf2h(1:n_int_I,:,n),n_int_I)

          ! Ref - 2I contributions
          if (nac <= 4 .and. cfg%n2I > 0 &
               .and. cfg%off2I(n) /= cfg%off2I(n+1)) then
             call rdm_batch(&
                  n,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                  cfg%n2I,cfg%n0h,&       ! no. bra and ket confs
                  cfg%conf2I,cfg%sop2I,&  ! bra confs and SOPs
                  cfg%n2h,cfg%off2I,&     ! no. bra hole confs and offsets
                  cfg%csfs2I,cfg%csfs0h,& ! bra and ket CSF offsets
                  csfdim,nroots,vec,rho,cfg%m2c,.false.)
          endif

          ! Ref - 2E contributions
          if (nac == 0 .and. cfg%n2E > 0 &
               .and. cfg%off2E(n) /= cfg%off2E(n+1)) then
             call rdm_batch(&
                  n,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                  cfg%n2E,cfg%n0h,&       ! no. bra and ket confs
                  cfg%conf2E,cfg%sop2E,&  ! bra confs and SOPs
                  cfg%n2h,cfg%off2E,&     ! no. bra hole confs and offsets
                  cfg%csfs2E,cfg%csfs0h,& ! bra and ket CSF offsets
                  csfdim,nroots,vec,rho,cfg%m2c,.false.)
          endif

          ! Ref - 1I1E contributions
          if (nac <= 2 .and. cfg%n1I1E > 0 &
               .and. cfg%off1I1E(n) /= cfg%off1I1E(n+1)) then
             call rdm_batch(&
                  n,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                  cfg%n1I1E,cfg%n0h,&        ! no. bra and ket confs
                  cfg%conf1I1E,cfg%sop1I1E,& ! bra confs and SOPs
                  cfg%n2h,cfg%off1I1E,&      ! no. bra hole confs and offsets
                  cfg%csfs1I1E,cfg%csfs0h,&  ! bra and ket CSF offsets
                  csfdim,nroots,vec,rho,cfg%m2c,.false.)
          endif
          
       enddo
       
    enddo
    
    return
    
  end subroutine rdm_0h_2h

!######################################################################
! rdm_1h_1h: Calculation of the 1I-1I, 1I-1E and 1E-1E contributions
!            to the MRCI 1-RDMs
!######################################################################
  subroutine rdm_1h_1h(cfg,csfdim,nroots,vec,rho)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! No. CSFs
    integer(is), intent(in)  :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in)  :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! Eigenvectors
    real(dp), intent(in)     :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout)  :: rho(nmo,nmo,nroots)

    ! Number of open shells preceding each MO
    integer(is)              :: nbefore(nmo)

    ! Working arrays
    integer(ib), allocatable :: kconf_int(:,:)
    integer(ib), allocatable :: ksop_int(:,:)
    integer(ib)              :: kconf_full(n_int,2)
    integer(ib)              :: ksop_full(n_int,2)
    integer(is), parameter   :: maxexci=1
    integer(is)              :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)              :: n_int_I,knopen,knsp,kn,bn,nac,nac1
    integer(is)              :: ikconf
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I

    allocate(kconf_int(n_int_I,2))
    allocate(ksop_int(n_int_I,2))
    kconf_int=0_ib
    ksop_int=0_ib
    
!----------------------------------------------------------------------
! Compute the 1H-1H contributions to the 1-RDMs
!----------------------------------------------------------------------
    ! Loop over ket 1-hole configurations
    do kn=1,cfg%n1h
    
       ! Loop over bra 1-hole configurations
       do bn=1,cfg%n1h
       
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
             do ikconf=cfg%off1I(kn),cfg%off1I(kn+1)-1

                ! Ket 1I configuration
                kconf_int=cfg%conf1I(1:n_int_I,:,ikconf)
                ksop_int=cfg%sop1I(1:n_int_I,:,ikconf)
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
                if (nac <= 3 .and. bn >= kn &
                     .and. cfg%off1I(bn) /= cfg%off1I(bn+1)) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n1I,cfg%n1I,&       ! no. bra and ket confs
                        cfg%conf1I,cfg%sop1I,&  ! bra confs and SOPs
                        cfg%n1h,cfg%off1I,&     ! no. bra hole confs and offsets
                        cfg%csfs1I,cfg%csfs1I,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.true.)
                endif
                
                ! 1I - 1E contributions
                if (nac <= 1 .and. cfg%n1E > 0 &
                     .and. cfg%off1E(bn) /= cfg%off1E(bn+1)) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n1E,cfg%n1I,&       ! no. bra and ket confs
                        cfg%conf1E,cfg%sop1E,&  ! bra confs and SOPs
                        cfg%n1h,cfg%off1E,&     ! no. bra hole confs and offsets
                        cfg%csfs1E,cfg%csfs1I,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.false.)
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

          ! Skip duplicate 1E - 1E  elements
          if (bn < kn) cycle

          ! Cycle if there are no 1E configurations
          if (cfg%n1E == 0) cycle
          
          ! Loop over ket 1E configurations
          do ikconf=cfg%off1E(kn),cfg%off1E(kn+1)-1
             
             ! Ket 1E configuration
             kconf_full=cfg%conf1E(:,:,ikconf)
             ksop_full=cfg%sop1E(:,:,ikconf)
             
             ! Number of open shells in the ket 1E configuration
             knopen=sop_nopen(ksop_full,n_int)
                
             ! Number of ket 1E CSFs
             knsp=ncsfs(knopen)

             ! Get the number of open shells preceding each ket
             ! conf MO
             call nobefore(ksop_full,nbefore)
             
             ! 1E - 1E matrix contributions
             call rdm_batch(bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                  cfg%n1E,cfg%n1E,&       ! no. bra and ket confs
                  cfg%conf1E,cfg%sop1E,&  ! bra confs and SOPs
                  cfg%n1h,cfg%off1E,&     ! no. bra hole confs and offsets
                  cfg%csfs1E,cfg%csfs1E,& ! bra and ket CSF offsets
                  csfdim,nroots,vec,rho,cfg%m2c,.true.)
             
          enddo
             
       enddo

    enddo
          
    return
    
  end subroutine rdm_1h_1h

!######################################################################
! rdm_2h_1h: Calculation of the 2I-1I, 2I-1E, 2E-1I, 2E-1E, 1I1E-1I,
!            and 1I1E-1E contributions to the MRCI 1-RDMs
!######################################################################
  subroutine rdm_2h_1h(cfg,csfdim,nroots,vec,rho)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! No. CSFs
    integer(is), intent(in)  :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in)  :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! Eigenvectors
    real(dp), intent(in)     :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout)  :: rho(nmo,nmo,nroots)

    ! Number of open shells preceding each MO
    integer(is)              :: nbefore(nmo)

    ! Working arrays
    integer(ib), allocatable :: kconf_int(:,:)
    integer(ib), allocatable :: ksop_int(:,:)
    integer(ib)              :: kconf_full(n_int,2)
    integer(ib)              :: ksop_full(n_int,2)
    integer(ib)              :: bconf1h_full(n_int,2)
    integer(is), parameter   :: maxexci=1
    integer(is)              :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)              :: n_int_I,knopen,knsp,kn,bn,nac,nac1
    integer(is)              :: ikconf
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I

    allocate(kconf_int(n_int_I,2))
    allocate(ksop_int(n_int_I,2))
    kconf_int=0_ib
    ksop_int=0_ib

!----------------------------------------------------------------------
! Calculate the 2H-1H contributions to the 1-RDMs
!----------------------------------------------------------------------
    ! Loop over ket 2-hole configurations
    do kn=1,cfg%n2h
       
       ! Loop over bra 1-hole configurations
       do bn=1,cfg%n1h
          
          ! Bra 1-hole configuration in the full MO space
          bconf1h_full=0_ib
          bconf1h_full(1:n_int_I,:)=cfg%conf1h(1:n_int_I,:,bn)

          ! Number of creation and annihilation operators linking the
          ! bra and ket 1- and 2-hole configurations
          nac1=n_create_annihilate(cfg%conf2h(1:n_int_I,:,kn),&
               cfg%conf1h(1:n_int_I,:,bn),n_int_I)

          ! Cycle if the full configurations generated by the hole
          ! configurations cannot be interacting
          if (nac1 > 5) cycle

          !
          ! Ket: 2I
          ! Bra: 1I and 1E
          !
          if (cfg%n2I > 0) then

             ! Loop over ket 2I configurations
             do ikconf=cfg%off2I(kn),cfg%off2I(kn+1)-1

                ! Ket 2I configuration
                kconf_int=cfg%conf2I(1:n_int_I,:,ikconf)
                ksop_int=cfg%sop2I(1:n_int_I,:,ikconf)
                kconf_full=0_ib
                ksop_full=0_ib
                kconf_full(1:n_int_I,:)=kconf_int
                ksop_full(1:n_int_I,:)=ksop_int
                
                ! Number of creation and annihilation operators linking
                ! the ket 2I and bra 1-hole configurations
                nac=n_create_annihilate(kconf_int,&
                     cfg%conf1h(1:n_int_I,:,bn),n_int_I)
                
                ! Cycle if the the bra 1-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 2I configuration
                if (nac > 3) cycle

                ! Number of open shells in the ket 2I configuration
                knopen=sop_nopen(ksop_int(1:n_int_I,:),n_int_I)
                
                ! Number of ket 2I CSFs
                knsp=ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! 2I - 1I matrix contributions
                if (nac <= 3 &
                     .and. cfg%off1I(bn) /= cfg%off1I(bn+1)) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n1I,cfg%n2I,&       ! no. bra and ket confs
                        cfg%conf1I,cfg%sop1I,&  ! bra confs and SOPs
                        cfg%n1h,cfg%off1I,&     ! no. bra hole confs and offsets
                        cfg%csfs1I,cfg%csfs2I,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.false.)
                endif
                
                ! 2I - 1E matrix contributions
                if (nac <= 1 &
                     .and. cfg%off1E(bn) /= cfg%off1E(bn+1)) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n1E,cfg%n2I,&       ! no. bra and ket confs
                        cfg%conf1E,cfg%sop1E,&  ! bra confs and SOPs
                        cfg%n1h,cfg%off1E,&     ! no. bra hole confs and offsets
                        cfg%csfs1E,cfg%csfs2I,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.false.)
                endif
                
             enddo
             
          endif

          !
          ! Ket: 2E
          ! Bra: 1I and 1E
          !
          if (cfg%n2E > 0 .and. nac1 <= 1) then
             
             ! Loop over ket 2E configurations
             do ikconf=cfg%off2E(kn),cfg%off2E(kn+1)-1

                ! Ket 2E configuration
                kconf_full=cfg%conf2E(:,:,ikconf)
                ksop_full=cfg%sop2E(:,:,ikconf)

                ! Number of open shells in the ket 2E configuration
                knopen=sop_nopen(ksop_full,n_int)
                
                ! Number of ket 2E CSFs
                knsp=ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! 2E - 1I matrix contributions
                if (cfg%n1I > 0 &
                     .and. cfg%off1I(bn) /= cfg%off1I(bn+1)) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n1I,cfg%n2E,&       ! no. bra and ket confs
                        cfg%conf1I,cfg%sop1I,&  ! bra confs and SOPs
                        cfg%n1h,cfg%off1I,&     ! no. bra hole confs and offsets
                        cfg%csfs1I,cfg%csfs2E,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.false.)
                endif
                   
                ! 2E - 1E matrix contributions
                if (cfg%off1E(bn) /= cfg%off1E(bn+1)) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n1E,cfg%n2E,&       ! no. bra and ket confs
                        cfg%conf1E,cfg%sop1E,&  ! bra confs and SOPs
                        cfg%n1h,cfg%off1E,&     ! no. bra hole confs and offsets
                        cfg%csfs1E,cfg%csfs2E,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.false.)
                endif
                
             enddo

          endif

          !
          ! Ket: 1I1E
          ! Bra: 1I and 1E
          !
          if (cfg%n1I1E > 0) then
             
             ! Loop over ket 1I1E configurations
             do ikconf=cfg%off1I1E(kn),cfg%off1I1E(kn+1)-1

                ! Ket 1I1E configuration
                kconf_full=cfg%conf1I1E(:,:,ikconf)
                ksop_full=cfg%sop1I1E(:,:,ikconf)

                ! Number of creation and annihilation operators linking
                ! the ket 1I1E and bra 1-hole configurations
                nac=n_create_annihilate(kconf_full,bconf1h_full,n_int)

                ! Cycle if the the bra 1-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 1I1E configuration
                if (nac > 3) cycle

                ! Number of open shells in the ket 1I1E configuration
                knopen=sop_nopen(ksop_full,n_int)
                
                ! Number of ket 1I1E CSFs
                knsp=ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! 1I1E - 1I matrix contributions
                if (cfg%off1I(bn) /= cfg%off1I(bn+1)) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n1I,cfg%n1I1E,&       ! no. bra and ket confs
                        cfg%conf1I,cfg%sop1I,&    ! bra confs and SOPs
                        cfg%n1h,cfg%off1I,&       ! no. bra hole confs and offsets
                        cfg%csfs1I,cfg%csfs1I1E,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.false.)
                endif

                ! 1I1E - 1E matrix contributions
                if (cfg%off1E(bn) /= cfg%off1E(bn+1)) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n1E,cfg%n1I1E,&       ! no. bra and ket confs
                        cfg%conf1E,cfg%sop1E,&    ! bra confs and SOPs
                        cfg%n1h,cfg%off1E,&       ! no. bra hole confs and offsets
                        cfg%csfs1E,cfg%csfs1I1E,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.false.)
                endif
                
             enddo

          endif
                
       enddo
       
    enddo
          
    return
    
  end subroutine rdm_2h_1h

!######################################################################
! rdm_2h_2h: Calculation of the 2I-2I, 2I-2E, 2I-1I1E, 2E-2E, 2I-1I1E,
!            and 1I1E-1I1E contributions to the MRCI 1-RDMs
!######################################################################
  subroutine rdm_2h_2h(cfg,csfdim,nroots,vec,rho)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! No. CSFs
    integer(is), intent(in)  :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in)  :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! Eigenvectors
    real(dp), intent(in)     :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout)  :: rho(nmo,nmo,nroots)

    ! Number of open shells preceding each MO
    integer(is)              :: nbefore(nmo)

    ! Working arrays
    integer(ib), allocatable :: kconf_int(:,:)
    integer(ib), allocatable :: ksop_int(:,:)
    integer(ib)              :: kconf_full(n_int,2)
    integer(ib)              :: ksop_full(n_int,2)
    integer(ib)              :: bconf2h_full(n_int,2)
    integer(is), parameter   :: maxexci=1
    integer(is)              :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)              :: n_int_I,knopen,knsp,kn,bn,nac,nac1
    integer(is)              :: ikconf

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I

    allocate(kconf_int(n_int_I,2))
    allocate(ksop_int(n_int_I,2))
    kconf_int=0_ib
    ksop_int=0_ib

!----------------------------------------------------------------------
! 2H-2H contributions to the 1-RDMs
!----------------------------------------------------------------------
    ! Loop over ket 2-hole configurations
    do kn=1,cfg%n2h
    
       ! Loop over bra 2-hole configurations
       !do bn=kn,cfg%n2h
       do bn=1,cfg%n2h
       
          ! Bra 2-hole configuration in the full MO space
          bconf2h_full=0_ib
          bconf2h_full(1:n_int_I,:)=cfg%conf2h(1:n_int_I,:,bn)

          ! Number of creation and annihilation operators linking the
          ! bra and ket 2-hole configurations
          nac1=n_create_annihilate(cfg%conf2h(1:n_int_I,:,kn),&
               cfg%conf2h(1:n_int_I,:,bn),n_int_I)

          ! Cycle if the full configurations generated by the hole
          ! configurations cannot be interacting
          if (nac1 > 6) cycle

          !
          ! Ket: 2I
          ! Bra: 2I, 2E and 1I1E
          !
          if (cfg%n2I > 0) then
             
             ! Loop over ket 2I configurations
             do ikconf=cfg%off2I(kn),cfg%off2I(kn+1)-1

                ! Ket 2I configuration
                kconf_int=cfg%conf2I(1:n_int_I,:,ikconf)
                ksop_int=cfg%sop2I(1:n_int_I,:,ikconf)
                kconf_full=0_ib
                ksop_full=0_ib
                kconf_full(1:n_int_I,:)=kconf_int
                ksop_full(1:n_int_I,:)=ksop_int
          
                ! Number of creation and annihilation operators linking
                ! the ket 2I and bra 2-hole configurations
                nac=n_create_annihilate(kconf_int,&
                     cfg%conf2h(1:n_int_I,:,bn),n_int_I)
          
                ! Cycle if the the bra 2-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 2I configuration
                if (nac > 4) cycle

                ! Number of open shells in the ket 2I configuration
                knopen=sop_nopen(ksop_int(1:n_int_I,:),n_int_I)
                
                ! Number of ket 2I CSFs
                knsp=ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! 2I - 2I matrix elements
                if (nac <= 4 .and. bn >= kn &
                     .and. cfg%off2I(bn) /= cfg%off2I(bn+1)) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n2I,cfg%n2I,&       ! no. bra and ket confs
                        cfg%conf2I,cfg%sop2I,&  ! bra confs and SOPs
                        cfg%n2h,cfg%off2I,&     ! no. bra hole confs and offsets
                        cfg%csfs2I,cfg%csfs2I,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.true.)
                endif
          
                ! 2I - 2E matrix elements
                if (nac == 0 &
                     .and. cfg%n2E > 0 &
                     .and. cfg%off2E(bn) /= cfg%off2E(bn+1)) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n2E,cfg%n2I,&       ! no. bra and ket confs
                        cfg%conf2E,cfg%sop2E,&  ! bra confs and SOPs
                        cfg%n2h,cfg%off2E,&     ! no. bra hole confs and offsets
                        cfg%csfs2E,cfg%csfs2I,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.false.)
                endif
          
                ! 2I - 1I1E matrix elements
                if (nac <= 2 &
                     .and. cfg%n1I1E > 0 &
                     .and. cfg%off1I1E(bn) /= cfg%off1I1E(bn+1)) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n1I1E,cfg%n2I,&        ! no. bra and ket confs
                        cfg%conf1I1E,cfg%sop1I1E,& ! bra confs and SOPs
                        cfg%n2h,cfg%off1I1E,&      ! no. bra hole confs and offsets
                        cfg%csfs1I1E,cfg%csfs2I,&  ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.false.)
                endif
                
             enddo
                
          endif

          !
          ! Ket: 2E
          ! Bra: 2E and 1I1E
          !
          if (cfg%n2E > 0 .and. nac1 <= 5) then

             ! Loop over ket 2E configurations
             do ikconf=cfg%off2E(kn),cfg%off2E(kn+1)-1

                ! Ket 2E configuration
                kconf_full=cfg%conf2E(:,:,ikconf)
                ksop_full=cfg%sop2E(:,:,ikconf)
                kconf_int=0_ib
                ksop_int=0_ib
                kconf_int(1:n_int_I,:)=kconf_full(1:n_int_I,:)
                ksop_int(1:n_int_I,:)=ksop_full(1:n_int_I,:)
          
                ! Number of creation and annihilation operators linking
                ! the ket 2E and bra 2-hole configurations
                nac=n_create_annihilate(kconf_full,bconf2h_full,n_int)
          
                ! Cycle if the the bra 2-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 2E configuration
                if (nac > 4) cycle

                ! Number of open shells in the ket 2E configuration
                knopen=sop_nopen(ksop_full,n_int)
                
                ! Number of ket 2E CSFs
                knsp=ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! 2E - 2E matrix elements
                if (cfg%off2E(bn) /= cfg%off2E(bn+1) .and. bn >= kn &
                     .and.nac1 <= 2) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n2E,cfg%n2E,&       ! no. bra and ket confs
                        cfg%conf2E,cfg%sop2E,&  ! bra confs and SOPs
                        cfg%n2h,cfg%off2E,&     ! no. bra hole confs and offsets
                        cfg%csfs2E,cfg%csfs2E,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.true.)
                endif
          
                ! 2E - 1I1E matrix elements
                if (cfg%off1I1E(bn) /= cfg%off1I1E(bn+1) &
                     .and. cfg%n1I1E /= 0) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n1I1E,cfg%n2E,&        ! no. bra and ket confs
                        cfg%conf1I1E,cfg%sop1I1E,& ! bra confs and SOPs
                        cfg%n2h,cfg%off1I1E,&      ! no. bra hole confs and offsets
                        cfg%csfs1I1E,cfg%csfs2E,&  ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.false.)
                endif
                
             enddo

          endif

          !
          ! Ket: 1I1E
          ! Bra: 1I1E
          !
          if (cfg%n1I1E > 0) then
          
             ! Loop over ket 1I1E configurations
             do ikconf=cfg%off1I1E(kn),cfg%off1I1E(kn+1)-1
                
                ! Ket 1I1E configuration
                kconf_full=cfg%conf1I1E(:,:,ikconf)
                ksop_full=cfg%sop1I1E(:,:,ikconf)
                kconf_int=0_ib
                ksop_int=0_ib
                kconf_int(1:n_int_I,:)=kconf_full(1:n_int_I,:)
                ksop_int(1:n_int_I,:)=ksop_full(1:n_int_I,:)
          
                ! Number of creation and annihilation operators linking
                ! the ket 2E and bra 2-hole configurations
                nac=n_create_annihilate(kconf_full,bconf2h_full,n_int)
          
                ! Cycle if the the bra 2-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 1I1E configuration
                if (nac > 4) cycle
          
                ! Number of open shells in the ket 1I1E configuration
                knopen=sop_nopen(ksop_full,n_int)
                
                ! Number of ket 1I1E CSFs
                knsp=ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! 1I1E - 1I1E matrix elements
                if (cfg%off1I1E(bn) /= cfg%off1I1E(bn+1) .and. bn >= kn) then
                   call rdm_batch(&
                        bn,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                        cfg%n1I1E,cfg%n1I1E,&       ! no. bra and ket confs
                        cfg%conf1I1E,cfg%sop1I1E,&  ! bra confs and SOPs
                        cfg%n2h,cfg%off1I1E,&       ! no. bra hole confs and offsets
                        cfg%csfs1I1E,cfg%csfs1I1E,& ! bra and ket CSF offsets
                        csfdim,nroots,vec,rho,cfg%m2c,.true.)
                endif
                
             enddo

          endif
          
       enddo

    enddo
          
    return
    
  end subroutine rdm_2h_2h
    
!######################################################################
! rdm_batch: Computes all the contributions to the 1-RDMs from a
!            ket conf and a single class (1I, 2I, etc.) of confs
!            generated by a single bra hole conf
!######################################################################
  subroutine rdm_batch(bn,ikconf,kconf,ksop,knopen,knsp,knbefore,&
       nbconf,nkconf,bconfs,bsops,nh,boffset,bcsfs,kcsfs,&
       csfdim,nroots,vec,rho,m2c,same_class)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    ! Index of the bra hole configuration
    integer(is), intent(in) :: bn

    ! Index of the ket conf
    integer(is), intent(in) :: ikconf
    
    ! Dimensions
    integer(is), intent(in) :: csfdim,nroots,nbconf,nkconf
    
    ! Ket conf and SOP
    integer(ib), intent(in) :: kconf(n_int,2),ksop(n_int,2)
    integer(is), intent(in) :: knopen,knsp,knbefore(nmo)

    ! Bra confs and SOPs
    integer(ib), intent(in) :: bconfs(n_int,2,nbconf)
    integer(ib), intent(in) :: bsops(n_int,2,nbconf)
    
    ! Bra configuration offset array
    integer(is), intent(in) :: nh
    integer(is), intent(in) :: boffset(nh+1)

    ! CSF offsets
    integer(is), intent(in) :: bcsfs(nbconf+1),kcsfs(nkconf+1)

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: rho(nmo,nmo,nroots)

    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! Same bra and ket configuration class?
    logical, intent(in)     :: same_class
    
    ! Working arrays
    integer(ib)             :: bconf(n_int,2),bsop(n_int,2)
    integer(is), parameter  :: maxexci=1
    integer(is)             :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: ibconf,ista,ibcsf,ikcsf,bomega,komega
    integer(is)             :: nexci,bnopen,bnsp,i,a
    real(dp)                :: bcoe,kcoe,prod
    
    ! Loop over the bra confs generated by the hole conf
    do ibconf=boffset(bn),boffset(bn+1)-1

       ! Bra configuration in the full space
       bconf=bconfs(:,:,ibconf)
       bsop=bsops(:,:,ibconf)

       ! Compute the excitation degree between the two
       ! configurations
       nexci=exc_degree_conf(kconf,bconf,n_int)
       
       ! Cycle if the excitation degree is not equal to 1
       if (nexci /= 1) cycle

       ! Number of open shells in the bra configuration
       bnopen=sop_nopen(bsop,n_int)
    
       ! Number of bra CSFs
       bnsp=ncsfs(bnopen)

       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(kconf,bconf,n_int,hlist(1),plist(1),1)

       ! Get the spin-coupling coefficients
       spincp(1:knsp,1:bnsp)=spincp_coeff(knsp,bnsp,ksop,plist(1),&
            hlist(1),knopen,knbefore)

       ! Idices of the 1-RDM elements
       i=m2c(hlist(1))
       a=m2c(plist(1))

       ! Loop over roots
       do ista=1,nroots

          ! Loop over bra CSFs
          bomega=0
          do ibcsf=bcsfs(ibconf),bcsfs(ibconf+1)-1
             bomega=bomega+1
             bcoe=vec(ibcsf,ista)

             ! Loop over ket CSFs
             komega=0
             do ikcsf=kcsfs(ikconf),kcsfs(ikconf+1)-1
                komega=komega+1
                kcoe=vec(ikcsf,ista)

                ! Cycle duplicate CSF pairs
                if (same_class .and. ibcsf < ikcsf) cycle
                
                ! Contribution to the 1-RDM
                prod=kcoe*bcoe*spincp(komega,bomega)
                rho(i,a,ista)=rho(i,a,ista)+prod
                rho(a,i,ista)=rho(a,i,ista)+prod
                
             enddo
             
          enddo
          
       enddo
       
    enddo
       
    return
    
  end subroutine rdm_batch
    
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
