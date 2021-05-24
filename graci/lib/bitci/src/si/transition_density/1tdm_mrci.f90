!**********************************************************************
! Routines for the calculation of 1-TDMs for MRCI wavefunctions
!**********************************************************************
module tdm

  use constants
  
  implicit none

  ! Spin-coupling coefficients
  real(dp), allocatable, private :: spincp(:,:)
  
contains

!######################################################################
! 1tdm_mrci: Master routine for the calculation of MRCI 1-TDMs for
!            a specified set of bra and ket states
!######################################################################
!            Note that the 1-TDMS have the 'Canonical' MO indexing
!######################################################################
  subroutine tdm_mrci(cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,&
       vecK,npairs,rhoij,Bmap,Kmap)

    use constants
    use bitglobal
    use conftype
    use iomod
    use timing
    use utils
    
    implicit none

    ! MRCI configuration derived types
    type(mrcfg), intent(in) :: cfgB,cfgK

    ! Dimensions
    integer(is), intent(in) :: csfdimB,csfdimK,nvecB,nvecK,npairs

    ! Eigenvectors
    real(dp), intent(in)    :: vecB(csfdimB,nvecB)
    real(dp), intent(in)    :: vecK(csfdimK,nvecK)

    ! 1-TDMs
    real(dp), intent(out)   :: rhoij(nmo,nmo,npairs)

    ! Bra-ket pair to eigenvector mapping arrays
    integer(is), intent(in) :: Bmap(npairs),Kmap(npairs)
    
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
    rhoij=0.0d0

    allocate(spincp(ncsfs(nomax),ncsfs(nomax)))
    spincp=0.0d0

!----------------------------------------------------------------------
! (1) Ref - Ref contributions to the 1-TDMs
!----------------------------------------------------------------------
    call tdm_0h_0h(cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,vecK,&
         npairs,rhoij,Bmap,Kmap)

!----------------------------------------------------------------------
! (2) Ref - 1H contributions to the 1-RDMs
!----------------------------------------------------------------------
    call tdm_0h_1h(cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,vecK,&
         npairs,rhoij,Bmap,Kmap)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(spincp)
    
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'tdm_mrci')
    
    return
    
  end subroutine tdm_mrci

!######################################################################
! tdm_mrci_0h_0h: Calculation of the Ref-Ref contributions to the MRCI
!                 1-TDMs
!######################################################################
  subroutine tdm_0h_0h(cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,&
       vecK,npairs,rhoij,Bmap,Kmap)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! MRCI configuration derived types
    type(mrcfg), intent(in) :: cfgB,cfgK

    ! Dimensions
    integer(is), intent(in) :: csfdimB,csfdimK,nvecB,nvecK,npairs

    ! Eigenvectors
    real(dp), intent(in)    :: vecB(csfdimB,nvecB)
    real(dp), intent(in)    :: vecK(csfdimK,nvecK)

    ! 1-TDMs
    real(dp), intent(inout) :: rhoij(nmo,nmo,npairs)

    ! Bra-ket pair to eigenvector mapping arrays
    integer(is), intent(in) :: Bmap(npairs),Kmap(npairs)
    
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
    integer(is)             :: i,a,ipair,Bindx,Kindx
    integer(is)             :: ikcsf,ibcsf,komega,bomega
    real(dp)                :: kcoe,bcoe
    real(dp)                :: prod

!----------------------------------------------------------------------
! Contributions from bra and ket reference space CSFs
!----------------------------------------------------------------------
    n_int_I=cfgK%n_int_I

    ! Loop over ket configurations
    do ikconf=1,cfgK%n0h-1

       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_I,:)=cfgK%conf0h(:,:,ikconf)
       ksop_full(1:n_int_I,:)=cfgK%sop0h(:,:,ikconf)

       ! Number of open shells in the ket configuration
       knopen=sop_nopen(cfgK%sop0h(:,:,ikconf),n_int_I)

       ! Number of ket CSFs
       knsp=ncsfs(knopen)
       
       ! Get the number of open shells preceding each ket conf MO
       call nobefore(ksop_full,nbefore)

       ! Loop over bra configurations
       do ibconf=1,cfgB%n0h

          ! Compute the excitation degree between the two
          ! configurations
          nexci=exc_degree_conf(cfgK%conf0h(:,:,ikconf),&
               cfgB%conf0h(:,:,ibconf),n_int_I)

          ! Cycle if the excitation degree is not equal to 1
          if (nexci /= 1) cycle

          ! Number of open shells in the bra configuration
          bnopen=sop_nopen(cfgB%sop0h(:,:,ibconf),n_int_I)

          ! Number of bra CSFs
          bnsp=ncsfs(bnopen)

          ! Get the indices of the MOs involved in the excitation
          hlist=0
          plist=0
          call get_exci_indices(cfgK%conf0h(:,:,ikconf),&
               cfgB%conf0h(:,:,ibconf),n_int_I,hlist(1:nexci),&
               plist(1:nexci),nexci)

          ! Get the spin-coupling coefficients
          spincp(1:knsp,1:bnsp)=spincp_coeff(knsp,bnsp,ksop_full,&
               plist(1),hlist(1),knopen,nbefore)

          ! Idices of the 1-TDM elements
          i=cfgB%m2c(hlist(1))
          a=cfgB%m2c(plist(1))

          ! Loop over 1-TDMs
          do ipair=1,npairs

             ! Bra and ket eigenvector indices
             Bindx=Bmap(ipair)
             Kindx=Kmap(ipair)

             ! Loop over bra CSFs
             bomega=0
             do ibcsf=cfgB%csfs0h(ibconf),cfgB%csfs0h(ibconf+1)-1
                bomega=bomega+1
                bcoe=vecB(ibcsf,Bindx)
                
                ! Loop over ket CSFs
                komega=0
                do ikcsf=cfgK%csfs0h(ikconf),cfgK%csfs0h(ikconf+1)-1
                   komega=komega+1
                   kcoe=vecK(ikcsf,Kindx)

                   ! Contribution to the 1-TDM
                   prod=kcoe*bcoe*spincp(komega,bomega)
                   rhoij(a,i,ipair)=rhoij(a,i,ipair)+prod
                   
                enddo
                
             enddo
             
          enddo
          
       enddo
          
    enddo
       
    return
    
  end subroutine tdm_0h_0h

!######################################################################
! tdm_mrci_0h_1h: Calculation of the Ref-1I and Ref-1E contributions
!                 to the MRCI 1-TDMs
!######################################################################
  subroutine tdm_0h_1h(cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,&
       vecK,npairs,rhoij,Bmap,Kmap)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! MRCI configuration derived types
    type(mrcfg), intent(in) :: cfgB,cfgK

    ! Dimensions
    integer(is), intent(in) :: csfdimB,csfdimK,nvecB,nvecK,npairs

    ! Eigenvectors
    real(dp), intent(in)    :: vecB(csfdimB,nvecB)
    real(dp), intent(in)    :: vecK(csfdimK,nvecK)

    ! 1-TDMs
    real(dp), intent(inout) :: rhoij(nmo,nmo,npairs)

    ! Bra-ket pair to eigenvector mapping arrays
    integer(is), intent(in) :: Bmap(npairs),Kmap(npairs)

    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)

    ! Working arrays
    integer(ib)             :: kconf_full(n_int,2)
    integer(ib)             :: ksop_full(n_int,2)
    integer(ib)             :: bconf_full(n_int,2)
    integer(ib)             :: bsop_full(n_int,2)
    
    ! Everything else
    integer(is)             :: ikconf,ibconf,n,nac,nexci,n_int_I
    integer(is)             :: knsp,knopen,bnsp,bnopen

!----------------------------------------------------------------------
! Contributions of the ket reference and bra 1I & 1E CSFs
!----------------------------------------------------------------------
    n_int_I=cfgK%n_int_I

    ! Loop over ket reference configurations
    do ikconf=1,cfgK%n0h

       ! Number of open shells in the ket configuration
       knopen=sop_nopen(cfgK%sop0h(:,:,ikconf),n_int_I)
       
       ! Number of ket CSFs
       knsp=ncsfs(knopen)

       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_I,:)=cfgK%conf0h(:,:,ikconf)
       ksop_full(1:n_int_I,:)=cfgK%sop0h(:,:,ikconf)
       
       ! Get the number of open shells preceding each ket conf MO
       call nobefore(ksop_full,nbefore)

       ! Loop over bra 1-hole configurations
       do n=1,cfgB%n1h

          ! Number of creation and annihilation operators linking the
          ! reference and 1-hole configurations
          nac=n_create_annihilate(cfgK%conf0h(1:n_int_I,:,ikconf), &
               cfgB%conf1h(1:n_int_I,:,n),n_int_I)

          ! Ket Ref - bra 1I contributions
          if (nac <= 3 .and. cfgB%n1I > 0 &
               .and. cfgB%off1I(n) /= cfgB%off1I(n+1)) then
             call tdm_batch(&
                  n,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                  cfgB%n1I,cfgK%n0h,&       ! no. bra and ket confs
                  cfgB%conf1I,cfgB%sop1I,&  ! bra confs and SOPs
                  cfgB%n1h,cfgB%off1I,&     ! no. bra hole confs and offsets
                  cfgB%csfs1I,cfgK%csfs0h,& ! bra and ket CSF offsets
                  csfdimB,csfdimK,nvecB,nvecK,&
                  vecB,vecK,npairs,rhoij,Bmap,Kmap,&
                  cfgB%m2c,.false.)
          endif

          ! Ket Ref - bra 1E contributions
          if (nac <= 1 .and. cfgB%n1E > 0 &
               .and. cfgB%off1E(n) /= cfgB%off1E(n+1)) then
             call tdm_batch(&
                  n,ikconf,kconf_full,ksop_full,knopen,knsp,nbefore,&
                  cfgB%n1E,cfgK%n0h,&       ! no. bra and ket confs
                  cfgB%conf1E,cfgB%sop1E,&  ! bra confs and SOPs
                  cfgB%n1h,cfgB%off1E,&     ! no. bra hole confs and offsets
                  cfgB%csfs1E,cfgK%csfs0h,& ! bra and ket CSF offsets
                  csfdimB,csfdimK,nvecB,nvecK,&
                  vecB,vecK,npairs,rhoij,Bmap,Kmap,&
                  cfgB%m2c,.false.)
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Contributions of the ket 1I & 1E and bra reference CSFs
!----------------------------------------------------------------------
    
    return
    
  end subroutine tdm_0h_1h

!######################################################################
! tdm_batch: Computes all the contributions to the 1-TDMs from a
!            ket conf and a single class (1I, 2I, etc.) of confs
!            generated by a single bra hole conf
!######################################################################
  subroutine tdm_batch(bn,ikconf,kconf,ksop,knopen,knsp,knbefore,&
       nbconf,nkconf,bconfs,bsops,nh,boffset,bcsfs,kcsfs,&
       csfdimB,csfdimK,nvecB,nvecK,vecB,vecK,npairs,rhoij,Bmap,Kmap,&
       m2c,same_class)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    ! Index of the bra hole configuration
    integer(is), intent(in) :: bn

    ! Index of the ket conf
    integer(is), intent(in) :: ikconf
    
    ! Dimensions
    integer(is), intent(in) :: csfdimB,csfdimK,nvecB,nvecK,npairs
    integer(is), intent(in) :: nbconf,nkconf
    
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
    real(dp), intent(in)    :: vecB(csfdimB,nvecB)
    real(dp), intent(in)    :: vecK(csfdimK,nvecK)
    
    ! 1-TDMs
    real(dp), intent(inout) :: rhoij(nmo,nmo,npairs)

    ! Bra-ket pair to eigenvector mapping arrays
    integer(is), intent(in) :: Bmap(npairs),Kmap(npairs)
    
    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! Same bra and ket configuration class?
    logical, intent(in)     :: same_class
    
    ! Working arrays
    integer(ib)             :: bconf(n_int,2),bsop(n_int,2)
    integer(is), parameter  :: maxexci=1
    integer(is)             :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: ibconf,ipair,ibcsf,ikcsf,bomega,komega
    integer(is)             :: Bindx,Kindx
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

       ! Loop over 1-TDMs
       do ipair=1,npairs

          ! Bra and ket eigenvector indices
          Bindx=Bmap(ipair)
          Kindx=Kmap(ipair)
          
          ! Loop over bra CSFs
          bomega=0
          do ibcsf=bcsfs(ibconf),bcsfs(ibconf+1)-1
             bomega=bomega+1
             bcoe=vecB(ibcsf,Bindx)

             ! Loop over ket CSFs
             komega=0
             do ikcsf=kcsfs(ikconf),kcsfs(ikconf+1)-1
                komega=komega+1
                kcoe=vecK(ikcsf,Kindx)

                ! Cycle duplicate CSF pairs
                if (same_class .and. ibcsf < ikcsf) cycle
                
                ! Contribution to the 1-TDM
                prod=kcoe*bcoe*spincp(komega,bomega)
                rhoij(a,i,ipair)=rhoij(a,i,ipair)+prod
                                
             enddo
             
          enddo
          
       enddo
       
    enddo
       
    return
    
  end subroutine tdm_batch
    
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
  
end module tdm
