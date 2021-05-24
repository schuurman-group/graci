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
