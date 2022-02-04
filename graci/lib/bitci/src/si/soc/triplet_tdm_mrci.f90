!**********************************************************************
! Routines for the calculation of triplet TDMS for MRCI wavefunctions
!**********************************************************************
module triplet_tdm

  use constants

  implicit none

  private

  public :: triplet_tdm_mrci
  
  ! Triplet spin-coupling coefficients
  real(dp), allocatable, private :: scp(:)

contains

!######################################################################
! triplet_tdm_mrci: Master routine for the calculation of MRCI triplet
!                   TDMs for a specified set of bra and ket states
!######################################################################
!                   Note that the TDMS have the 'Canonical' MO indexing
!######################################################################
  subroutine triplet_tdm_mrci(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,&
       nvecK,vecB,vecK,npairs,Tij,Bmap,Kmap)

    use constants
    use bitglobal
    use conftype
    use iomod
    use timing
    use utils
    
    implicit none

    ! Component of the triplet spin tensor operator
    integer(is), intent(in) :: kval
    
    ! MRCI configuration derived types
    type(mrcfg), intent(in) :: cfgB,cfgK

    ! Dimensions
    integer(is), intent(in) :: csfdimB,csfdimK,nvecB,nvecK,npairs

    ! Eigenvectors
    real(dp), intent(in)    :: vecB(csfdimB,nvecB)
    real(dp), intent(in)    :: vecK(csfdimK,nvecK)

    ! Triplet TDMs
    real(dp), intent(out)   :: Tij(nmo,nmo,npairs)

    ! Bra-ket pair to eigenvector mapping arrays
    integer(is), intent(in) :: Bmap(npairs),Kmap(npairs)

    ! Everything else
    integer(is)             :: arrdim
    
    ! Timing variables
    real(dp)                :: tcpu_start,tcpu_end,twall_start,&
                               twall_end

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)

!----------------------------------------------------------------------
! Sanity check on the triplet spin tensor operator component
!----------------------------------------------------------------------
    if (kval /= 0 .and. kval /= 1) then
       errmsg='Error in triplet_tdm_mrc: nonsensical triplet spin ' &
            //'tensor operator component'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    Tij=0.0d0

    arrdim=maxval(ncsfsB(0:nomax))*maxval(ncsfsK(0:nomax))
    allocate(scp(arrdim))
    scp=0.0d0

!----------------------------------------------------------------------
! (1) Ref - Ref contributions to the triplet TDMs
!----------------------------------------------------------------------
    call Tij_0h_0h(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,&
         vecK,npairs,Tij,Bmap,Kmap)
    
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'tdm_mrci')

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(scp)
    
    STOP
    
    return
    
  end subroutine triplet_tdm_mrci
  
!######################################################################
! Tij_0h_0h: Calculation of the Ref-Ref contributions to the MRCI
!            triplet TDMs
!######################################################################
  subroutine Tij_0h_0h(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,&
       vecB,vecK,npairs,Tij,Bmap,Kmap)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! Component of the triplet spin tensor operator
    integer(is), intent(in) :: kval
    
    ! MRCI configuration derived types
    type(mrcfg), intent(in) :: cfgB,cfgK

    ! Dimensions
    integer(is), intent(in) :: csfdimB,csfdimK,nvecB,nvecK,npairs

    ! Eigenvectors
    real(dp), intent(in)    :: vecB(csfdimB,nvecB)
    real(dp), intent(in)    :: vecK(csfdimK,nvecK)

    ! 1-TDMs
    real(dp), intent(inout) :: Tij(nmo,nmo,npairs)

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
    integer(is)             :: ikcsf,ibcsf
    integer(is)             :: counter
    real(dp)                :: kcoe,bcoe
    real(dp)                :: prod

!----------------------------------------------------------------------
! Contributions from bra and ket reference space CSFs
!----------------------------------------------------------------------
    n_int_I=cfgK%n_int_I

    ! Loop over ket configurations
    do ikconf=1,cfgK%n0h

       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_I,:)=cfgK%conf0h(:,:,ikconf)
       ksop_full(1:n_int_I,:)=cfgK%sop0h(:,:,ikconf)

       ! Number of open shells in the ket configuration
       knopen=sop_nopen(cfgK%sop0h(:,:,ikconf),n_int_I)

       ! Number of ket CSFs
       knsp=ncsfsK(knopen)
       
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
          bnsp=ncsfsB(bnopen)

          ! Get the indices of the MOs involved in the excitation
          hlist=0
          plist=0
          call get_exci_indices(cfgK%conf0h(:,:,ikconf),&
               cfgB%conf0h(:,:,ibconf),n_int_I,hlist(1:nexci),&
               plist(1:nexci),nexci)

          ! Get the spin-coupling coefficients
          ! k <-> kval !
          scp(1:knsp*bnsp)=triplet_scc(kval,knsp,bnsp,ksop_full,&
               plist(1),hlist(1),knopen,nbefore)

          ! Idices of the TDM elements
          i=cfgB%m2c(hlist(1))
          a=cfgB%m2c(plist(1))

          ! Loop over TDMs
          do ipair=1,npairs

             ! Bra and ket eigenvector indices
             Bindx=Bmap(ipair)
             Kindx=Kmap(ipair)

             ! Initialise the spin-coupling coefficient counter
             counter=0

             ! Loop over ket CSFs
             do ikcsf=cfgK%csfs0h(ikconf),cfgK%csfs0h(ikconf+1)-1
                kcoe=vecK(ikcsf,Kindx)
                
                ! Loop over bra CSFs
                do ibcsf=cfgB%csfs0h(ibconf),cfgB%csfs0h(ibconf+1)-1
                   bcoe=vecB(ibcsf,Bindx)
                   counter=counter+1
                   
                   ! Contribution to the 1-TDM
                   prod=kcoe*bcoe*scp(counter)
                   Tij(a,i,ipair)=Tij(a,i,ipair)+prod
                   
                enddo
                
             enddo
             
          enddo
             
       enddo
          
    enddo
    
    return
    
  end subroutine Tij_0h_0h
    
!######################################################################
! triplet_scc: Given a SOP and pair of creation/annihilation operator
!              indices, returns the complete set of triplet
!              spin-coupling coefficients for the kval'th component
!              of the triplet spin tensor operator
!######################################################################
  function triplet_scc(kval,knsp,bnsp,sop,ac,ia,nopen,nbefore)

    use constants
    use bitglobal
    use pattern_indices, only: get_icase
    use pattern_indices_k1
    use bitstrings
    use iomod
    
    implicit none

    ! Component of the triplet spin tensor operator
    integer(is), intent(in) :: kval
    
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
    real(dp)                :: triplet_scc(knsp*bnsp)
    
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
    select case(kval)
    case(0)
       print*,''
       print*,'pattern_index_k0 needs writing'
       stop
    case(1)
       pattern=pattern_index_k1(sop,ac,ia,nc,na,nopen,icase)
    end select
    
!----------------------------------------------------------------------
! Fill in the array of spin-coupling coefficients
!----------------------------------------------------------------------
    triplet_scc=spincp(pattern:pattern+knsp*bnsp-1)
    
    return

  end function triplet_scc
  
!######################################################################
  
end module triplet_tdm
