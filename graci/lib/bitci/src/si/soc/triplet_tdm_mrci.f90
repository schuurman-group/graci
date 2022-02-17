!**********************************************************************
! Routines for the calculation of triplet TDMS for MRCI wavefunctions
!**********************************************************************
module triplet_tdm

  use constants

  implicit none

  private

  public :: triplet_tdm_mrci
  
  ! Triplet spin-coupling coefficients
  real(dp), allocatable, private :: scc(:)

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

    arrdim=maxval(cfgB%ncsfs(0:nomax))*maxval(cfgK%ncsfs(0:nomax))
    allocate(scc(arrdim))
    scc=0.0d0

!----------------------------------------------------------------------
! (0) On-diagonal elements T_ii
!----------------------------------------------------------------------
    call Tii_all(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,&
         vecK,npairs,Tij,Bmap,Kmap)
    
!----------------------------------------------------------------------
! (1) Ref - Ref contributions to the triplet TDMs
!----------------------------------------------------------------------
    call Tij_0h_0h(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,&
         vecK,npairs,Tij,Bmap,Kmap)

!----------------------------------------------------------------------
! (2) Ref - 1H contributions to the triplet TDMs
!----------------------------------------------------------------------
    call Tij_0h_1h(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,&
         vecK,npairs,Tij,Bmap,Kmap)

!----------------------------------------------------------------------
! (3) 1H - 1H contributions to the triplet TDMs
!----------------------------------------------------------------------
    call Tij_1h_1h(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,&
         vecK,npairs,Tij,Bmap,Kmap)

!----------------------------------------------------------------------
! (4) 2H - 1H contributions to the triplet TDMs
!----------------------------------------------------------------------
    call Tij_2h_1h(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,&
         vecK,npairs,Tij,Bmap,Kmap)

!----------------------------------------------------------------------
! (5) 2H - 2H contributions to the triplet TDMs
!----------------------------------------------------------------------
    call Tij_2h_2h(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,vecB,&
         vecK,npairs,Tij,Bmap,Kmap)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(scc)
    
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'triplet_tdm_mrci')

    return
    
  end subroutine triplet_tdm_mrci

!######################################################################
! Tii_all: Calculation of all contributions to the on-diagonal
!          elements of the MRCI triplet TDMs
!######################################################################
! Note that these elements are only non-zero for the k=0 component of
! triplet spin tensor operator
!######################################################################
  subroutine Tii_all(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,&
       vecB,vecK,npairs,Tij,Bmap,Kmap)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use dethash
    
    implicit none

    ! Component of the triplet spin tensor operator
    integer(is), intent(in)  :: kval
    
    ! MRCI configuration derived types
    type(mrcfg), intent(in)  :: cfgB,cfgK

    ! Dimensions
    integer(is), intent(in)  :: csfdimB,csfdimK,nvecB,nvecK,npairs

    ! Eigenvectors
    real(dp), intent(in)     :: vecB(csfdimB,nvecB)
    real(dp), intent(in)     :: vecK(csfdimK,nvecK)

    ! 1-TDMs
    real(dp), intent(inout)  :: Tij(nmo,nmo,npairs)

    ! Bra-ket pair to eigenvector mapping arrays
    integer(is), intent(in)  :: Bmap(npairs),Kmap(npairs)

    ! Hash table
    type(dhtbl)              :: h
    integer(is)              :: initial_size
    integer(ib), allocatable :: key(:,:)
        
    ! Open shell MO indices
    integer(is)              :: socc(nmo)
    integer(is)              :: nsocc
    
    ! Everything else
    integer(is)              :: n_int_I
    integer(is)              :: n,ibconf,ikconf,nopen,ibcsf,ikcsf
    integer(is)              :: iopen,i,i1,ipair,Bindx,Kindx,ioff
    real(dp)                 :: kcoe,bcoe,prod
    
!----------------------------------------------------------------------
! Return if k != 0. This could be checked in the calling routine, but
! this seems like a nice safety measure
!----------------------------------------------------------------------
    if (kval /= 0) return

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfgK%n_int_I

    allocate(key(n_int,2))
    
!----------------------------------------------------------------------
! Ref contributions
!----------------------------------------------------------------------
    ! Initialise the hash table
    initial_size=2**15
    call h%initialise_table(initial_size)
    
    ! Insert the ket ref confs and their indices into the hash table
    key=0_ib
    do ikconf=1,cfgK%n0h
       key(1:n_int_I,:)=cfgK%conf0h(:,:,ikconf)
       call h%insert_key(key,ikconf)
    enddo

    ! Loop over bra ref confs
    do ibconf=1,cfgB%n0h

       ! Is this also a ket ref conf?
       key(1:n_int_I,:)=cfgB%conf0h(:,:,ibconf)
       if (h%key_exists(key)) then

          ! Retrieve the ket conf index
          ikconf=h%get_value(key)

          ! Number of open shells
          nopen=sop_nopen(cfgB%sop0h(:,:,ibconf),n_int_I)
          
          ! Cycle if there are no open shells
          if (nopen == 0) cycle
          
          ! Open shell indices
          call sop_socc_list(cfgB%sop0h(:,:,ibconf),n_int_I,socc,&
               nmo,nsocc)

          ! Compute the contributions to the TDMs
          call ondiag_batch(ibconf,ikconf,&
               cfgB%n0h,cfgK%n0h,&
               cfgB%csfs0h,cfgK%csfs0h,&
               csfdimB,csfdimK,nvecB,nvecK,&
               vecB,vecK,npairs,Bmap,Kmap,nopen,socc,cfgB%m2c,Tij)
          
       endif
       
    enddo
       
    ! Delete the hash table
    call h%delete_table

!----------------------------------------------------------------------
! 1I contributions
!----------------------------------------------------------------------
    if (cfgK%n1I > 0 .and. cfgB%n1I > 0) then
           
       ! Initialise the hash table
       initial_size=max(2**15,cfgK%n1I)
       call h%initialise_table(initial_size)
       
       ! Insert the ket 1I confs and their indices into the hash table
       do n=1,cfgK%n1h
          do ikconf=cfgK%off1I(n),cfgK%off1I(n+1)-1
             key=cfgK%conf1I(:,:,ikconf)
             call h%insert_key(key,ikconf)
          enddo
       enddo
       
       ! Loop over bra 1I confs
       do n=1,cfgB%n1h
          do ibconf=cfgB%off1I(n),cfgB%off1I(n+1)-1

             ! Is this also a ket 1I conf?
             key=cfgB%conf1I(:,:,ibconf)
             if (h%key_exists(key)) then

                ! Retrieve the ket conf index
                ikconf=h%get_value(key)

                ! Number of open shells
                nopen=sop_nopen(cfgB%sop1I(:,:,ibconf),n_int)
                
                ! Cycle if there are no open shells
                if (nopen == 0) cycle
                
                ! Open shell indices
                call sop_socc_list(cfgB%sop1I(:,:,ibconf),n_int,socc,&
                     nmo,nsocc)

                ! Compute the contributions to the TDMs
                call ondiag_batch(ibconf,ikconf,&
                     cfgB%n1I,cfgK%n1I,&
                     cfgB%csfs1I,cfgK%csfs1I,&
                     csfdimB,csfdimK,nvecB,nvecK,&
                     vecB,vecK,npairs,Bmap,Kmap,nopen,socc,cfgB%m2c,Tij)
                
             endif
                
          enddo
       enddo
          
       ! Delete the hash table
       call h%delete_table

    endif

!----------------------------------------------------------------------
! 1E contributions
!----------------------------------------------------------------------
    if (cfgK%n1E > 0 .and. cfgB%n1E > 0) then

       ! Initialise the hash table
       initial_size=max(2**15,cfgK%n1E)
       call h%initialise_table(initial_size)
       
       ! Insert the ket 1E confs and their indices into the hash table
       do n=1,cfgK%n1h
          do ikconf=cfgK%off1E(n),cfgK%off1E(n+1)-1
             key=cfgK%conf1E(:,:,ikconf)
             call h%insert_key(key,ikconf)
          enddo
       enddo

       ! Loop over bra 1E confs
       do n=1,cfgB%n1h
          do ibconf=cfgB%off1E(n),cfgB%off1E(n+1)-1

             ! Is this also a ket 1E conf?
             key=cfgB%conf1E(:,:,ibconf)
             if (h%key_exists(key)) then

                ! Retrieve the ket conf index
                ikconf=h%get_value(key)

                ! Number of open shells
                nopen=sop_nopen(cfgB%sop1E(:,:,ibconf),n_int)
                
                ! Cycle if there are no open shells
                if (nopen == 0) cycle
                
                ! Open shell indices
                call sop_socc_list(cfgB%sop1E(:,:,ibconf),n_int,socc,&
                     nmo,nsocc)

                ! Compute the contributions to the TDMs
                call ondiag_batch(ibconf,ikconf,&
                     cfgB%n1E,cfgK%n1E,&
                     cfgB%csfs1E,cfgK%csfs1E,&
                     csfdimB,csfdimK,nvecB,nvecK,&
                     vecB,vecK,npairs,Bmap,Kmap,nopen,socc,cfgB%m2c,Tij)
                
             endif
                
          enddo
       enddo
          
       ! Delete the hash table
       call h%delete_table

    endif

!----------------------------------------------------------------------
! 2I contributions
!----------------------------------------------------------------------
    if (cfgK%n2I > 0 .and. cfgB%n2I > 0) then
           
       ! Initialise the hash table
       initial_size=max(2**15,cfgK%n2I)
       call h%initialise_table(initial_size)
       
       ! Insert the ket 2I confs and their indices into the hash table
       do n=1,cfgK%n2h
          do ikconf=cfgK%off2I(n),cfgK%off2I(n+1)-1
             key=cfgK%conf2I(:,:,ikconf)
             call h%insert_key(key,ikconf)
          enddo
       enddo
       
       ! Loop over bra 2I confs
       do n=1,cfgB%n2h
          do ibconf=cfgB%off2I(n),cfgB%off2I(n+1)-1

             ! Is this also a ket 2I conf?
             key=cfgB%conf2I(:,:,ibconf)
             if (h%key_exists(key)) then

                ! Retrieve the ket conf index
                ikconf=h%get_value(key)

                ! Number of open shells
                nopen=sop_nopen(cfgB%sop2I(:,:,ibconf),n_int)
                
                ! Cycle if there are no open shells
                if (nopen == 0) cycle
                
                ! Open shell indices
                call sop_socc_list(cfgB%sop2I(:,:,ibconf),n_int,socc,&
                     nmo,nsocc)

                ! Compute the contributions to the TDMs
                call ondiag_batch(ibconf,ikconf,&
                     cfgB%n2I,cfgK%n2I,&
                     cfgB%csfs2I,cfgK%csfs2I,&
                     csfdimB,csfdimK,nvecB,nvecK,&
                     vecB,vecK,npairs,Bmap,Kmap,nopen,socc,cfgB%m2c,Tij)
                
             endif
                
          enddo
       enddo
          
       ! Delete the hash table
       call h%delete_table

    endif

!----------------------------------------------------------------------
! 2E contributions
!----------------------------------------------------------------------
    if (cfgK%n2E > 0 .and. cfgB%n2E > 0) then
           
       ! Initialise the hash table
       initial_size=max(2**15,cfgK%n2E)
       call h%initialise_table(initial_size)
       
       ! Insert the ket 2E confs and their indices into the hash table
       do n=1,cfgK%n2h
          do ikconf=cfgK%off2E(n),cfgK%off2E(n+1)-1
             key=cfgK%conf2E(:,:,ikconf)
             call h%insert_key(key,ikconf)
          enddo
       enddo
       
       ! Loop over bra 2E confs
       do n=1,cfgB%n2h
          do ibconf=cfgB%off2E(n),cfgB%off2E(n+1)-1

             ! Is this also a ket 2E conf?
             key=cfgB%conf2E(:,:,ibconf)
             if (h%key_exists(key)) then

                ! Retrieve the ket conf index
                ikconf=h%get_value(key)

                ! Number of open shells
                nopen=sop_nopen(cfgB%sop2E(:,:,ibconf),n_int)
                
                ! Cycle if there are no open shells
                if (nopen == 0) cycle
                
                ! Open shell indices
                call sop_socc_list(cfgB%sop2E(:,:,ibconf),n_int,socc,&
                     nmo,nsocc)

                ! Compute the contributions to the TDMs
                call ondiag_batch(ibconf,ikconf,&
                     cfgB%n2E,cfgK%n2E,&
                     cfgB%csfs2E,cfgK%csfs2E,&
                     csfdimB,csfdimK,nvecB,nvecK,&
                     vecB,vecK,npairs,Bmap,Kmap,nopen,socc,cfgB%m2c,Tij)
                
             endif
                
          enddo
       enddo
          
       ! Delete the hash table
       call h%delete_table

    endif

!----------------------------------------------------------------------
! 1I1E contributions
!----------------------------------------------------------------------
    if (cfgK%n1I1E > 0 .and. cfgB%n1I1E > 0) then
           
       ! Initialise the hash table
       initial_size=max(2**15,cfgK%n1I1E)
       call h%initialise_table(initial_size)
       
       ! Insert the ket 1I1E confs and their indices into the hash
       ! table
       do n=1,cfgK%n2h
          do ikconf=cfgK%off1I1E(n),cfgK%off1I1E(n+1)-1
             key=cfgK%conf1I1E(:,:,ikconf)
             call h%insert_key(key,ikconf)
          enddo
       enddo
       
       ! Loop over bra 1I1E confs
       do n=1,cfgB%n2h
          do ibconf=cfgB%off1I1E(n),cfgB%off1I1E(n+1)-1

             ! Is this also a ket 1I1E conf?
             key=cfgB%conf1I1E(:,:,ibconf)
             if (h%key_exists(key)) then

                ! Retrieve the ket conf index
                ikconf=h%get_value(key)

                ! Number of open shells
                nopen=sop_nopen(cfgB%sop1I1E(:,:,ibconf),n_int)
                
                ! Cycle if there are no open shells
                if (nopen == 0) cycle
                
                ! Open shell indices
                call sop_socc_list(cfgB%sop1I1E(:,:,ibconf),n_int,&
                     socc,nmo,nsocc)

                ! Compute the contributions to the TDMs
                call ondiag_batch(ibconf,ikconf,&
                     cfgB%n1I1E,cfgK%n1I1E,&
                     cfgB%csfs1I1E,cfgK%csfs1I1E,&
                     csfdimB,csfdimK,nvecB,nvecK,&
                     vecB,vecK,npairs,Bmap,Kmap,nopen,socc,cfgB%m2c,Tij)
                
             endif
                
          enddo
       enddo
          
       ! Delete the hash table
       call h%delete_table

    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(key)
    
    return
    
  end subroutine Tii_all

!######################################################################
! ondiag_batch: Computes a batch of on-diagonal triplet TDM matrix
!               elements for a given pair of (identical) bra and ket
!               confs
!######################################################################
  subroutine ondiag_batch(ibconf,ikconf,nbconf,nkconf,bcsfs,kcsfs,&
       csfdimB,csfdimK,nvecB,nvecK,vecB,vecK,npairs,Bmap,Kmap,nopen,&
       socc,m2c,Tij)

    use constants
    use bitglobal
    
    implicit none

    ! Index of the bra and ket conf
    integer(is), intent(in) :: ibconf,ikconf

    ! Bra and ket CSF offsets
    integer(is), intent(in) :: nbconf,nkconf
    integer(is), intent(in) :: bcsfs(nbconf+1),kcsfs(nkconf+1)

    ! Eigenvectors
    integer(is), intent(in) :: csfdimB,csfdimK,nvecB,nvecK
    real(dp), intent(in)    :: vecB(csfdimB,nvecB)
    real(dp), intent(in)    :: vecK(csfdimK,nvecK)

    ! Bra-ket pair to eigenvector mapping arrays
    integer(is), intent(in) :: npairs
    integer(is), intent(in) :: Bmap(npairs),Kmap(npairs)

    ! Open shells
    integer(is), intent(in) :: nopen
    integer(is), intent(in) :: socc(nmo)
    
    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! Triplet TDMs
    real(dp), intent(inout) :: Tij(nmo,nmo,npairs)

    ! Everything else
    integer(is)             :: ipair,Bindx,Kindx
    integer(is)             :: ioff,iopen,i,i1,ikcsf,ibcsf
    real(dp)                :: bcoe,kcoe,prod
    
    ! Loop over TDM pairs
    do ipair=1,npairs

       ! Bra and ket eigenvector indices
       Bindx=Bmap(ipair)
       Kindx=Kmap(ipair)

       ! Initialise the spincp offset
       ioff=offspincp(6+nopen)

       ! Loop over open shells
       do iopen=1,nopen

          ! MO index
          i=socc(iopen)
          i1=m2c(i)

          ! Loop over CSF pairs
          do ikcsf=kcsfs(ikconf),kcsfs(ikconf+1)-1
             kcoe=vecK(ikcsf,Kindx)
                
             ! Loop over bra CSFs
             do ibcsf=bcsfs(ibconf),bcsfs(ibconf+1)-1
                bcoe=vecB(ibcsf,Bindx)

                ! Contribution to the 1-TDM
                prod=kcoe*bcoe*spincp(ioff)
                Tij(i1,i1,ipair)=Tij(i1,i1,ipair)+prod
                      
                ! Increment the spincp offset
                ioff=ioff+1
                
             enddo

          enddo
             
       enddo
          
    enddo
       
    return
    
  end subroutine ondiag_batch
  
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
       knsp=cfgK%ncsfs(knopen)
       
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
          bnsp=cfgB%ncsfs(bnopen)

          ! Get the indices of the MOs involved in the excitation
          hlist=0
          plist=0
          call get_exci_indices(cfgK%conf0h(:,:,ikconf),&
               cfgB%conf0h(:,:,ibconf),n_int_I,hlist(1:nexci),&
               plist(1:nexci),nexci)

          ! Get the spin-coupling coefficients
          scc(1:knsp*bnsp)=triplet_scc(kval,knsp,bnsp,ksop_full,&
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
                   prod=kcoe*bcoe*scc(counter)
                   Tij(a,i,ipair)=Tij(a,i,ipair)+prod
                   
                enddo
                
             enddo
             
          enddo
             
       enddo
          
    enddo
    
    return
    
  end subroutine Tij_0h_0h

!######################################################################
! Tij_0h_1h: Calculation of the Ref-1I and Ref-1E contributions to the
!            MRCI triplet TDMs
!######################################################################
  subroutine Tij_0h_1h(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,&
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

    ! Working arrays
    integer(ib)             :: kconf_full(n_int,2)
    integer(ib)             :: ksop_full(n_int,2)
    
    ! Everything else
    integer(is)             :: ikconf,n,nac,nexci,n_int_I
    integer(is)             :: knsp,knopen
    logical                 :: transpose
    
!----------------------------------------------------------------------
! Contributions of the ket reference and bra 1I & 1E CSFs
!----------------------------------------------------------------------
    n_int_I=cfgK%n_int_I

    ! Work with elements < Bra | T_ij^(1,k) | Ket >
    transpose=.false.
    
    ! Loop over ket reference configurations
    do ikconf=1,cfgK%n0h

       ! Number of open shells in the ket configuration
       knopen=sop_nopen(cfgK%sop0h(:,:,ikconf),n_int_I)
       
       ! Number of ket CSFs
       knsp=cfgK%ncsfs(knopen)

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
             call triplet_tdm_batch(&
                  kval,n,ikconf,kconf_full,ksop_full,&
                  knopen,knsp,nbefore,&
                  cfgB%ncsfs,&              ! bra ncsfs array
                  cfgB%n1I,cfgK%n0h,&       ! no. bra and ket confs
                  cfgB%conf1I,cfgB%sop1I,&  ! bra confs and SOPs
                  cfgB%n1h,cfgB%off1I,&     ! no. bra hole confs and offsets
                  cfgB%csfs1I,cfgK%csfs0h,& ! bra and ket CSF offsets
                  csfdimB,csfdimK,nvecB,nvecK,&
                  vecB,vecK,npairs,Tij,Bmap,Kmap,&
                  cfgB%m2c,transpose)
          endif

          ! Ket Ref - bra 1E contributions
          if (nac <= 1 .and. cfgB%n1E > 0 &
               .and. cfgB%off1E(n) /= cfgB%off1E(n+1)) then
             call triplet_tdm_batch(&
                  kval,n,ikconf,kconf_full,ksop_full,&
                  knopen,knsp,nbefore,&
                  cfgB%ncsfs,&              ! bra ncsfs array
                  cfgB%n1E,cfgK%n0h,&       ! no. bra and ket confs
                  cfgB%conf1E,cfgB%sop1E,&  ! bra confs and SOPs
                  cfgB%n1h,cfgB%off1E,&     ! no. bra hole confs and offsets
                  cfgB%csfs1E,cfgK%csfs0h,& ! bra and ket CSF offsets
                  csfdimB,csfdimK,nvecB,nvecK,&
                  vecB,vecK,npairs,Tij,Bmap,Kmap,&
                  cfgB%m2c,transpose)
          endif
          
       enddo
          
    enddo

!----------------------------------------------------------------------
! Contributions of the bra reference and ket 1I & 1E CSFs
!----------------------------------------------------------------------
    ! Work with elements < Ket | T_ji^(1,-k) | Bra >
    transpose=.true.

    ! Loop over ket reference configurations
    do ikconf=1,cfgB%n0h

       ! Number of open shells in the ket configuration
       knopen=sop_nopen(cfgB%sop0h(:,:,ikconf),n_int_I)
       
       ! Number of ket CSFs
       knsp=cfgB%ncsfs(knopen)

       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_I,:)=cfgB%conf0h(:,:,ikconf)
       ksop_full(1:n_int_I,:)=cfgB%sop0h(:,:,ikconf)
       
       ! Get the number of open shells preceding each ket conf MO
       call nobefore(ksop_full,nbefore)

       ! Loop over bra 1-hole configurations
       do n=1,cfgK%n1h

          ! Number of creation and annihilation operators linking the
          ! reference and 1-hole configurations
          nac=n_create_annihilate(cfgB%conf0h(1:n_int_I,:,ikconf), &
               cfgK%conf1h(1:n_int_I,:,n),n_int_I)

          ! Ket Ref - bra 1I contributions
          if (nac <= 3 .and. cfgK%n1I > 0 &
               .and. cfgK%off1I(n) /= cfgK%off1I(n+1)) then
             call triplet_tdm_batch(&
                  kval,n,ikconf,kconf_full,ksop_full,&
                  knopen,knsp,nbefore,&
                  cfgK%ncsfs,&              ! bra ncsfs array
                  cfgK%n1I,cfgB%n0h,&       ! no. bra and ket confs
                  cfgK%conf1I,cfgK%sop1I,&  ! bra confs and SOPs
                  cfgK%n1h,cfgK%off1I,&     ! no. bra hole confs and offsets
                  cfgK%csfs1I,cfgB%csfs0h,& ! bra and ket CSF offsets
                  csfdimK,csfdimB,nvecK,nvecB,&
                  vecK,vecB,npairs,Tij,Kmap,Bmap,&
                  cfgK%m2c,transpose)
          endif
          
          ! Ket Ref - bra 1E contributions
          if (nac <= 1 .and. cfgK%n1E > 0 &
               .and. cfgK%off1E(n) /= cfgK%off1E(n+1)) then
             call triplet_tdm_batch(&
                  kval,n,ikconf,kconf_full,&
                  ksop_full,knopen,knsp,nbefore,&
                  cfgB%ncsfs,&              ! bra ncsfs array
                  cfgK%n1E,cfgB%n0h,&       ! no. bra and ket confs
                  cfgK%conf1E,cfgK%sop1E,&  ! bra confs and SOPs
                  cfgK%n1h,cfgK%off1E,&     ! no. bra hole confs and offsets
                  cfgK%csfs1E,cfgB%csfs0h,& ! bra and ket CSF offsets
                  csfdimK,csfdimB,nvecK,nvecB,&
                  vecK,vecB,npairs,Tij,Kmap,Bmap,&
                  cfgK%m2c,transpose)
          endif
          
       enddo
       
    enddo
       
    return
    
  end subroutine Tij_0h_1h

!######################################################################
! Tij_1h_1h: Calculation of the 1I-1I, 1I-1E and 1E-1E contributions
!            to the MRCI triplet TDMs
!######################################################################
  subroutine Tij_1h_1h(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,&
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

    ! Triplet TDMs
    real(dp), intent(inout) :: Tij(nmo,nmo,npairs)

    ! Bra-ket pair to eigenvector mapping arrays
    integer(is), intent(in) :: Bmap(npairs),Kmap(npairs)

    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)

    ! Working arrays
    integer(ib), allocatable :: kconf_int(:,:)
    integer(ib), allocatable :: ksop_int(:,:)
    integer(ib)              :: kconf_full(n_int,2)
    integer(ib)              :: ksop_full(n_int,2)
    integer(is), parameter   :: maxexci=1
    integer(is)              :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: ikconf,kn,bn,nac,nac1,nexci,n_int_I
    integer(is)             :: knsp,knopen
    logical                 :: transpose

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfgK%n_int_I

    allocate(kconf_int(n_int_I,2))
    allocate(ksop_int(n_int_I,2))
    kconf_int=0_ib
    ksop_int=0_ib

!----------------------------------------------------------------------
! Compute the 1H-1H contributions to the triplet TDMs
!----------------------------------------------------------------------
    ! Loop over ket 1-hole configurations
    do kn=1,cfgK%n1h
    
       ! Loop over bra 1-hole configurations
       do bn=1,cfgB%n1h

          ! Number of creation and annihilation operators linking the
          ! bra and ket 1-hole configurations
          nac1=n_create_annihilate(cfgK%conf1h(1:n_int_I,:,kn),&
               cfgB%conf1h(1:n_int_I,:,bn),n_int_I)

          ! Cycle if the full configurations generated by the 1-hole
          ! configurations cannot be interacting wrt the singlet
          ! excitation operators E_a^i
          if (nac1 > 4) cycle

          !
          ! Ket 1I, bra 1I and 1E contributions
          !
          ! Work with elements < Bra | T_ij^(1,k) | Ket >
          transpose=.false.
          if (cfgK%n1I > 0) then

             ! Loop over ket 1I configurations
             do ikconf=cfgK%off1I(kn),cfgK%off1I(kn+1)-1

                ! Ket 1I configuration
                kconf_int=cfgK%conf1I(1:n_int_I,:,ikconf)
                ksop_int=cfgK%sop1I(1:n_int_I,:,ikconf)
                kconf_full=0_ib
                ksop_full=0_ib
                kconf_full(1:n_int_I,:)=kconf_int
                ksop_full(1:n_int_I,:)=ksop_int

                ! Number of creation and annihilation operators linking
                ! the ket 1I and bra 1-hole configurations
                nac=n_create_annihilate(kconf_int,&
                     cfgB%conf1h(1:n_int_I,:,bn),n_int_I)

                ! Cycle if the the bra 1-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 1I configuration wrt the triplet excitation
                ! operators T_ij^(1,k)
                if (nac > 3) cycle

                ! Number of open shells in the ket 1I configuration
                knopen=sop_nopen(ksop_int(1:n_int_I,:),n_int_I)
                
                ! Number of ket 1I CSFs
                knsp=cfgK%ncsfs(knopen)

                ! Get the number of open shells preceding each ket
                ! conf MO
                call nobefore(ksop_full,nbefore)

                ! Ket 1I - bra 1I contributions
                if (nac <= 3 .and. &
                     cfgB%off1I(bn) /= cfgB%off1I(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,&
                        ksop_full,knopen,knsp,nbefore,&
                        cfgB%ncsfs,&              ! bra ncsfs array
                        cfgB%n1I,cfgK%n1I,&       ! no. bra and ket confs
                        cfgB%conf1I,cfgB%sop1I,&  ! bra confs and SOPs
                        cfgB%n1h,cfgB%off1I,&     ! no. bra hole confs and offsets
                        cfgB%csfs1I,cfgK%csfs1I,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif

                ! Ket 1I - bra 1E contributions
                if (nac <= 1 &
                     .and. cfgB%off1E(bn) /= cfgB%off1E(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,&
                        ksop_full,knopen,knsp,nbefore,&
                        cfgB%ncsfs,&              ! bra ncsfs array
                        cfgB%n1E,cfgK%n1I,&       ! no. bra and ket confs
                        cfgB%conf1E,cfgB%sop1E,&  ! bra confs and SOPs
                        cfgB%n1h,cfgB%off1E,&     ! no. bra hole confs and offsets
                        cfgB%csfs1E,cfgK%csfs1I,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif
                
             enddo
                
          endif

          !
          ! Ket 1E - bra 1I contributions
          !
          ! Work with elements < Ket | T_ji^(1,-k) | Bra >
          transpose=.true.
          if (cfgB%n1I > 0) then

             ! Loop over ket 1I configurations
             do ikconf=cfgB%off1I(bn),cfgB%off1I(bn+1)-1

                ! Ket 1I configuration
                kconf_int=cfgB%conf1I(1:n_int_I,:,ikconf)
                ksop_int=cfgB%sop1I(1:n_int_I,:,ikconf)
                kconf_full=0_ib
                ksop_full=0_ib
                kconf_full(1:n_int_I,:)=kconf_int
                ksop_full(1:n_int_I,:)=ksop_int

                ! Number of creation and annihilation operators linking
                ! the ket 1I and bra 1-hole configurations
                nac=n_create_annihilate(kconf_int,&
                     cfgK%conf1h(1:n_int_I,:,kn),n_int_I)

                ! Cycle if the the bra 1-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 1I configuration wrt the singlet excitation
                ! operators E_a^i
                if (nac > 3) cycle

                ! Number of open shells in the ket 1I configuration
                knopen=sop_nopen(ksop_int(1:n_int_I,:),n_int_I)
                
                ! Number of ket 1I CSFs
                knsp=cfgB%ncsfs(knopen)

                ! Get the number of open shells preceding each ket
                ! conf MO
                call nobefore(ksop_full,nbefore)

                ! Ket 1I - bra 1E contributions
                if (nac <= 1 &
                     .and. cfgK%off1E(kn) /= cfgK%off1E(kn+1)) then
                   call triplet_tdm_batch(&
                        kval,kn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgK%ncsfs,&              ! bra ncsfs array
                        cfgK%n1E,cfgB%n1I,&       ! no. bra and ket confs
                        cfgK%conf1E,cfgK%sop1E,&  ! bra confs and SOPs
                        cfgK%n1h,cfgK%off1E,&     ! no. bra hole confs and offsets
                        cfgK%csfs1E,cfgB%csfs1I,& ! bra and ket CSF offsets
                        csfdimK,csfdimB,nvecK,nvecB,&
                        vecK,vecB,npairs,Tij,Kmap,Bmap,&
                        cfgK%m2c,transpose)
                endif
                
             enddo
                
          endif

          !
          ! Ket 1E, bra 1E contributions
          !
          ! Work with elements < Bra | T_ij^(1,k) | Ket >
          transpose=.false.
          
          ! Cycle if the bra 1-hole configuration doesn't generate
          ! any 1E configurations
          if (cfgB%off1E(bn) == cfgB%off1E(bn+1)) cycle
          
          ! Cycle if the the bra and ket 1-hole configurations
          ! cannot generate interacting 1E configurations wrt the
          ! singlet excitation operators E_a^i
          if (nac1 > 2) cycle

          ! Loop over ket 1E configurations
          do ikconf=cfgK%off1E(kn),cfgK%off1E(kn+1)-1

             ! Ket 1E configuration
             kconf_full=cfgK%conf1E(:,:,ikconf)
             ksop_full=cfgK%sop1E(:,:,ikconf)
             
             ! Number of open shells in the ket 1E configuration
             knopen=sop_nopen(ksop_full,n_int)

             ! Number of ket 1E CSFs
             knsp=cfgK%ncsfs(knopen)

             ! Get the number of open shells preceding each ket
             ! conf MO
             call nobefore(ksop_full,nbefore)

             ! Ket 1E - bra 1E matrix contributions
             call triplet_tdm_batch(&
                  kval,bn,ikconf,kconf_full,ksop_full,&
                  knopen,knsp,nbefore,&
                  cfgB%ncsfs,&              ! bra ncsfs array
                  cfgB%n1E,cfgK%n1E,&       ! no. bra and ket confs
                  cfgB%conf1E,cfgB%sop1E,&  ! bra confs and SOPs
                  cfgB%n1h,cfgB%off1E,&     ! no. bra hole confs and offsets
                  cfgB%csfs1E,cfgK%csfs1E,& ! bra and ket CSF offsets
                  csfdimB,csfdimK,nvecB,nvecK,&
                  vecB,vecK,npairs,Tij,Bmap,Kmap,&
                  cfgB%m2c,transpose)
             
          enddo
          
       enddo
       
    enddo
    
    return
    
  end subroutine Tij_1h_1h

!######################################################################
! Tij_2h_1h: Calculation of the 2I-1I, 2I-1E, 2E-1I, 2E-1E, 1I1E-1I,
!            and 1I1E-1E contributions to the MRCI triplet TDMs
!######################################################################
  subroutine Tij_2h_1h(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,&
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

    ! Triplet TDMs
    real(dp), intent(inout) :: Tij(nmo,nmo,npairs)

    ! Bra-ket pair to eigenvector mapping arrays
    integer(is), intent(in) :: Bmap(npairs),Kmap(npairs)

    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)

    ! Working arrays
    integer(ib), allocatable :: kconf_int(:,:)
    integer(ib), allocatable :: ksop_int(:,:)
    integer(ib)              :: kconf_full(n_int,2)
    integer(ib)              :: ksop_full(n_int,2)
    integer(ib)              :: bconf1h_full(n_int,2)
    integer(is), parameter   :: maxexci=1
    integer(is)              :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: ikconf,kn,bn,nac,nac1,nexci,n_int_I
    integer(is)             :: knsp,knopen
    logical                 :: transpose

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfgK%n_int_I

    allocate(kconf_int(n_int_I,2))
    allocate(ksop_int(n_int_I,2))
    kconf_int=0_ib
    ksop_int=0_ib

!----------------------------------------------------------------------
! Calculate the ket 2H - bra 1H contributions to the triplet TDMs
!----------------------------------------------------------------------
    ! Work with elements < Bra | T_ij^(1,k) | Ket >
    transpose=.false.

    ! Loop over ket 2-hole configurations
    do kn=1,cfgK%n2h

       ! Loop over bra 1-hole configurations
       do bn=1,cfgB%n1h

          ! Bra 1-hole configuration in the full MO space
          bconf1h_full=0_ib
          bconf1h_full(1:n_int_I,:)=cfgB%conf1h(1:n_int_I,:,bn)

          ! Number of creation and annihilation operators linking the
          ! bra and ket 1- and 2-hole configurations
          nac1=n_create_annihilate(cfgK%conf2h(1:n_int_I,:,kn),&
               cfgB%conf1h(1:n_int_I,:,bn),n_int_I)

          ! Cycle if the full configurations generated by the hole
          ! configurations cannot be interacting
          if (nac1 > 5) cycle

          !
          ! Ket: 2I
          ! Bra: 1I and 1E
          !
          if (cfgK%n2I > 0) then

             ! Loop over ket 2I configurations
             do ikconf=cfgK%off2I(kn),cfgK%off2I(kn+1)-1

                ! Ket 2I configuration
                kconf_int=cfgK%conf2I(1:n_int_I,:,ikconf)
                ksop_int=cfgK%sop2I(1:n_int_I,:,ikconf)
                kconf_full=0_ib
                ksop_full=0_ib
                kconf_full(1:n_int_I,:)=kconf_int
                ksop_full(1:n_int_I,:)=ksop_int

                ! Number of creation and annihilation operators linking
                ! the ket 2I and bra 1-hole configurations
                nac=n_create_annihilate(kconf_int,&
                     cfgB%conf1h(1:n_int_I,:,bn),n_int_I)
                
                ! Cycle if the the bra 1-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 2I configuration
                if (nac > 3) cycle

                ! Number of open shells in the ket 2I configuration
                knopen=sop_nopen(ksop_int(1:n_int_I,:),n_int_I)
                
                ! Number of ket 2I CSFs
                knsp=cfgK%ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! Ket 2I - bra 1I contributions
                if (nac <= 3 &
                     .and. cfgB%off1I(bn) /= cfgB%off1I(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&              ! bra ncsfs array
                        cfgB%n1I,cfgK%n2I,&       ! no. bra and ket confs
                        cfgB%conf1I,cfgB%sop1I,&  ! bra confs and SOPs
                        cfgB%n1h,cfgB%off1I,&     ! no. bra hole confs and offsets
                        cfgB%csfs1I,cfgK%csfs2I,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif

                ! Ket 2I - bra 1E contributions
                if (nac <= 1 &
                     .and. cfgB%off1E(bn) /= cfgB%off1E(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&              ! bra ncsfs array
                        cfgB%n1E,cfgK%n2I,&       ! no. bra and ket confs
                        cfgB%conf1E,cfgB%sop1E,&  ! bra confs and SOPs
                        cfgB%n1h,cfgB%off1E,&     ! no. bra hole confs and offsets
                        cfgB%csfs1E,cfgK%csfs2I,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif
                   
             enddo

          endif

          !
          ! Ket: 2E
          ! Bra: 1I and 1E
          !
          if (cfgK%n2E > 0 .and. nac1 <= 1) then
             
             ! Loop over ket 2E configurations
             do ikconf=cfgK%off2E(kn),cfgK%off2E(kn+1)-1

                ! Ket 2E configuration
                kconf_full=cfgK%conf2E(:,:,ikconf)
                ksop_full=cfgK%sop2E(:,:,ikconf)

                ! Number of open shells in the ket 2E configuration
                knopen=sop_nopen(ksop_full,n_int)
                
                ! Number of ket 2E CSFs
                knsp=cfgK%ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! Ket 2E - bra 1I contributions
                if (cfgB%n1I > 0 &
                     .and. cfgB%off1I(bn) /= cfgB%off1I(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&              ! bra ncsfs array
                        cfgB%n1I,cfgK%n2E,&       ! no. bra and ket confs
                        cfgB%conf1I,cfgB%sop1I,&  ! bra confs and SOPs
                        cfgB%n1h,cfgB%off1I,&     ! no. bra hole confs and offsets
                        cfgB%csfs1I,cfgK%csfs2E,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif

                ! Ket 2E - bra 1E contributions
                if (cfgB%off1E(bn) /= cfgB%off1E(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&              ! bra ncsfs array
                        cfgB%n1E,cfgK%n2E,&       ! no. bra and ket confs
                        cfgB%conf1E,cfgB%sop1E,&  ! bra confs and SOPs
                        cfgB%n1h,cfgB%off1E,&     ! no. bra hole confs and offsets
                        cfgB%csfs1E,cfgK%csfs2E,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif
                   
             enddo

          endif

          !
          ! Ket: 1I1E
          ! Bra: 1I and 1E
          !
          if (cfgK%n1I1E > 0) then
             
             ! Loop over ket 1I1E configurations
             do ikconf=cfgK%off1I1E(kn),cfgK%off1I1E(kn+1)-1

                ! Ket 1I1E configuration
                kconf_full=cfgK%conf1I1E(:,:,ikconf)
                ksop_full=cfgK%sop1I1E(:,:,ikconf)

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
                knsp=cfgK%ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! Ket 1I1E - bra 1I contributions
                if (cfgB%off1I(bn) /= cfgB%off1I(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&                ! bra ncsfs array
                        cfgB%n1I,cfgK%n1I1E,&       ! no. bra and ket confs
                        cfgB%conf1I,cfgB%sop1I,&    ! bra confs and SOPs
                        cfgB%n1h,cfgB%off1I,&       ! no. bra hole confs and offsets
                        cfgB%csfs1I,cfgK%csfs1I1E,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif

                ! Ket 1I1E - bra 1E contributions
                if (cfgB%off1E(bn) /= cfgB%off1E(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&                ! bra ncsfs array
                        cfgB%n1E,cfgK%n1I1E,&       ! no. bra and ket confs
                        cfgB%conf1E,cfgB%sop1E,&    ! bra confs and SOPs
                        cfgB%n1h,cfgB%off1E,&       ! no. bra hole confs and offsets
                        cfgB%csfs1E,cfgK%csfs1I1E,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif
                
             enddo

          endif
                
       enddo
       
    enddo

!----------------------------------------------------------------------
! Calculate the bra 2H - ket 1H contributions to the triplet TDMs
!----------------------------------------------------------------------
    ! Work with elements < Ket | T_ji^(1,-k) | Bra >
    transpose=.true.

    ! Loop over ket 2-hole configurations
    do kn=1,cfgB%n2h
       
       ! Loop over bra 1-hole configurations
       do bn=1,cfgK%n1h

          ! Bra 1-hole configuration in the full MO space
          bconf1h_full=0_ib
          bconf1h_full(1:n_int_I,:)=cfgK%conf1h(1:n_int_I,:,bn)

          ! Number of creation and annihilation operators linking the
          ! bra and ket 1- and 2-hole configurations
          nac1=n_create_annihilate(cfgB%conf2h(1:n_int_I,:,kn),&
               cfgK%conf1h(1:n_int_I,:,bn),n_int_I)

          ! Cycle if the full configurations generated by the hole
          ! configurations cannot be interacting
          if (nac1 > 5) cycle

          !
          ! Ket: 2I
          ! Bra: 1I and 1E
          !
          if (cfgB%n2I > 0) then

             ! Loop over ket 2I configurations
             do ikconf=cfgB%off2I(kn),cfgB%off2I(kn+1)-1

                ! Ket 2I configuration
                kconf_int=cfgB%conf2I(1:n_int_I,:,ikconf)
                ksop_int=cfgB%sop2I(1:n_int_I,:,ikconf)
                kconf_full=0_ib
                ksop_full=0_ib
                kconf_full(1:n_int_I,:)=kconf_int
                ksop_full(1:n_int_I,:)=ksop_int

                ! Number of creation and annihilation operators linking
                ! the ket 2I and bra 1-hole configurations
                nac=n_create_annihilate(kconf_int,&
                     cfgK%conf1h(1:n_int_I,:,bn),n_int_I)
                
                ! Cycle if the the bra 1-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 2I configuration
                if (nac > 3) cycle

                ! Number of open shells in the ket 2I configuration
                knopen=sop_nopen(ksop_int(1:n_int_I,:),n_int_I)
                
                ! Number of ket 2I CSFs
                knsp=cfgB%ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! Ket 2I - bra 1I contributions
                if (nac <= 3 &
                     .and. cfgK%off1I(bn) /= cfgK%off1I(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgK%ncsfs,&              ! bra ncsfs array
                        cfgK%n1I,cfgB%n2I,&       ! no. bra and ket confs
                        cfgK%conf1I,cfgK%sop1I,&  ! bra confs and SOPs
                        cfgK%n1h,cfgK%off1I,&     ! no. bra hole confs and offsets
                        cfgK%csfs1I,cfgB%csfs2I,& ! bra and ket CSF offsets
                        csfdimK,csfdimB,nvecK,nvecB,&
                        vecK,vecB,npairs,Tij,Kmap,Bmap,&
                        cfgK%m2c,transpose)
                endif

                ! Ket 2I - bra 1E contributions
                if (nac <= 1 &
                     .and. cfgK%off1E(bn) /= cfgK%off1E(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgK%ncsfs,&              ! bra ncsfs array
                        cfgK%n1E,cfgB%n2I,&       ! no. bra and ket confs
                        cfgK%conf1E,cfgK%sop1E,&  ! bra confs and SOPs
                        cfgK%n1h,cfgK%off1E,&     ! no. bra hole confs and offsets
                        cfgK%csfs1E,cfgB%csfs2I,& ! bra and ket CSF offsets
                        csfdimK,csfdimB,nvecK,nvecB,&
                        vecK,vecB,npairs,Tij,Kmap,Bmap,&
                        cfgK%m2c,transpose)
                endif
                   
             enddo
                
          endif

          !
          ! Ket: 2E
          ! Bra: 1I and 1E
          !
          if (cfgB%n2E > 0 .and. nac1 <= 1) then
             
             ! Loop over ket 2E configurations
             do ikconf=cfgB%off2E(kn),cfgB%off2E(kn+1)-1

                ! Ket 2E configuration
                kconf_full=cfgB%conf2E(:,:,ikconf)
                ksop_full=cfgB%sop2E(:,:,ikconf)

                ! Number of open shells in the ket 2E configuration
                knopen=sop_nopen(ksop_full,n_int)
                
                ! Number of ket 2E CSFs
                knsp=cfgB%ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! Ket 2E - bra 1I contributions
                if (cfgK%n1I > 0 &
                     .and. cfgK%off1I(bn) /= cfgK%off1I(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgK%ncsfs,&              ! bra ncsfs array
                        cfgK%n1I,cfgB%n2E,&       ! no. bra and ket confs
                        cfgK%conf1I,cfgK%sop1I,&  ! bra confs and SOPs
                        cfgK%n1h,cfgK%off1I,&     ! no. bra hole confs and offsets
                        cfgK%csfs1I,cfgB%csfs2E,& ! bra and ket CSF offsets
                        csfdimK,csfdimB,nvecK,nvecB,&
                        vecK,vecB,npairs,Tij,Kmap,Bmap,&
                        cfgK%m2c,transpose)
                endif

                ! Ket 2E - bra 1E contributions
                if (cfgK%off1E(bn) /= cfgK%off1E(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgK%ncsfs,&              ! bra ncsfs array
                        cfgK%n1E,cfgB%n2E,&       ! no. bra and ket confs
                        cfgK%conf1E,cfgK%sop1E,&  ! bra confs and SOPs
                        cfgK%n1h,cfgK%off1E,&     ! no. bra hole confs and offsets
                        cfgK%csfs1E,cfgB%csfs2E,& ! bra and ket CSF offsets
                        csfdimK,csfdimB,nvecK,nvecB,&
                        vecK,vecB,npairs,Tij,Kmap,Bmap,&
                        cfgK%m2c,transpose)
                endif
                   
             enddo

          endif

          !
          ! Ket: 1I1E
          ! Bra: 1I and 1E
          !
          if (cfgB%n1I1E > 0) then
             
             ! Loop over ket 1I1E configurations
             do ikconf=cfgB%off1I1E(kn),cfgB%off1I1E(kn+1)-1

                ! Ket 1I1E configuration
                kconf_full=cfgB%conf1I1E(:,:,ikconf)
                ksop_full=cfgB%sop1I1E(:,:,ikconf)

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
                knsp=cfgB%ncsfs(knopen)

                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)

                ! Ket 1I1E - bra 1I contributions
                if (cfgK%off1I(bn) /= cfgK%off1I(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgK%ncsfs,&                ! bra ncsfs array
                        cfgK%n1I,cfgB%n1I1E,&       ! no. bra and ket confs
                        cfgK%conf1I,cfgK%sop1I,&    ! bra confs and SOPs
                        cfgK%n1h,cfgK%off1I,&       ! no. bra hole confs and offsets
                        cfgK%csfs1I,cfgB%csfs1I1E,& ! bra and ket CSF offsets
                        csfdimK,csfdimB,nvecK,nvecB,&
                        vecK,vecB,npairs,Tij,Kmap,Bmap,&
                        cfgK%m2c,transpose)
                endif

                ! Ket 1I1E - bra 1E contributions
                if (cfgK%off1E(bn) /= cfgK%off1E(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgK%ncsfs,&              ! bra ncsfs array
                        cfgK%n1E,cfgB%n1I1E,&       ! no. bra and ket confs
                        cfgK%conf1E,cfgK%sop1E,&    ! bra confs and SOPs
                        cfgK%n1h,cfgK%off1E,&       ! no. bra hole confs and offsets
                        cfgK%csfs1E,cfgB%csfs1I1E,& ! bra and ket CSF offsets
                        csfdimK,csfdimB,nvecK,nvecB,&
                        vecK,vecB,npairs,Tij,Kmap,Bmap,&
                        cfgK%m2c,transpose)
                endif
                   
             enddo

          endif
             
       enddo

    enddo
    
    return
    
  end subroutine Tij_2h_1h

!######################################################################
! Tij_2h_2h: Calculation of the 2I-2I, 2I-2E, 2I-1I1E, 2E-2E, 2I-1I1E,
!            and 1I1E-1I1E contributions to the MRCI triplet TDMs
!######################################################################
  subroutine Tij_2h_2h(kval,cfgB,cfgK,csfdimB,csfdimK,nvecB,nvecK,&
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

    ! Working arrays
    integer(ib), allocatable :: kconf_int(:,:)
    integer(ib), allocatable :: ksop_int(:,:)
    integer(ib)              :: kconf_full(n_int,2)
    integer(ib)              :: ksop_full(n_int,2)
    integer(ib)              :: bconf2h_full(n_int,2)
    integer(ib)              :: kconf2h_full(n_int,2)
    integer(is), parameter   :: maxexci=1
    integer(is)              :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: ikconf,kn,bn,nac,nac1,nexci,n_int_I
    integer(is)             :: knsp,knopen
    logical                 :: transpose

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfgK%n_int_I

    allocate(kconf_int(n_int_I,2))
    allocate(ksop_int(n_int_I,2))
    kconf_int=0_ib
    ksop_int=0_ib

!----------------------------------------------------------------------
! 2H-2H contributions to the 1-RDMs
!----------------------------------------------------------------------
    ! Loop over ket 2-hole configurations
    do kn=1,cfgK%n2h

       ! Ket 2-hole configuration in the full MO space
       kconf2h_full=0_ib
       kconf2h_full(1:n_int_I,:)=cfgK%conf2h(1:n_int_I,:,kn)
       
       ! Loop over bra 2-hole configurations
       do bn=1,cfgB%n2h

          ! Bra 2-hole configuration in the full MO space
          bconf2h_full=0_ib
          bconf2h_full(1:n_int_I,:)=cfgB%conf2h(1:n_int_I,:,bn)

          ! Number of creation and annihilation operators linking the
          ! bra and ket 2-hole configurations
          nac1=n_create_annihilate(cfgK%conf2h(1:n_int_I,:,kn),&
               cfgB%conf2h(1:n_int_I,:,bn),n_int_I)

          ! Cycle if the full configurations generated by the hole
          ! configurations cannot be interacting
          if (nac1 > 6) cycle

          !
          ! Ket: 2I
          ! Bra: 2I, 2E and 1I1E
          !
          ! Work with elements < Bra | T_ij^(1,k) | Ket >
          transpose=.false.
          if (cfgK%n2I > 0) then
             
             ! Loop over ket 2I configurations
             do ikconf=cfgK%off2I(kn),cfgK%off2I(kn+1)-1
          
                ! Ket 2I configuration
                kconf_int=cfgK%conf2I(1:n_int_I,:,ikconf)
                ksop_int=cfgK%sop2I(1:n_int_I,:,ikconf)
                kconf_full=0_ib
                ksop_full=0_ib
                kconf_full(1:n_int_I,:)=kconf_int
                ksop_full(1:n_int_I,:)=ksop_int
          
                ! Number of creation and annihilation operators linking
                ! the ket 2I and bra 2-hole configurations
                nac=n_create_annihilate(kconf_int,&
                     cfgB%conf2h(1:n_int_I,:,bn),n_int_I)
          
                ! Cycle if the the bra 2-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 2I configuration
                if (nac > 4) cycle
          
                ! Number of open shells in the ket 2I configuration
                knopen=sop_nopen(ksop_int(1:n_int_I,:),n_int_I)
                
                ! Number of ket 2I CSFs
                knsp=cfgK%ncsfs(knopen)
          
                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)
          
                ! Ket 2I - bra 2I contributions
                if (nac <= 4 .and. cfgB%n2I > 0 .and. &
                     cfgB%off2I(bn) /= cfgB%off2I(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&              ! bra ncsfs array
                        cfgB%n2I,cfgK%n2I,&       ! no. bra and ket confs
                        cfgB%conf2I,cfgB%sop2I,&  ! bra confs and SOPs
                        cfgB%n2h,cfgB%off2I,&     ! no. bra hole confs and offsets
                        cfgB%csfs2I,cfgK%csfs2I,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif
                
                ! Ket 2I - bra 2E contributions
                if (nac == 0 .and. cfgB%n2E > 0 &
                     .and. cfgB%off2E(bn) /= cfgB%off2E(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&              ! bra ncsfs array
                        cfgB%n2E,cfgK%n2I,&       ! no. bra and ket confs
                        cfgB%conf2E,cfgB%sop2E,&  ! bra confs and SOPs
                        cfgB%n2h,cfgB%off2E,&     ! no. bra hole confs and offsets
                        cfgB%csfs2E,cfgK%csfs2I,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif
          
                ! Ket 2I - bra 1I1E contributions
                if (nac <= 2 &
                     .and. cfgB%n1I1E > 0 &
                     .and. cfgB%off1I1E(bn) /= cfgB%off1I1E(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&                 ! bra ncsfs array
                        cfgB%n1I1E,cfgK%n2I,&        ! no. bra and ket confs
                        cfgB%conf1I1E,cfgB%sop1I1E,& ! bra confs and SOPs
                        cfgB%n2h,cfgB%off1I1E,&      ! no. bra hole confs and offsets
                        cfgB%csfs1I1E,cfgK%csfs2I,&  ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif
                   
             enddo
          
          endif

          !
          ! Bra: 2I
          ! Ket: 2E and 1I1E
          !
          ! Work with elements < Ket | T_ji^(1,-k) | Bra >
          transpose=.true.
          if (cfgB%n2I > 0) then
             
             ! Loop over ket 2I configurations
             do ikconf=cfgB%off2I(bn),cfgB%off2I(bn+1)-1
          
                ! Ket 2I configuration
                kconf_int=cfgB%conf2I(1:n_int_I,:,ikconf)
                ksop_int=cfgB%sop2I(1:n_int_I,:,ikconf)
                kconf_full=0_ib
                ksop_full=0_ib
                kconf_full(1:n_int_I,:)=kconf_int
                ksop_full(1:n_int_I,:)=ksop_int
          
                ! Number of creation and annihilation operators linking
                ! the ket 2I and bra 2-hole configurations
                nac=n_create_annihilate(kconf_int,&
                     cfgK%conf2h(1:n_int_I,:,kn),n_int_I)
          
                ! Cycle if the the bra 2-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 2I configuration
                if (nac > 4) cycle
                
                ! Number of open shells in the ket 2I configuration
                knopen=sop_nopen(ksop_int(1:n_int_I,:),n_int_I)
                
                ! Number of ket 2I CSFs
                knsp=cfgB%ncsfs(knopen)
          
                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)
          
                ! Ket 2I - bra 2E contributions
                if (nac == 0 &
                     .and. cfgK%n2E > 0 &
                     .and. cfgK%off2E(kn) /= cfgK%off2E(kn+1)) then
                   call triplet_tdm_batch(&
                        kval,kn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgK%ncsfs,&              ! bra ncsfs array
                        cfgK%n2E,cfgB%n2I,&       ! no. bra and ket confs
                        cfgK%conf2E,cfgK%sop2E,&  ! bra confs and SOPs
                        cfgK%n2h,cfgK%off2E,&     ! no. bra hole confs and offsets
                        cfgK%csfs2E,cfgB%csfs2I,& ! bra and ket CSF offsets
                        csfdimK,csfdimB,nvecK,nvecB,&
                        vecK,vecB,npairs,Tij,Kmap,Bmap,&
                        cfgK%m2c,transpose)
                endif
          
                ! Ket 2I - bra 1I1E contributions
                if (nac <= 2 &
                     .and. cfgK%n1I1E > 0 &
                     .and. cfgK%off1I1E(kn) /= cfgK%off1I1E(kn+1)) then
                   call triplet_tdm_batch(&
                        kval,kn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgK%ncsfs,&                 ! bra ncsfs array
                        cfgK%n1I1E,cfgB%n2I,&        ! no. bra and ket confs
                        cfgK%conf1I1E,cfgK%sop1I1E,& ! bra confs and SOPs
                        cfgK%n2h,cfgK%off1I1E,&      ! no. bra hole confs and offsets
                        cfgK%csfs1I1E,cfgB%csfs2I,&  ! bra and ket CSF offsets
                        csfdimK,csfdimB,nvecK,nvecB,&
                        vecK,vecB,npairs,Tij,Kmap,Bmap,&
                        cfgK%m2c,transpose)
                endif
                   
             enddo
          
          endif

          !
          ! Ket: 2E
          ! Bra: 2E and 1I1E
          !
          ! Work with elements < Bra | T_ij^(1,k) | Ket >
          transpose=.false.
          if (cfgK%n2E > 0 .and. nac1 <= 5) then
          
             ! Loop over ket 2E configurations
             do ikconf=cfgK%off2E(kn),cfgK%off2E(kn+1)-1
          
                ! Ket 2E configuration
                kconf_full=cfgK%conf2E(:,:,ikconf)
                ksop_full=cfgK%sop2E(:,:,ikconf)
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
                knsp=cfgK%ncsfs(knopen)
          
                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)
          
                ! Ket 2E - bra 2E contributions
                if (cfgB%off2E(bn) /= cfgB%off2E(bn+1) .and. &
                     nac1 <= 2 .and. cfgB%n2E /= 0) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&              ! bra ncsfs array
                        cfgB%n2E,cfgK%n2E,&       ! no. bra and ket confs
                        cfgB%conf2E,cfgB%sop2E,&  ! bra confs and SOPs
                        cfgB%n2h,cfgB%off2E,&     ! no. bra hole confs and offsets
                        cfgB%csfs2E,cfgK%csfs2E,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif
          
                ! Ket 2E - bra 1I1E matrix elements
                if (cfgB%off1I1E(bn) /= cfgB%off1I1E(bn+1) &
                     .and. cfgB%n1I1E /= 0) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&                 ! bra ncsfs array
                        cfgB%n1I1E,cfgK%n2E,&        ! no. bra and ket confs
                        cfgB%conf1I1E,cfgB%sop1I1E,& ! bra confs and SOPs
                        cfgB%n2h,cfgB%off1I1E,&      ! no. bra hole confs and offsets
                        cfgB%csfs1I1E,cfgK%csfs2E,&  ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif
                   
             enddo
             
          endif

          !
          ! Bra: 2E
          ! Ket: 1I1E
          !
          ! Work with elements < Ket | T_ji^(1,-k) | Bra >
          transpose=.true.
          if (cfgB%n2E > 0 .and. nac1 <= 5) then
          
             ! Loop over ket 2E configurations
             do ikconf=cfgB%off2E(bn),cfgB%off2E(bn+1)-1
          
                ! Ket 2E configuration
                kconf_full=cfgB%conf2E(:,:,ikconf)
                ksop_full=cfgB%sop2E(:,:,ikconf)
                kconf_int=0_ib
                ksop_int=0_ib
                kconf_int(1:n_int_I,:)=kconf_full(1:n_int_I,:)
                ksop_int(1:n_int_I,:)=ksop_full(1:n_int_I,:)
          
                ! Number of creation and annihilation operators linking
                ! the ket 2E and bra 2-hole configurations
                nac=n_create_annihilate(kconf_full,kconf2h_full,n_int)
          
                ! Cycle if the the bra 2-hole configuration cannot
                ! generate configurations that interact with the
                ! ket 2E configuration
                if (nac > 4) cycle
          
                ! Number of open shells in the ket 2E configuration
                knopen=sop_nopen(ksop_full,n_int)
                
                ! Number of ket 2E CSFs
                knsp=cfgB%ncsfs(knopen)
          
                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)
          
                ! Ket 2E - bra 1I1E matrix elements
                if (cfgK%off1I1E(kn) /= cfgK%off1I1E(kn+1) &
                     .and. cfgK%n1I1E /= 0) then
                   call triplet_tdm_batch(&
                        kval,kn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgK%ncsfs,&                 ! bra ncsfs array
                        cfgK%n1I1E,cfgB%n2E,&        ! no. bra and ket confs
                        cfgK%conf1I1E,cfgK%sop1I1E,& ! bra confs and SOPs
                        cfgK%n2h,cfgK%off1I1E,&      ! no. bra hole confs and offsets
                        cfgK%csfs1I1E,cfgB%csfs2E,&  ! bra and ket CSF offsets
                        csfdimK,csfdimB,nvecK,nvecB,&
                        vecK,vecB,npairs,Tij,Kmap,Bmap,&
                        cfgK%m2c,transpose)
                endif
                
             enddo
          
          endif

          !
          ! Ket: 1I1E
          ! Bra: 1I1E
          !
          ! Work with elements < Bra | T_ij^(1,k) | Ket >
          transpose=.false.
          if (cfgK%n1I1E > 0) then
          
             ! Loop over ket 1I1E configurations
             do ikconf=cfgK%off1I1E(kn),cfgK%off1I1E(kn+1)-1
          
                ! Ket 1I1E configuration
                kconf_full=cfgK%conf1I1E(:,:,ikconf)
                ksop_full=cfgK%sop1I1E(:,:,ikconf)
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
                knsp=cfgK%ncsfs(knopen)
          
                ! Get the number of open shells preceding each ket conf MO
                call nobefore(ksop_full,nbefore)
          
                ! Ket 1I1E - bra 1I1E contributions
                if (cfgB%n1I1E > 0 .and. &
                     cfgB%off1I1E(bn) /= cfgB%off1I1E(bn+1)) then
                   call triplet_tdm_batch(&
                        kval,bn,ikconf,kconf_full,ksop_full,&
                        knopen,knsp,nbefore,&
                        cfgB%ncsfs,&                  ! bra ncsfs array
                        cfgB%n1I1E,cfgK%n1I1E,&       ! no. bra and ket confs
                        cfgB%conf1I1E,cfgB%sop1I1E,&  ! bra confs and SOPs
                        cfgB%n2h,cfgB%off1I1E,&       ! no. bra hole confs and offsets
                        cfgB%csfs1I1E,cfgK%csfs1I1E,& ! bra and ket CSF offsets
                        csfdimB,csfdimK,nvecB,nvecK,&
                        vecB,vecK,npairs,Tij,Bmap,Kmap,&
                        cfgB%m2c,transpose)
                endif
                   
             enddo
          
          endif
          
       enddo

    enddo
       
    return
    
  end subroutine Tij_2h_2h
    
!######################################################################
! triplet_tdm_batch: Computes all the contributions to the triplet
!                    TDMs (component k=kval) from a ket conf and a
!                    single class (1I, 2I, etc.) of confs generated by
!                    a single bra hole conf
!######################################################################
  subroutine triplet_tdm_batch(kval,bn,ikconf,kconf,ksop,knopen,knsp,&
       knbefore,bncsfs,nbconf,nkconf,bconfs,bsops,nh,boffset,bcsfs,&
       kcsfs,csfdimB,csfdimK,nvecB,nvecK,vecB,vecK,npairs,Tij,Bmap,&
       Kmap,m2c,transpose)

    use constants
    use bitglobal
    use mrciutils

    implicit none

    ! Component of the triplet spin tensor operator
    integer(is), intent(in) :: kval

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

    ! No. bra CSFs as a function of the no. open shells
    integer(is), intent(in) :: bncsfs(0:nocase2)
    
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
    
    ! Triplet TDMs
    real(dp), intent(inout) :: Tij(nmo,nmo,npairs)

    ! Bra-ket pair to eigenvector mapping arrays
    integer(is), intent(in) :: Bmap(npairs),Kmap(npairs)
    
    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! Are we work with
    ! - < Ket | T_ji^(1,-k) | Bra >
    ! instead of
    ! + < Bra | T_ij^(1,k) | Ket> ?
    logical, intent(in)     :: transpose
    
    ! Working arrays
    integer(ib)             :: bconf(n_int,2),bsop(n_int,2)
    integer(is), parameter  :: maxexci=1
    integer(is)             :: hlist(maxexci),plist(maxexci)
    
    ! Everything else
    integer(is)             :: ibconf,ipair,ibcsf,ikcsf
    integer(is)             :: Bindx,Kindx
    integer(is)             :: nexci,bnopen,bnsp,i,a,kval1
    integer(is)             :: counter
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
       bnsp=bncsfs(bnopen)

       ! Get the indices of the MOs involved in the excitation
       hlist=0
       plist=0
       call get_exci_indices(kconf,bconf,n_int,hlist(1),plist(1),1)
       
       ! If we are working with - < Ket | T_ji^(1,-k) | Bra >,
       ! then use -k in the evaluation of the spin-coupling
       ! coefficients
       if (transpose) then
          kval1=-kval
       else
          kval1=kval
       endif

       ! Get the spin-coupling coefficients
       scc(1:knsp*bnsp)=triplet_scc(kval1,knsp,bnsp,ksop,plist(1),&
            hlist(1),knopen,knbefore)
                 
       ! Idices of the triplet TDM elements
       i=m2c(hlist(1))
       a=m2c(plist(1))

       ! Loop over triplet TDMs
       do ipair=1,npairs

          ! Bra and ket eigenvector indices
          Bindx=Bmap(ipair)
          Kindx=Kmap(ipair)

          ! Initialise the spin-coupling coefficient counter
          counter=0
          
          ! Loop over ket CSFs
          do ikcsf=kcsfs(ikconf),kcsfs(ikconf+1)-1
             kcoe=vecK(ikcsf,Kindx)
             
             ! Loop over bra CSFs
             do ibcsf=bcsfs(ibconf),bcsfs(ibconf+1)-1
                bcoe=vecB(ibcsf,Bindx)
                counter=counter+1
                
                ! Contribution to the 1-TDM
                prod=kcoe*bcoe*scc(counter)
                if (transpose) then
                   Tij(i,a,ipair)=Tij(i,a,ipair)-prod
                else
                   Tij(a,i,ipair)=Tij(a,i,ipair)+prod
                endif
                   
             enddo
             
          enddo
          
       enddo
       
    enddo
       
    return
    
  end subroutine triplet_tdm_batch
    
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
    use pattern_indices_k0
    use pattern_indices_kp1
    use pattern_indices_km1
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
       pattern=pattern_index_k0(sop,ac,ia,nc,na,nopen,icase)
    case(-1)
       pattern=pattern_index_km1(sop,ac,ia,nc,na,nopen,icase)
    case(1)
       pattern=pattern_index_kp1(sop,ac,ia,nc,na,nopen,icase)
    end select
    
!----------------------------------------------------------------------
! Fill in the array of spin-coupling coefficients
!----------------------------------------------------------------------
    select case(kval)
    case(0)
       triplet_scc=spincp(pattern:pattern+knsp*bnsp-1)
    case(-1)
       triplet_scc=spincp_minus(pattern:pattern+knsp*bnsp-1)
    case(1)
       triplet_scc=spincp_plus(pattern:pattern+knsp*bnsp-1)
    end select
       
    return

  end function triplet_scc
  
!######################################################################
  
end module triplet_tdm
