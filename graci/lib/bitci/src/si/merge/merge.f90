!**********************************************************************
! Routines for the merging the reference spaces of a pair of bra and
! ket configuration derived types
!**********************************************************************
module merge

  implicit none

contains

!######################################################################
! merge_ref_space: given a pair of bra and ket configuration derived
!                  types, re-structures the conf and SOP bit strings
!                  to correspond to the union of the two reference
!                  spaces
!######################################################################
  subroutine merge_ref_space(cfgB,cfgK)

    use constants
    use bitglobal
    use conftype
    use iomod
    
    implicit none

    ! Bra and ket MRCI configuration derived types
    type(mrcfg), intent(inout) :: cfgB,cfgK

    ! Merged internal & external MO information
    integer(is)                :: nmoI,nmoE
    integer(is)                :: Ilist(nmo),Elist(nmo)
    integer(is)                :: m2c(nmo),c2m(nmo)
    
!----------------------------------------------------------------------
! Do nothing if the two reference spaces are the same
!----------------------------------------------------------------------
    if (sameref(cfgB,cfgK)) return

!----------------------------------------------------------------------
! Get the indices of the MOs in the union of the bra and ket reference
! spaces
!----------------------------------------------------------------------
    call internal_union(cfgB,cfgK,nmoI,nmoE,Ilist,Elist,m2c,c2m)

!----------------------------------------------------------------------
! Rearrange the bra and ket configuration bit strings to correspond
! to the new internal-external partitioning
!----------------------------------------------------------------------
    
    
    STOP

    
    return
    
  end subroutine merge_ref_space
  
!######################################################################

  function sameref(cfgB,cfgK)

    use constants
    use conftype
    
    implicit none

    ! Function result
    logical                 :: sameref

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfgB,cfgK
    
    ! Everything else
    integer                 :: n,k,i

!----------------------------------------------------------------------
! Different number ref space confs across all irreps <-> different ref
! spaces    
!----------------------------------------------------------------------
    if (cfgB%nR /= cfgK%nR) then
       sameref=.false.
       return
    endif

!----------------------------------------------------------------------
! Different values of n_int_I <-> different ref spaces
!----------------------------------------------------------------------
    if (cfgB%n_int_I /= cfgK%n_int_I) then
       sameref=.false.
       return
    endif
    
!----------------------------------------------------------------------
! Same number of ref space confs across all irreps: check for
! differences
!----------------------------------------------------------------------
    sameref=.true.

    ! Loop over all ref confs
    do n=1,cfgB%nR

       ! Different bra and ref confs?
       do i=1,2
          do k=1,cfgB%n_int_I
             if (cfgB%confR(k,i,n) /= cfgK%confR(k,i,n)) then
                sameref=.false.
                return
             endif
          enddo
       enddo
          
    enddo
    
    return
    
  end function sameref
  
!######################################################################
! internal_union: Given a pair of bra and ket configuration derived
!                 types, determines the union of the sets of internal
!                 MO indices
!######################################################################
  subroutine internal_union(cfgB,cfgK,nmoI,nmoE,Ilist,Elist,m2c,c2m)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use bitutils
    
    implicit none

    ! Bra and ket MRCI configuration derived types
    type(mrcfg), intent(inout) :: cfgB,cfgK

    ! Merged internal & external MO information
    integer(is), intent(out)   :: nmoI,nmoE
    integer(is), intent(out)   :: Ilist(nmo),Elist(nmo)
    integer(is), intent(out)   :: m2c(nmo),c2m(nmo)

    ! Working arrays
    integer(is)                :: nmoI_B,nmoE_B
    integer(is)                :: nmoI_K,nmoE_K
    integer(is)                :: Ilist_B(nmo),Elist_B(nmo)
    integer(is)                :: Ilist_K(nmo),Elist_K(nmo)
    integer(is)                :: m2c_B(nmo),c2m_B(nmo)
    integer(is)                :: m2c_K(nmo),c2m_K(nmo)
    integer(ib)                :: iint(n_int),iext(n_int)
    integer(ib), allocatable   :: conf0h_B(:,:,:),conf0h_K(:,:,:)
    
    ! Everything else
    integer(is)                :: i,k,n

!----------------------------------------------------------------------
! Bra and ket reference space configurations
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(conf0h_B(n_int,2,cfgB%n0h))
    allocate(conf0h_K(n_int,2,cfgK%n0h))
    conf0h_B=0_ib
    conf0h_K=0_ib

    ! Bra ref confs
    conf0h_B(1:cfgB%n_int_I,:,:)=cfgB%conf0h(:,:,:)

    ! Ket ref confs
    conf0h_K(1:cfgK%n_int_I,:,:)=cfgK%conf0h(:,:,:)
    
!----------------------------------------------------------------------
! Put the bra and ket configuration bit strings into the canonical MO
! ordering
!----------------------------------------------------------------------
    ! Bra
    call canonical_ordering(cfgB%m2c,conf0h_B,cfgB%n0h)
    
    ! Ket
    call canonical_ordering(cfgK%m2c,conf0h_K,cfgK%n0h)

!----------------------------------------------------------------------
! Compute the bit string encodings of the union of the bra and ket
! internal MOs
!----------------------------------------------------------------------
    ! Initialisation
    iint=0_ib

    ! Loop over bra configurations
    do i=1,cfgB%n0h

       ! Fill in the iint bit string
       do k=1,n_int
          iint(k)=ior(iint(k),conf0h_B(k,1,i))
          iint(k)=ior(iint(k),conf0h_B(k,2,i))
       enddo
       
    enddo

    ! Loop over ket configurations
    do i=1,cfgK%n0h

       ! Fill in the iint bit string
       do k=1,n_int
          iint(k)=ior(iint(k),conf0h_K(k,1,i))
          iint(k)=ior(iint(k),conf0h_K(k,2,i))
       enddo
       
    enddo

!----------------------------------------------------------------------
! Compute the bit string encodings of the complement of the union of
! the sets bra and ket internal MOs. That is, the merged external space
!----------------------------------------------------------------------
    ! Flip the bits in iint
    do k=1,n_int
       iext(k)=not(iint(k))
    enddo

    ! Unset any unused bits at the end of the array
    n=mod(nmo,64)
    if (n /= 0) iext(n_int)=ibits(iext(n_int),0,n)

!----------------------------------------------------------------------
! New numbers of internal and external MOs
!----------------------------------------------------------------------
    ! Internal MOs
    nmoI=0
    do k=1,n_int
       nmoI=nmoI+popcnt(iint(k))
    enddo
    
    ! External MOs
    nmoE=nmo-nmoI

!----------------------------------------------------------------------
! Get the lists of merged internal and external MOs
!----------------------------------------------------------------------
    ! Internal MOs
    Ilist=0
    call list_from_bitstring(iint,Ilist,nmo)
    
    ! External MOs
    Elist=0
    call list_from_bitstring(iext,Elist,nmo)

!----------------------------------------------------------------------
! Construct the merged canonical-to-MRCI MO index mapping array
!----------------------------------------------------------------------
    call get_mo_mapping(nmoI,nmoE,Ilist,Elist,m2c,c2m)
    
    return
    
  end subroutine internal_union
    
!######################################################################
  
end module merge
