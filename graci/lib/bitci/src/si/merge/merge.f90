!**********************************************************************
! Routines for the merging the reference spaces of a pair of bra and
! ket configuration derived types
!**********************************************************************
! IMPORTANT: We currently do not account for 1E, 2E and 1I1E confs
!            becoming R confs in the unified ref space. This is fine
!            for, e.g., valence -> CVS TDMs, but needs to be changed
!            before SOC terms can be computed.
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
    integer(is)                :: n_int_I
    
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
    call rearrange_confs(cfgB,cfgK,nmoI,nmoE,Ilist,Elist,m2c,c2m)


    
    STOP

    
    
    return
    
  end subroutine merge_ref_space
  
!######################################################################
! sameref: returns true if the bra and ket reference spaces are the
!          same
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
! internal_union: given a pair of bra and ket configuration derived
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
    type(mrcfg), intent(in)  :: cfgB,cfgK

    ! Merged internal & external MO information
    integer(is), intent(out) :: nmoI,nmoE
    integer(is), intent(out) :: Ilist(nmo),Elist(nmo)
    integer(is), intent(out) :: m2c(nmo),c2m(nmo)

    ! Working arrays
    integer(is)              :: nmoI_B,nmoE_B
    integer(is)              :: nmoI_K,nmoE_K
    integer(is)              :: Ilist_B(nmo),Elist_B(nmo)
    integer(is)              :: Ilist_K(nmo),Elist_K(nmo)
    integer(is)              :: m2c_B(nmo),c2m_B(nmo)
    integer(is)              :: m2c_K(nmo),c2m_K(nmo)
    integer(ib)              :: iint(n_int),iext(n_int)
    integer(ib), allocatable :: confR_B(:,:,:),confR_K(:,:,:)
    
    ! Everything else
    integer(is)              :: i,k,n

!----------------------------------------------------------------------
! Bra and ket reference space configurations
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(confR_B(n_int,2,cfgB%nR))
    allocate(confR_K(n_int,2,cfgK%nR))
    confR_B=0_ib
    confR_K=0_ib

    ! Bra ref confs
    confR_B(1:cfgB%n_int_I,:,:)=cfgB%confR(:,:,:)

    ! Ket ref confs
    confR_K(1:cfgK%n_int_I,:,:)=cfgK%confR(:,:,:)
    
!----------------------------------------------------------------------
! Put the bra and ket configuration bit strings into the canonical MO
! ordering
!----------------------------------------------------------------------
    ! Bra
    call reorder_confs(cfgB%m2c,confR_B,cfgB%nR)
    
    ! Ket
    call reorder_confs(cfgK%m2c,confR_K,cfgK%nR)

!----------------------------------------------------------------------
! Compute the bit string encodings of the union of the bra and ket
! internal MOs
!----------------------------------------------------------------------
    ! Initialisation
    iint=0_ib

    ! Loop over bra configurations
    do i=1,cfgB%nR

       ! Fill in the iint bit string
       do k=1,n_int
          iint(k)=ior(iint(k),confR_B(k,1,i))
          iint(k)=ior(iint(k),confR_B(k,2,i))
       enddo
       
    enddo

    ! Loop over ket configurations
    do i=1,cfgK%nR

       ! Fill in the iint bit string
       do k=1,n_int
          iint(k)=ior(iint(k),confR_K(k,1,i))
          iint(k)=ior(iint(k),confR_K(k,2,i))
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
! rearrange_confs: rearranges the bra and ket configurations to
!                  correspond to the new internal-external MO
!                  partitioning
!######################################################################
  subroutine rearrange_confs(cfgB,cfgK,nmoI,nmoE,Ilist,Elist,m2c,c2m)

    use constants
    use bitglobal
    use conftype
    use mrciutils

    implicit none

    ! Bra and ket MRCI configuration derived types
    type(mrcfg), intent(inout) :: cfgB,cfgK

    ! Merged internal & external MO information
    integer(is), intent(in)    :: nmoI,nmoE
    integer(is), intent(in)    :: Ilist(nmo),Elist(nmo)
    integer(is), intent(in)    :: m2c(nmo),c2m(nmo)

    ! New configuration classes
    integer(is)                :: n1I_B,n2I_B,n1E_B,n2E_B,n1I1E_B
    integer(is)                :: n1I_K,n2I_K,n1E_K,n2E_K,n1I1E_K
    integer(is), allocatable   :: iBclass(:),iKclass(:)
        
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(iBclass(cfgB%confdim))
    ibclass=0
    
    allocate(iKclass(cfgK%confdim))
    ikclass=0
    
!----------------------------------------------------------------------
! Fill in the concatenated arrays of configuration bit strings
!----------------------------------------------------------------------
    ! Bra
    call cfgB%concatenate_arrays

    ! Ket
    call cfgK%concatenate_arrays

!----------------------------------------------------------------------
! Put the configuration bit strings into the new internal-external
! ordering
!----------------------------------------------------------------------
    ! Bra
    call reorder_confs(cfgB%m2c,cfgB%confall,cfgB%confdim)
    call reorder_confs(c2m,cfgB%confall,cfgB%confdim)

    ! Ket
    call reorder_confs(cfgK%m2c,cfgK%confall,cfgK%confdim)
    call reorder_confs(c2m,cfgK%confall,cfgK%confdim)

!----------------------------------------------------------------------
! Determine the number of configurations in each class
!----------------------------------------------------------------------
    ! Bra
    call class_dimensions(cfgB,cfgB%confall,cfgB%confdim,nmoI,nelB,&
         n1I_B,n2I_B,n1E_B,n2E_B,n1I1E_B,iBclass)

    ! Ket
    call class_dimensions(cfgK,cfgK%confall,cfgK%confdim,nmoI,nelK,&
         n1I_K,n2I_K,n1E_K,n2E_K,n1I1E_K,iKclass)

!----------------------------------------------------------------------
! Fill in the new configuration and offset arrays
!----------------------------------------------------------------------
    
    
    return
    
  end subroutine rearrange_confs
    
!######################################################################
! class_dimensions: determines the number of configurations in each
!                   class (1I, 1E, etc.) as well as the integer class
!                   labels of each:
!                   0 <-> R
!                   1 <-> 1I
!                   2 <-> 2I
!                   3 <-> 1E
!                   4 <-> 2E
!                   5 <-> 1I1E
!######################################################################
  subroutine class_dimensions(cfg,conf,confdim,nmoI,nelectrons,n1I,&
       n2I,n1E,n2E,n1I1E,iclass)

    use constants
    use bitglobal
    use conftype
    use iomod
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg
    
    ! Dimensions
    integer(is), intent(in)  :: confdim
    integer(is), intent(in)  :: nmoI,nelectrons
    integer(is), intent(out) :: n1I,n2I,n1E,n2E,n1I1E

    ! Configurations
    integer(ib), intent(in)  :: conf(n_int,2,confdim)

    ! Configuration class integer labels
    integer(is), intent(out) :: iclass(confdim)

    ! Everything else
    integer(is)              :: counter,iconf,next,k

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    iclass=-1
    counter=0

    n1I=cfg%n1I
    n2I=cfg%n2I
    n1E=0
    n2E=0
    n1I1E=0
    
!----------------------------------------------------------------------
! Reference configurations: immutable number of
!----------------------------------------------------------------------
    ! Loop over ref confs
    do iconf=1,cfg%n0h

       ! Increment the counter
       counter=counter+1

       ! Set the class integer label
       iclass(counter)=0
       
    enddo
    
!----------------------------------------------------------------------
! 1I configurations: immutable number of
!----------------------------------------------------------------------
    ! Loop over 1I confs
    do iconf=1,cfg%n1I

       ! Increment the counter
       counter=counter+1

       ! Set the class integer label
       iclass(counter)=1

    enddo
       
!----------------------------------------------------------------------
! 2I configurations: immutable number of
!----------------------------------------------------------------------
    ! Loop over 2I confs
    do iconf=1,cfg%n2I

       ! Increment the counter
       counter=counter+1

       ! Set the class integer label
       iclass(counter)=2

    enddo

!----------------------------------------------------------------------
! 1E configurations: possible change to 1I configurations (and Ref?)
!----------------------------------------------------------------------
    ! Loop over 1E confs
    do iconf=1,cfg%n1E
    
       ! Increment the counter
       counter=counter+1
    
       ! Number of external electrons
       next=nexternal(conf(:,:,counter),nmoI,nelectrons)
       
       ! Fill in the configuration class
       if (next == 1) then
          n1E=n1E+1
          iclass(counter)=3
       else if (next == 0) then
          n1I=n1I+1
          iclass(counter)=1
       else
          errmsg='Nonsensical 1E next value in class dimensions'
          call error_control
       endif
          
    enddo

!----------------------------------------------------------------------
! 2E configurations: possible change to 1I1E and 2I configurations
!                    (and Ref?)
!----------------------------------------------------------------------
    ! Loop over 2E confs
    do iconf=1,cfg%n2E
    
       ! Increment the counter
       counter=counter+1
    
       ! Number of external electrons
       next=nexternal(conf(:,:,counter),nmoI,nelectrons)

       ! Fill in the configuration class
       if (next == 2) then
          n2E=n2E+1
          iclass(counter)=4
       else if (next == 1) then
          n1I1E=n1I1E+1
          iclass(counter)=5
       else if (next == 0) then
          n2I=n2I+1
          iclass(counter)=2
       else
          errmsg='Nonsensical 2E next value in class dimensions'
          call error_control
       endif
          
    enddo

!----------------------------------------------------------------------
! 1I1E configurations: possible change to 2I configurations (and Ref?)
!----------------------------------------------------------------------
    ! Loop over 1I1E confs
    do iconf=1,cfg%n1I1E
    
       ! Increment the counter
       counter=counter+1
    
       ! Number of external electrons
       next=nexternal(conf(:,:,counter),nmoI,nelectrons)

       ! Fill in the configuration class
       if (next == 1) then
          n1I1E=n1I1E+1
          iclass(counter)=5
       else if (next == 0) then
          n2I=n2I+1
          iclass(counter)=2
       else
          errmsg='Nonsensical 1I1E next value in class dimensions'
          call error_control
       endif
          
    enddo

    return
    
  end subroutine class_dimensions

!######################################################################
! nexternal: given a configuration bit string, returns the number of
!            external electrons
!######################################################################
  function nexternal(conf,nmoI,nelectrons)

    use constants
    use bitglobal
    
    implicit none

    ! Function result
    integer(is)             :: nexternal

    ! Configuration bit string
    integer(ib), intent(in) :: conf(n_int,2)

    ! Number of internal MOs
    integer(is), intent(in) :: nmoI

    ! Number of electrons
    integer(is), intent(in) :: nelectrons
    
    ! Everything else
    integer(is)             :: i,k,n,n_int_I,ninternal
    integer(ib)             :: work(n_int,2)
    
    !
    ! Initialise the work array
    !
    work=conf

    !
    ! Number of 64 bit integers required to encode the internal
    ! MOs
    !
    n_int_I=(nmoI-1)/64+1

    !
    ! Unset the bits corresponding to the external MOs in the
    ! n_int_I'th block of the bit string
    !
    n=mod(nmoI,64)
    if (n /= 0) then
       work(n_int_I,1)=ibits(work(n_int_I,1),0,n)
       work(n_int_I,2)=ibits(work(n_int_I,2),0,n)
    endif

    !
    ! Count the number of internal electrons
    !
    ninternal=0
    do k=1,n_int_I
       ninternal=ninternal+popcnt(work(k,1))
       ninternal=ninternal+popcnt(work(k,2))
    enddo

    !
    ! Number of external electrons
    !
    nexternal=nelectrons-ninternal

    return
    
  end function nexternal
  
!######################################################################
  
end module merge
