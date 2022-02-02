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
  subroutine merge_ref_space(cfgB,cfgK,nspB,nspK)

    use constants
    use bitglobal
    use conftype
    use iomod
    
    implicit none

    ! Bra and ket MRCI configuration derived types
    type(mrcfg), intent(inout) :: cfgB,cfgK

    ! Numbers of bra and ket CSFs as a function of the number of
    integer(is), intent(in)    :: nspB(0:nocase2),nspK(0:nocase2)
    
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
! to the new internal-external partitioning. Here we will refill the
! bra and ket MRCI configuration derived types.
!----------------------------------------------------------------------
    call rearrange_confs(cfgB,cfgK,nmoI,nmoE,Ilist,Elist,m2c,c2m,&
         nspB,nspK)
    
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
  subroutine rearrange_confs(cfgB,cfgK,nmoI,nmoE,Ilist,Elist,m2c,c2m,&
       nspB,nspK)

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

    ! Numbers of bra and ket CSFs as a function of the number of
    integer(is), intent(in)    :: nspB(0:nocase2),nspK(0:nocase2)
    
    ! Sets of bra and ket ref confs that are not in the intersection
    ! of the two sets
    integer(is)                :: nUc_B,nUc_K
    integer(is), allocatable   :: Uc_B(:),Uc_K(:)
    
    ! New configuration classes
    integer(is)                :: n1I_B,n2I_B,n1E_B,n2E_B,n1I1E_B
    integer(is)                :: n1I_K,n2I_K,n1E_K,n2E_K,n1I1E_K
    integer(is), allocatable   :: iBclass(:),iKclass(:)
        
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(Uc_B(cfgB%n0h))
    Uc_B=0

    allocate(Uc_K(cfgB%n0h))
    Uc_K=0
    
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
! Determine the complement of the intersection of the sets of bra and
! ket reference configurations
!----------------------------------------------------------------------
    call new_ref_confs(cfgB,cfgK,cfgB%n0h,cfgK%n0h,Uc_B,Uc_K,nUc_B,&
         nUc_K)
    
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
! Save some memory by deallocating the concatenated configuration
! arrays now that we are done with them
!----------------------------------------------------------------------
    deallocate(cfgB%confall, cfgB%sopall, cfgB%csfsall)
    deallocate(cfgK%confall, cfgK%sopall, cfgK%csfsall)
    
!----------------------------------------------------------------------
! Fill in the new configuration and offset arrays
!----------------------------------------------------------------------
    ! Bra
    call refill_cfg(cfgB,cfgB%confdim,n1I_B,n2I_B,n1E_B,n2E_B,n1I1E_B,&
         iBclass,nmoI,nmoE,m2c,c2m,nspB)

    ! Ket
    call refill_cfg(cfgK,cfgK%confdim,n1I_K,n2I_K,n1E_K,n2E_K,n1I1E_K,&
         iKclass,nmoI,nmoE,m2c,c2m,nspK)
    
    return
    
  end subroutine rearrange_confs

!######################################################################
! new_ref_confs: determines the sets of bra and ket ref confs that are
!                *not* in the intersection of the two sets
!######################################################################
  subroutine new_ref_confs(cfgB,cfgK,n0h_B,n0h_K,Uc_B,Uc_K,nUc_B,nUc_K)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! Bra and ket MRCI configuration derived types
    type(mrcfg), intent(inout) :: cfgB,cfgK

    ! Dimensions of the bra and ket reference spaces
    integer(is), intent(in)    :: n0h_B,n0h_K
    
    ! Sets of bra and ket ref confs that are not in the intersection
    ! of the two sets
    integer(is), intent(out)   :: nUc_B,nUc_K
    integer(is), intent(out)   :: Uc_B(:),Uc_K(:)
    
    return
    
  end subroutine new_ref_confs
  
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
! Reference configurations: cannot change class
!----------------------------------------------------------------------
    ! Loop over ref confs
    do iconf=1,cfg%n0h

       ! Increment the counter
       counter=counter+1

       ! Set the class integer label
       iclass(counter)=0
       
    enddo
    
!----------------------------------------------------------------------
! 1I configurations: cannot change class
!----------------------------------------------------------------------
    ! Loop over 1I confs
    do iconf=1,cfg%n1I

       ! Increment the counter
       counter=counter+1

       ! Set the class integer label
       iclass(counter)=1

    enddo
       
!----------------------------------------------------------------------
! 2I configurations: cannot change class
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
          errmsg='Nonsensical 1E next value in class_dimensions'
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
          errmsg='Nonsensical 2E next value in class_dimensions'
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
          errmsg='Nonsensical 1I1E next value in class_dimensions'
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
! refill_cfg: refills a given MRCI configuration derived type with
!             the confs, SOPs, offsets, and creation operator indices
!             corresponding to the new, merged MO indexing
!######################################################################
  subroutine refill_cfg(cfg,confdim,n1I,n2I,n1E,n2E,n1I1E,iclass,&
       nmoI,nmoE,m2c,c2m,nsp)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use iomod
    
    implicit none

    ! MRCI configuration derived types
    type(mrcfg), intent(inout) :: cfg
    
    ! Dimensions
    integer(is), intent(in)    :: confdim,n1I,n2I,n1E,n2E,n1I1E

    ! Configuration class integer labels
    integer(is), intent(in)    :: iclass(confdim)

    ! Merged internal & external MO information
    integer(is), intent(in)    :: nmoI,nmoE
    integer(is), intent(in)    :: m2c(nmo),c2m(nmo)

    ! Numbers of CSFs as a function of the number of open shells
    ! (only needed to call cfg%initialise)
    integer(is), intent(in)    :: nsp(0:nocase2)
    
    ! I/O
    integer(is)                :: iscratch,scrnum
    character(len=255)         :: filename

    ! Working arrays
    integer(ib), allocatable   :: conf(:,:,:)
    integer(is), allocatable   :: ngen(:)
    integer(is), allocatable   :: offset(:)
    integer(is), allocatable   :: a1(:),a2(:,:)
    
    ! Everything else
    integer(is)                :: n_int_I
    integer(is)                :: n,ioff,counter,k,i,irrep
    integer(is)                :: counter1I,counter2I,counter1E,&
                                  counter2E,counter1I1E
    
!----------------------------------------------------------------------
! Open the new configuration scratch file
!----------------------------------------------------------------------
    ! Register the scratch file
    call scratch_name('mrciconf.merged',filename)
    call register_scratch_file(scrnum,filename)

    ! Open the scratch file
    iscratch=scrunit(scrnum)
    open(iscratch,file=filename,form='unformatted',status='unknown')
    
!----------------------------------------------------------------------
! Subspace dimensions
!----------------------------------------------------------------------
    n_int_I=(nmoI-1)/64+1

    write(iscratch) n_int_I
    write(iscratch) nmoI
    write(iscratch) nmoE
    write(iscratch) confdim
    write(iscratch) cfg%nR
    write(iscratch) cfg%n1h
    write(iscratch) cfg%n2h
    write(iscratch) cfg%n0h
    write(iscratch) n1I
    write(iscratch) n2I
    write(iscratch) n1E
    write(iscratch) n2E
    write(iscratch) n1I1E

!----------------------------------------------------------------------
! Ref confs across all irreps
!----------------------------------------------------------------------
    ! Fill in the work array
    allocate(conf(n_int,2,cfg%nR))
    conf=0_ib
    conf(1:cfg%n_int_I,:,:)=cfg%confR(:,:,:)

    ! Re-order to correspond to the new internal-external MO
    ! partitioning
    call reorder_confs(cfg%m2c,conf,cfg%nR)
    call reorder_confs(c2m,conf,cfg%nR)

    ! Write to disk
    write(iscratch) conf(1:n_int_I,:,:)

    ! Deallocate the work array
    deallocate(conf)

!----------------------------------------------------------------------
! 1-hole confs
!----------------------------------------------------------------------
    ! Fill in the work array
    allocate(conf(n_int,2,cfg%n1h))
    conf=0_ib
    conf(1:cfg%n_int_I,:,:)=cfg%conf1h(:,:,:)
    
    ! Re-order to correspond to the new internal-external MO
    ! partitioning
    call reorder_confs(cfg%m2c,conf,cfg%n1h)
    call reorder_confs(c2m,conf,cfg%n1h)
    
    ! Write to disk
    write(iscratch) conf(1:n_int_I,:,:)

    ! Deallocate the work array
    deallocate(conf)

!----------------------------------------------------------------------
! 1-hole annihilation operator indices
!----------------------------------------------------------------------
    ! Allocate the work array
    allocate(a1(cfg%n1h))
    a1=0
    
    ! Convert the annihilation operator indices to the new MO ordering
    a1=cfg%a1h
    do k=1,cfg%n1h
       i=cfg%m2c(a1(k))
       a1(k)=c2m(i)
    enddo

    ! Write to disk
    write(iscratch) a1
    
    ! Deallocate the work array
    deallocate(a1)
    
!----------------------------------------------------------------------
! 1-hole offsets
!----------------------------------------------------------------------
    write(iscratch) cfg%off1h
    
!----------------------------------------------------------------------
! 2-hole confs
!----------------------------------------------------------------------
    ! Fill in the work array
    allocate(conf(n_int,2,cfg%n2h))
    conf=0_ib
    conf(1:cfg%n_int_I,:,:)=cfg%conf2h(:,:,:)

    ! Re-order to correspond to the new internal-external MO
    ! partitioning
    call reorder_confs(cfg%m2c,conf,cfg%n2h)
    call reorder_confs(c2m,conf,cfg%n2h)

    ! Write to disk
    write(iscratch) conf(1:n_int_I,:,:)

    ! Deallocate the work array
    deallocate(conf)

!----------------------------------------------------------------------
! 2-hole annihilation operator indices
!----------------------------------------------------------------------
    ! Allocate the work array
    allocate(a2(2,cfg%n2h))
    a2=0
    
    ! Convert the annihilation operator indices to the new MO ordering
    a2=cfg%a2h
    do k=1,cfg%n2h
       i=cfg%m2c(a2(1,k))
       a2(1,k)=c2m(i)
       i=cfg%m2c(a2(2,k))
       a2(2,k)=c2m(i)
    enddo
    
    ! Write to disk
    write(iscratch) a2
    
    ! Deallocate the work array
    deallocate(a2)
    
!----------------------------------------------------------------------
! 2-hole offsets
!----------------------------------------------------------------------
    write(iscratch) cfg%off2h

!----------------------------------------------------------------------
! Ref confs for this irrep
!----------------------------------------------------------------------
    ! Fill in the work array    
    allocate(conf(n_int,2,cfg%n0h))
    conf=0_ib
    conf(1:cfg%n_int_I,:,:)=cfg%conf0h(:,:,:)
    
    ! Re-order to correspond to the new internal-external MO
    ! partitioning
    call reorder_confs(cfg%m2c,conf,cfg%n0h)
    call reorder_confs(c2m,conf,cfg%n0h)

    ! Write to disk
    write(iscratch) conf(1:n_int_I,:,:)

    ! Deallocate the work array
    deallocate(conf)
    
!----------------------------------------------------------------------
! 1I confs: these can arise due to old 1I and 1E confs
!----------------------------------------------------------------------
    if (n1I > 0) then
    
       ! Allocate work arrays
       allocate(a1(n1I), ngen(cfg%n1h), offset(cfg%n1h+1))
       a1=0; ngen=0; offset=0
    
       ! New 1I conf counter
       k=0

       !
       ! New 1I confs arising from old 1I and 1E confs
       !
       ! Initialise the conf counter
       counter1I=cfg%n0h

       counter1E=cfg%n0h+cfg%n1I+cfg%n2I

       ! Loop over 1-hole confs    
       do n=1,cfg%n1h
          
          if (cfg%n1I > 0) then
             
             ! Loop over old 1I confs generated by the 1-hole conf
             do ioff=cfg%off1I(n),cfg%off1I(n+1)-1
             
                ! Increment the conf counter
                counter1I=counter1I+1
                
                ! Are we at a new 1I conf?
                if (iclass(counter1I) == 1) then

                   ! Increment the new 1I conf counter
                   k=k+1
                   
                   ! No. new 1I confs generated by the 1-hole conf
                   ngen(n)=ngen(n)+1
                
                   ! Creation operator index
                   a1(k)=cfg%a1I(ioff)
                   
                endif
          
             enddo

          endif

          if (cfg%n1E > 0) then
             
             ! Loop over old 1E confs generated by the 1-hole conf
             do ioff=cfg%off1E(n),cfg%off1E(n+1)-1
             
                ! Increment the conf counter
                counter1E=counter1E+1
             
                ! Are we at a new 1I conf?
                if (iclass(counter1E) == 1) then
                   
                   ! Increment the new 1I conf counter
                   k=k+1
                
                   ! No. new 1I confs generated by the 1-hole conf
                   ngen(n)=ngen(n)+1
                
                   ! Creation operator index
                   a1(k)=cfg%a1E(ioff)
                   
                endif

             enddo
             
          endif
          
       enddo

       ! Convert the creation operator indices to the new MO ordering
       do k=1,n1I
          i=cfg%m2c(a1(k))
          a1(k)=c2m(i)
       enddo
    
       ! Fill in the new 1I offset array
       offset(1)=1
       do n=2,cfg%n1h+1
          offset(n)=offset(n-1)+ngen(n-1)
       enddo
       
       ! Write the creation operator indices and conf offsets to disk
       write(iscratch) a1
       write(iscratch) offset

       ! Deallocate work arrays
       deallocate(a1, ngen, offset)
       
    endif
       
!----------------------------------------------------------------------
! 2I confs: these can arise due to old 2I, 2E and 1I1E confs
!----------------------------------------------------------------------
    if (n2I > 0) then

       ! Allocate work arrays
       allocate(a2(2,n2I), ngen(cfg%n2h), offset(cfg%n2h+1))
       a2=0; ngen=0; offset=0

       ! New 2I conf counter
       k=0

       !
       ! New 2I confs arising from old 2I, 2E and 1I1E confs
       !
       ! Initialise the conf counters
       counter2I=cfg%n0h+cfg%n1I
       counter2E=cfg%n0h+cfg%n1I+cfg%n2I+cfg%n1E
       counter1I1E=cfg%n0h+cfg%n1I+cfg%n2I+cfg%n1E+cfg%n2E
          
       ! Loop over 2-hole confs    
       do n=1,cfg%n2h

          if (cfg%n2I > 0) then
             
             ! Loop over old 2I confs generated by the 2-hole conf
             do ioff=cfg%off2I(n),cfg%off2I(n+1)-1
                
                ! Increment the conf counter
                counter2I=counter2I+1

                ! Are we at a new 2I conf?
                if (iclass(counter2I) == 2) then
                   
                   ! Increment the new 2I conf counter
                   k=k+1

                   ! No. new 2I confs generated by the 2-hole conf
                   ngen(n)=ngen(n)+1
                
                   ! Creation operator indices
                   a2(:,k)=cfg%a2I(:,ioff)
                   
                endif
             
             enddo

          endif

          if (cfg%n2E > 0) then
          
             ! Loop over old 2E confs generated by the 2-hole conf
             do ioff=cfg%off2E(n),cfg%off2E(n+1)-1
                
                ! Increment the conf counter
                counter2E=counter2E+1
             
                ! Are we at a new 2I conf?
                if (iclass(counter2E) == 2) then
                   
                   ! Increment the new 2I conf counter
                   k=k+1
             
                   ! No. new 2I confs generated by the 2-hole conf
                   ngen(n)=ngen(n)+1
                   
                   ! Creation operator indices
                   a2(:,k)=cfg%a2E(:,ioff)
                   
                endif
                
             enddo
          
          endif
          
          if (cfg%n1I1E > 0) then
          
             ! Loop over old 1I1E confs generated by the 2-hole conf
             do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1
             
                ! Increment the conf counter
                counter1I1E=counter1I1E+1

                ! Are we at a new 2I conf?
                if (iclass(counter1I1E) == 2) then

                   ! Increment the new 2I conf counter
                   k=k+1
                
                   ! No. new 2I confs generated by the 2-hole conf
                   ngen(n)=ngen(n)+1
                   
                   ! Creation operator indices
                   a2(:,k)=cfg%a1I1E(:,ioff)
                
                endif
             
             enddo

          endif

       enddo
       
       ! Convert the creation operator indices to the new MO ordering
       do k=1,n2I
          i=cfg%m2c(a2(1,k))
          a2(1,k)=c2m(i)
          i=cfg%m2c(a2(2,k))
          a2(2,k)=c2m(i)
       enddo

       ! Fill in the new 2I offset array
       offset(1)=1
       do n=2,cfg%n2h+1
          offset(n)=offset(n-1)+ngen(n-1)
       enddo
    
       ! Write the creation operator indices and conf offsets to disk
       write(iscratch) a2
       write(iscratch) offset
              
       ! Deallocate work arrays
       deallocate(a2, ngen, offset)
       
    endif
    
!----------------------------------------------------------------------
! 1E confs: these can arise due to old 1E confs only
!----------------------------------------------------------------------
    if (n1E > 0) then

       ! Allocate work arrays
       allocate(a1(n1E), ngen(cfg%n1h), offset(cfg%n1h+1))
       a1=0; ngen=0; offset=0
    
       ! New 1E conf counter
       k=0

       !
       ! New 1E confs arising from old 1E confs
       !
       if (cfg%n1E > 0) then

          ! Initialise the conf counter
          counter1E=cfg%n0h+cfg%n1I+cfg%n2I

          ! Loop over 1-hole confs
          do n=1,cfg%n1h
             
             ! Loop over old 1E confs generated by the 1-hole conf
             do ioff=cfg%off1E(n),cfg%off1E(n+1)-1
             
                ! Increment the conf counter
                counter1E=counter1E+1

                ! Are we at a new 1E conf?
                if (iclass(counter1E) == 3) then

                   ! Increment the new 1E conf counter
                   k=k+1
                
                   ! No. new 1E confs generated by the 1-hole conf
                   ngen(n)=ngen(n)+1
                
                   ! Creation operator index
                   a1(k)=cfg%a1E(ioff)                 
                   
                endif
                   
             enddo

          enddo
                
       endif

       ! Convert the creation operator indices to the new MO ordering
       do k=1,n1E
          i=cfg%m2c(a1(k))
          a1(k)=c2m(i)
       enddo
    
       ! Fill in the new 1E offset array
       offset(1)=1
       do n=2,cfg%n1h+1
          offset(n)=offset(n-1)+ngen(n-1)
       enddo
    
       ! Write the creation operator indices and conf offsets to disk
       write(iscratch) a1
       write(iscratch) offset

       ! Deallocate work arrays
       deallocate(a1, ngen, offset)
       
    endif

!----------------------------------------------------------------------
! 2E confs: these can arise due to old 2E confs only
!----------------------------------------------------------------------
    if (n2E > 0) then

       ! Allocate work arrays
       allocate(a2(2,n2E), ngen(cfg%n2h), offset(cfg%n2h+1))
       a2=0; ngen=0; offset=0
    
       ! New 2E conf counter
       k=0

       !
       ! New 2E confs arising from old 2E confs
       !
       if (cfg%n2E > 0) then
          
          ! Initialise the conf counter
          counter2E=cfg%n0h+cfg%n1I+cfg%n2I+cfg%n1E

          ! Loop over 2-hole confs    
          do n=1,cfg%n2h
             
             ! Loop over old 2E confs generated by the 2-hole conf
             do ioff=cfg%off2E(n),cfg%off2E(n+1)-1
                
                ! Increment the conf counter
                counter2E=counter2E+1
             
                ! Are we at a new 2E conf?
                if (iclass(counter2E) == 4) then
             
                   ! Increment the new 2E conf counter
                   k=k+1
             
                   ! No. new 2E confs generated by the 2-hole conf
                   ngen(n)=ngen(n)+1
                   
                   ! Creation operator indices
                   a2(:,k)=cfg%a2E(:,ioff)
                   
                endif
                
             enddo
             
          enddo

       endif

       ! Convert the creation operator indices to the new MO ordering
       do k=1,n2E
          i=cfg%m2c(a2(1,k))
          a2(1,k)=c2m(i)
          i=cfg%m2c(a2(2,k))
          a2(2,k)=c2m(i)
       enddo

       ! Fill in the new 2E offset array
       offset(1)=1
       do n=2,cfg%n2h+1
          offset(n)=offset(n-1)+ngen(n-1)
       enddo
    
       ! Write the creation operator indices and conf offsets to disk
       write(iscratch) a2
       write(iscratch) offset
              
       ! Deallocate work arrays
       deallocate(a2, ngen, offset)
       
    endif

!----------------------------------------------------------------------
! 1I1E confs: these can arise due to old 2E and 1I1E confs
!----------------------------------------------------------------------
    if (n1I1E > 0) then

       ! Allocate work arrays
       allocate(a2(2,n1I1E), ngen(cfg%n2h), offset(cfg%n2h+1))
       a2=0; ngen=0; offset=0

       ! New 1I1E conf counter
       k=0

       !
       ! New 1I1E confs arising from old 2E confs
       !
       ! Initialise the conf counters
       counter2E=cfg%n0h+cfg%n1I+cfg%n2I+cfg%n1E
       counter1I1E=cfg%n0h+cfg%n1I+cfg%n2I+cfg%n1E+cfg%n2E
          
       ! Loop over 2-hole confs    
       do n=1,cfg%n2h

          if (cfg%n2E > 0) then
             
             ! Loop over old 2E confs generated by the 2-hole conf
             do ioff=cfg%off2E(n),cfg%off2E(n+1)-1
                
                ! Increment the conf counter
                counter2E=counter2E+1
                
                ! Are we at a new 1I1E conf?
                if (iclass(counter2E) == 5) then
                   
                   ! Increment the new 1I1E conf counter
                   k=k+1
             
                   ! No. new 1I1E confs generated by the 2-hole conf
                   ngen(n)=ngen(n)+1
                   
                   ! Creation operator indices
                   a2(:,k)=cfg%a2E(:,ioff)

                endif
                
             enddo

          endif

          if (cfg%n1I1E > 0) then
          
             ! Loop over old 1I1E confs generated by the 2-hole conf
             do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1

                ! Increment the conf counter
                counter1I1E=counter1I1E+1
                
                ! Are we at a new 1I1E conf?
                if (iclass(counter1I1E) == 5) then

                   ! Increment the new 1I1E conf counter
                   k=k+1

                   ! No. new 1I1E confs generated by the 2-hole conf
                   ngen(n)=ngen(n)+1
                
                   ! Creation operator indices
                   a2(:,k)=cfg%a1I1E(:,ioff)
                   
                endif
             
             enddo

          endif

       enddo

       ! Convert the creation operator indices to the new MO ordering
       do k=1,n1I1E
          i=cfg%m2c(a2(1,k))
          a2(1,k)=c2m(i)
          i=cfg%m2c(a2(2,k))
          a2(2,k)=c2m(i)
       enddo

       ! Fill in the new 1I1E offset array
       offset(1)=1
       do n=2,cfg%n2h+1
          offset(n)=offset(n-1)+ngen(n-1)
       enddo
    
       ! Write the creation operator indices and conf offsets to disk
       write(iscratch) a2
       write(iscratch) offset
       
       ! Deallocate work arrays
       deallocate(a2, ngen, offset)
       
    endif
    
!----------------------------------------------------------------------
! MO mapping arrays
!----------------------------------------------------------------------
    write(iscratch) m2c
    write(iscratch) c2m

!----------------------------------------------------------------------
! Number of CSFs as a function of the the number of open shells
!----------------------------------------------------------------------
    write(iscratch) nsp
    
!----------------------------------------------------------------------
! Close the new configuration scratch file
!----------------------------------------------------------------------
    close(iscratch)

!----------------------------------------------------------------------
! Delete and refill the MRCI configuration derived type
!----------------------------------------------------------------------
    ! Save the irrep number
    irrep=cfg%irrep

    ! Deallocate arrays
    call cfg%finalise

    ! Fill the derived type with the new configuration information
    call cfg%initialise(irrep,scrnum)
    
    return
    
  end subroutine refill_cfg
  
!######################################################################
  
end module merge
