!**********************************************************************
! Routines for the generation of hole configurations
!**********************************************************************
module holeconfs
  
  implicit none
  
contains

!######################################################################
! generate_hole_confs: generates all possible 1-hole and 2-hole
!                      configurations from a given set of reference
!                      space configurations
!######################################################################
  subroutine generate_hole_confs(cfgM,icvs)
    
    use constants
    use bitglobal
    use conftype

    implicit none
    
    ! MRCI configuration derived types for all irreps
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)
    
    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical                    :: lcvs

    ! Everything else
    integer(is)                :: modus,i

!----------------------------------------------------------------------
! Is this a CVS-MRCI calculation
!----------------------------------------------------------------------
    if (sum(icvs) > 0) then
       lcvs=.true.
    else
       lcvs=.false.
    endif

!----------------------------------------------------------------------
! 1-hole configurations
!----------------------------------------------------------------------
    ! First pass: determine the no. 1-hole configurations
    modus=0
    call builder_1hole(modus,cfgM(0),icvs,lcvs)

    ! Allocate and initialise arrays
    allocate(cfgM(0)%conf1h(cfgM(0)%n_int_I,2,cfgM(0)%n1h))
    allocate(cfgM(0)%off1h(cfgM(0)%nR+1))
    allocate(cfgM(0)%a1h(cfgM(0)%n1h))
    cfgM(0)%conf1h=0_ib
    cfgM(0)%off1h=0
    cfgM(0)%a1h=0

    ! Second pass: fill in the 1-hole configuration and offset arrays
    modus=1
    call builder_1hole(modus,cfgM(0),icvs,lcvs)

!----------------------------------------------------------------------
! 2-hole configurations
!----------------------------------------------------------------------
    ! First pass: determine the no. 2-hole configurations
    modus=0
    call builder_2hole(modus,cfgM(0),icvs,lcvs)
    
    ! Allocate and initialise arrays
    allocate(cfgM(0)%conf2h(cfgM(0)%n_int_I,2,cfgM(0)%n2h))
    allocate(cfgM(0)%off2h(cfgM(0)%nR+1))
    allocate(cfgM(0)%a2h(2,cfgM(0)%n2h))
    cfgM(0)%conf2h=0_ib
    cfgM(0)%off2h=0
    cfgM(0)%a2h=0

    ! Second pass: fill in the 2-hole configuration and offset arrays
    modus=1
    call builder_2hole(modus,cfgM(0),icvs,lcvs)

!----------------------------------------------------------------------
! Fill in the MRCI configuration derived types for the remaining irreps
!----------------------------------------------------------------------
    ! Loop over remaining irreps
    do i=1,nirrep-1

       ! No. 1-hole and 2-hole configurations
       cfgM(i)%n1h=cfgM(0)%n1h
       cfgM(i)%n2h=cfgM(0)%n2h

       ! 1-hole configurations and offsets
       allocate(cfgM(i)%conf1h(cfgM(i)%n_int_I,2,cfgM(i)%n1h))
       allocate(cfgM(i)%off1h(cfgM(i)%nR+1))
       allocate(cfgM(i)%a1h(cfgM(i)%n1h))
       cfgM(i)%conf1h=cfgM(0)%conf1h
       cfgM(i)%off1h=cfgM(0)%off1h
       cfgM(i)%a1h=cfgM(0)%a1h

       ! 2-hole configurations and offsets
       allocate(cfgM(i)%conf2h(cfgM(i)%n_int_I,2,cfgM(i)%n2h))
       allocate(cfgM(i)%off2h(cfgM(i)%nR+1))
       allocate(cfgM(i)%a2h(2,cfgM(i)%n2h))
       cfgM(i)%conf2h=cfgM(0)%conf2h
       cfgM(i)%off2h=cfgM(0)%off2h
       cfgM(i)%a2h=cfgM(0)%a2h
       
    enddo

    return
    
  end subroutine generate_hole_confs
    
!######################################################################
! builder_1hole: performs all the heavy lifting involved in the
!                generation of the 1-hole configurations
!######################################################################
  subroutine builder_1hole(modus,cfgM,icvs,lcvs)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use dethash
    use hparam
    
    implicit none

    ! Mode of operation: modus=0 <-> determine the number of hole
    !                                configurations
    !                    modus=1 <-> build all of the hole
    !                                configurations
    integer(is), intent(in)    :: modus
    
    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical, intent(in)        :: lcvs

    ! Lists of hole/particle indices linking ref confs
    integer(is), parameter     :: maxexci=1
    integer(is)                :: hlist(maxexci),plist(maxexci)

    ! Allowed annihilation operator indices
    integer(is)                :: iannihilate(nmo)

    ! Orbital classes
    integer(is)                :: socc(nmo),docc(nmo),unocc(nmo)
    integer(is)                :: nopen,nsocc,ndocc,nunocc

    ! 1-hole conf bit string
    integer(ib), allocatable   :: hconf(:,:)

    ! Number of 1-hole configurations generated by each
    ! reference configuration
    integer(is), allocatable   :: ngen(:)
    
    ! Everything else
    integer(is)                :: n_int_I,n,np,nexci,imo,i1
    integer(is)                :: counter

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfgM%n_int_I
    allocate(hconf(n_int_I,2))
    allocate(ngen(cfgM%nR))

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n1h=0
    counter=0
    ngen=0
    
!----------------------------------------------------------------------
! Generate the 1-hole configurations
!----------------------------------------------------------------------
    ! Loop over reference configurations
    do n=1,cfgM%nR

       ! Initialise the annihilation operator list
       iannihilate=1
       
       ! Loop over preceding reference configurations
       do np=1,n-1

          ! Determine the excitation degree between the to reference
          ! configurations
          nexci=exc_degree_conf(cfgM%confR(:,:,n),cfgM%confR(:,:,np),&
               n_int_I)

          ! Cycle if the excitation degree is greater than 1
          if (nexci > 1) cycle

          ! Determine the index of the annihilation operator and mark
          ! it as forbidden
          call get_exci_indices(cfgM%confR(:,:,n),cfgM%confR(:,:,np),&
               n_int_I,hlist,plist,1)
          
          iannihilate(hlist(1))=0

       enddo

       ! Get the lists of singly- and doubly-occupied MOs
       ! for the current reference configuration
       call sop_socc_list(cfgM%sopR(:,:,n),n_int_I,socc,nmo,nsocc)
       call sop_docc_list(cfgM%sopR(:,:,n),n_int_I,docc,nmo,ndocc)

       ! 1-hole configurations generated from annihilation of
       ! electrons in singly-occupied MOs
       do imo=1,nsocc

          ! MO index
          i1=socc(imo)

          ! Cycle if this is a CVS-MRCI calculation and we are creating
          ! a hole in a flagged core MO
          if (lcvs .and. icvs(cfgM%m2c(i1)) == 1) cycle
          
          ! Cycle if this annihilation operation will yield a
          ! duplicate 1-hole configuration
          if (iannihilate(i1) == 0) cycle

          ! Update the no. 1-hole configurations
          if (modus == 0) cfgM%n1h=cfgM%n1h+1

          ! Save the 1-hole configuration
          if (modus == 1) then
             counter=counter+1
             ngen(n)=ngen(n)+1
             hconf=annihilate_electron(cfgM%confR(:,:,n),n_int_I,i1)
             cfgM%conf1h(:,:,counter)=hconf
             cfgM%a1h(counter)=i1
          endif
          
       enddo

       ! 1-hole configurations generated from annihilation of
       ! electrons in doubly-occupied MOs
       do imo=1,ndocc

          ! MO index
          i1=docc(imo)

          ! Cycle if this is a CVS-MRCI calculation and we are creating
          ! a hole in a flagged core MO
          if (lcvs .and. icvs(cfgM%m2c(i1)) == 1) cycle
          
          ! Cycle if this annihilation operation will yield a
          ! duplicate 1-hole configuration
          if (iannihilate(i1) == 0) cycle

          ! Update the no. 1-hole configurations
          if (modus == 0) cfgM%n1h=cfgM%n1h+1

          ! Save the 1-hole configuration
          if (modus == 1) then
             counter=counter+1
             ngen(n)=ngen(n)+1
             hconf=annihilate_electron(cfgM%confR(:,:,n),n_int_I,i1)
             cfgM%conf1h(:,:,counter)=hconf
             cfgM%a1h(counter)=i1
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    if (modus == 1) then
       cfgM%off1h(1)=1
       do n=2,cfgM%nR+1
          cfgM%off1h(n)=cfgM%off1h(n-1)+ngen(n-1)
       enddo
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(hconf)
    deallocate(ngen)

    return
    
  end subroutine builder_1hole
    
!######################################################################
! builder_2hole: performs all the heavy lifting involved in the
!                generation of the 2-hole configurations
!######################################################################
  subroutine builder_2hole(modus,cfgM,icvs,lcvs)

    use constants
    use bitglobal
    use conftype
    use confinfo
    use mrciutils
    use dethash
    use hparam
    
    implicit none

    ! Mode of operation: modus=0 <-> determine the number of hole
    !                                configurations
    !                    modus=1 <-> build all of the hole
    !                                configurations
    integer(is), intent(in)    :: modus
    
    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical, intent(in)        :: lcvs
    
    ! 1-hole SOPs
    integer(ib), allocatable   :: sop1h(:,:,:)

    ! Lists of hole/particle indices linking ref confs
    integer(is), parameter     :: maxexci=2
    integer(is)                :: hlist(maxexci),plist(maxexci)
    
    ! Allowable annihilation operator indices
    integer(is)                :: iannihilate(nmo)
    
    ! Orbital classes
    integer(is)                :: socc(nmo),docc(nmo),unocc(nmo)
    integer(is)                :: nopen,nsocc,ndocc,nunocc

    ! Difference configuration information
    integer(is)                :: Dw(nmo,2)
    integer(is)                :: ndiff
    
    ! Number of 1-hole configurations generated by each
    ! reference configuration
    integer(is), allocatable   :: ngen(:)

     ! 1-hole conf bit string
    integer(ib), allocatable   :: hconf(:,:)
    
    ! Difference configuration information
    integer(is), allocatable   :: rholes(:,:)
    
    ! Everything else
    integer(is)                :: i,n,np,imo,i1,ioff,ia1h
    integer(is)                :: n_int_I,nexci,counter
    integer(is)                :: ndoccR
    integer(ib), allocatable   :: confI(:,:)
    logical                    :: allowed
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfgM%n_int_I
    allocate(hconf(n_int_I,2))
    hconf=0_ib
    
    allocate(sop1h(n_int_I,2,cfgM%n1h))
    sop1h=0_ib

    allocate(confI(n_int_I,2))
    confI=0_ib
    
    allocate(ngen(cfgM%nR))
    ngen=0

    allocate(rholes(3,cfgM%nR))
    rholes=0
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n2h=0
    counter=0

!----------------------------------------------------------------------
! Determine the number of doubly-occupied inactive MOs within the
! reference space
!----------------------------------------------------------------------
    ndoccR=ndocc_ref(cfgM)
    
!----------------------------------------------------------------------
! Generate the 1-hole SOPs
!----------------------------------------------------------------------    
    ! Loop over 1-hole configurations
    do i=1,cfgM%n1h

       ! Generate the next 1-hole SOP
       confI=cfgM%conf1h(:,:,i)
       sop1h(:,:,i)=conf_to_sop(confI,n_int_I)

    enddo

!----------------------------------------------------------------------
! Generate the 2-hole configurations
!----------------------------------------------------------------------
    ! Loop over reference configurations
    i=0
    do n=1,cfgM%nR

       ! Initialise the array of inter-ref-conf annihilation
       ! operator indices
       rholes=0
              
       ! Loop over preceding reference configurations
       do np=1,n-1

          ! Determine the excitation degree between the to reference
          ! configurations
          nexci=exc_degree_conf(cfgM%confR(:,:,n),cfgM%confR(:,:,np),&
               n_int_I)

          ! Cycle if the excitation degree is greater than 2
          if (nexci > 2) cycle

          ! Determine the indices of the creation/annihilation
          ! operators linking the two ref confs
          call get_exci_indices(cfgM%confR(:,:,n),cfgM%confR(:,:,np),&
               n_int_I,hlist(1:nexci),plist(1:nexci),nexci)

          ! Save the annihilation operator indices for use later
          ! in identifying duplicate 2-hole configurations
          rholes(1,np)=nexci
          if (nexci == 1) rholes(2,np)=hlist(1)
          if (nexci == 2) then
             rholes(2,np)=min(hlist(1),hlist(2))
             rholes(3,np)=max(hlist(1),hlist(2))
          endif
             
       enddo
          
       ! Loop over the 1-hole configurations generated by the current
       ! reference configuration
       do ioff=cfgM%off1h(n),cfgM%off1h(n+1)-1
          
          ! Increment the 1-hole conf counter
          i=i+1

          ! 1-hole annihilation operator index
          ia1h=cfgM%a1h(i)
          
          ! Initialise the annihilation operator list
          iannihilate=1
          
          ! Flag annihilation operator indices as forbidden based
          ! on the condition that the 2nd index has to be greater
          ! than or equal to the 1st index
          iannihilate(1:ia1h-1)=0
          
          ! Flag annihilation operator indices as forbidden based
          ! on the excitations linking the current ref conf to those
          ! preceding it
          do np=1,n-1
             if (rholes(1,np) == 1) then
                ! Single excitation
                iannihilate(rholes(2,np))=0
             else if (rholes(1,np) == 2) then
                ! Double excitation
                if (ia1h == rholes(2,np)) iannihilate(rholes(3,np))=0
             endif
          enddo
          
          ! Get the lists of singly- and doubly-occupied MOs
          ! for the current 1-hole configuration
          call sop_socc_list(sop1h(:,:,i),n_int_I,socc,nmo,nsocc)
          call sop_docc_list(sop1h(:,:,i),n_int_I,docc,nmo,ndocc)
          
          ! 2-hole configurations generated from annihilation of
          ! electrons in singly-occupied MOs
          do imo=1,nsocc
             
             ! MO index
             i1=socc(imo)

             ! Cycle if this is a CVS-MRCI calculation and we are creating
             ! a hole in a flagged core MO
             if (lcvs .and. icvs(cfgM%m2c(i1)) == 1) cycle
             
             ! Cycle if this annihilation operation will yield a
             ! duplicate 2-hole configuration
             if (iannihilate(i1) == 0) cycle
             
             ! Update the no. 2-hole configurations
             if (modus == 0) cfgM%n2h=cfgM%n2h+1

             ! Save the 2-hole configuration
             if (modus == 1) then
                counter=counter+1
                ngen(n)=ngen(n)+1
                hconf=annihilate_electron(cfgM%conf1h(:,:,i),n_int_I,i1)
                cfgM%conf2h(:,:,counter)=hconf
                cfgM%a2h(1,counter)=i1
                cfgM%a2h(2,counter)=ia1h
             endif
             
          enddo

          ! 2-hole configurations generated from annihilation of
          ! electrons in doubly-occupied MOs
          do imo=1,ndocc
             
             ! MO index
             i1=docc(imo)
             
             ! Cycle if this is a CVS-MRCI calculation and we are creating
             ! a hole in a flagged core MO
             if (lcvs .and. icvs(cfgM%m2c(i1)) == 1) cycle
             
             ! Cycle if this annihilation operation will yield a
             ! duplicate 2-hole configuration
             if (iannihilate(i1) == 0) cycle
             
             ! Update the no. 2-hole configurations
             if (modus == 0) cfgM%n2h=cfgM%n2h+1

             ! Save the 2-hole configuration
             if (modus == 1) then
                counter=counter+1
                ngen(n)=ngen(n)+1
                hconf=annihilate_electron(cfgM%conf1h(:,:,i),n_int_I,i1)
                cfgM%conf2h(:,:,counter)=hconf
                cfgM%a2h(1,counter)=i1
                cfgM%a2h(2,counter)=ia1h
             endif
             
          enddo
          
       enddo
          
    enddo

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    if (modus == 1) then
       cfgM%off2h(1)=1
       do n=2,cfgM%nR+1
          cfgM%off2h(n)=cfgM%off2h(n-1)+ngen(n-1)
       enddo
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(sop1h)
    deallocate(confI)
    deallocate(ngen)
    deallocate(rholes)
        
    return
    
  end subroutine builder_2hole
    
!######################################################################
! generate_1hole_1I_confs: generates 1-hole configurations derived
!                          from the 2-hole configurations by the
!                          application of internal creation
!                          operators
!######################################################################
  subroutine generate_1hole_1I_hash(conf1h1I,n1h1I,indx1h1I,cfgM,icvs)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use dethash
    use hparam
    use iomod
    
    implicit none

    ! 1H1I configurations
    integer(is), intent(out)   :: n1h1I
    integer(ib), allocatable   :: conf1h1I(:,:,:)
    integer(is), allocatable   :: indx1h1I(:,:)

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical                    :: lcvs

    ! Orbital classes
    integer(is), allocatable   :: socc(:),docc(:),unocc(:)
    integer(is)                :: nopen,nsocc,ndocc,nunocc

    ! Difference configuration information
    integer(is)                :: Dw(nmo,2)
    integer(is)                :: ndiff
    
    ! Hash table
    type(dhtbl)                :: h
    integer(is)                :: initial_size
    integer(ib)                :: key(n_int,2)
    integer(ib), allocatable   :: keyI(:,:)
    
    ! Storage the 1H1I configurations
    integer(is)                :: nold,nbuf,nrec
    integer(is), allocatable   :: ibuffer(:,:)
    integer(ib), allocatable   :: cbuffer(:,:,:)
    integer(is)                :: iscratch
    character(len=60)          :: buffile
    
    ! Everything else
    integer(is)                :: i,j,k,n,imo,i1
    integer(is)                :: n_int_I,nmoI,n1h,n2h
    
!----------------------------------------------------------------------
! Is this a CVS-MRCI calculation
!----------------------------------------------------------------------
    if (sum(icvs) > 0) then
       lcvs=.true.
    else
       lcvs=.false.
    endif
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    nmoI=cfgM%nmoI
    n_int_I=cfgM%n_int_I
    n1h=cfgM%n1h
    n2h=cfgM%n2h

    allocate(keyI(n_int_I,2))
    keyI=0_ib
    
    allocate(socc(nmoI),docc(nmoI),unocc(nmoI))
    socc=0; docc=0; unocc=0
    
    allocate(ibuffer(2,bufsize))
    ibuffer=0

    allocate(cbuffer(n_int_I,2,bufsize))
    cbuffer=0_ib

!----------------------------------------------------------------------
! Open the buffer file
!----------------------------------------------------------------------
    call freeunit(iscratch)
    call scratch_name('1h1Iconf',buffile)
    open(iscratch,file=buffile,form='unformatted',status='unknown')
    
!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------
    initial_size=min(n2h*(nel-1)/2,1024)
    call h%initialise_table(initial_size)

!----------------------------------------------------------------------
! Insert the 1-hole configurations
! These will be deleted before the keys are retrieved from the hash
! table, and serve only to filter out replications when the internal
! creation operators are applied to the 2-hole configurations
!----------------------------------------------------------------------
    ! Loop over 1-hole configurations
    do n=1,n1h

       ! Insert the 1-hole configuration
       key=0_ib
       key(1:n_int_I,:)=cfgM%conf1h(:,:,n)
       call h%insert_key(key)
       
    enddo

!----------------------------------------------------------------------
! Generate the 1-hole, 1I configurations
!----------------------------------------------------------------------
    ! Initialise the buffer variables
    nold=h%n_keys_stored
    nbuf=0
    nrec=0
    
    ! Loop over 2-hole configurations
    do n=1,n2h

       ! 1H1I configurations generated from the current
       ! 2-hole configuration
       do imo=1,nmoI

          ! Cycle if this is a CVS-MRCI calculation and we are creating
          ! an electron in a flagged core MO
          if (lcvs .and. icvs(cfgM%m2c(imo)) == 1) cycle
          
          ! Block index
          k=(imo-1)/n_bits+1
          
          ! Postion of the external MO within the kth block
          i=imo-(k-1)*n_bits-1          

          ! Cycle if this MO is doubly-occupied
          if (btest(cfgM%conf2h(k,2,n),i)) cycle

          ! Create the 1H1I configuration
          key=0_ib
          key(1:n_int_I,:)=cfgM%conf2h(:,:,n)
          if (btest(cfgM%conf2h(k,1,n),i)) then
             key(k,2)=ibset(key(k,2),i)
          else
             key(k,1)=ibset(key(k,1),i)
          endif

          ! Hash table insertion
          call h%insert_key(key)

          ! If this is a unique configuration, then save it
          if (h%n_keys_stored > nold) then
             keyI=key(1:n_int_I,:)
             call save_1h1I(nold,h%n_keys_stored,nbuf,n_int_I,&
                  ibuffer,cbuffer,iscratch,nrec,n,imo,keyI)
          endif
             
       enddo
       
    enddo

!----------------------------------------------------------------------
! Get the number of 1H1I configurations and allocate arrays
!----------------------------------------------------------------------
    ! Number of keys stored in the hash table minus those
    ! corresponding to the 1-hole configurations
    n1h1I=h%n_keys_stored-n1h

    ! Allocate arrays
    allocate(conf1h1I(n_int_I,2,n1h1I))
    conf1h1I=0_ib
    allocate(indx1h1I(2,n1h1I))
    indx1h1I=0
    
!----------------------------------------------------------------------
! Load the 1H1I configurations from the buffer
!----------------------------------------------------------------------
    call load_1h1I(indx1h1I,conf1h1I,n1h1I,n_int_I,ibuffer,cbuffer,&
         nrec,nbuf,iscratch)
    
!----------------------------------------------------------------------
! Delete the hash table
!----------------------------------------------------------------------
    call h%delete_table

!----------------------------------------------------------------------
! Close the buffer file
!----------------------------------------------------------------------
    close(iscratch)
 
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ibuffer)
    deallocate(cbuffer)
    deallocate(keyI)
    deallocate(socc,docc,unocc)
    
    return
    
  end subroutine generate_1hole_1I_hash
    
!######################################################################
! save_1h1I: buffered saving of the 1H1I configurations
!######################################################################
  subroutine save_1h1I(nold,n_keys_stored,nbuf,n_int_I,ibuffer,&
       cbuffer,iscratch,nrec,i2h,imo,conf)

    use bitglobal
    use constants
    
    implicit none

    integer(is), intent(out)   :: nold
    integer(ib), intent(in)    :: n_keys_stored
    integer(is), intent(inout) :: nbuf
    integer(is), intent(in)    :: n_int_I
    integer(is), intent(inout) :: ibuffer(2,bufsize)
    integer(ib), intent(inout) :: cbuffer(n_int_I,2,bufsize)
    integer(is), intent(in)    :: iscratch
    integer(is), intent(inout) :: nrec
    integer(is), intent(in)    :: i2h,imo
    integer(ib), intent(in)    :: conf(n_int_I,2)
    
    !
    ! Update the number of entries in the buffer
    !
    nbuf=nbuf+1

    !
    ! Number of confs generated up to this point
    !
    nold=n_keys_stored

    !
    ! Store the index-pair in the buffer
    !
    ibuffer(1,nbuf)=i2h
    ibuffer(2,nbuf)=imo

    !
    ! Store the configuration bit string in the buffer
    !
    cbuffer(:,:,nbuf)=conf
    
    !
    ! If the buffer is full, then write it to disk and reset
    ! the buffer
    !
    if (nbuf == bufsize) then

       ! Increment the record counter
       nrec=nrec+1

       ! Write the buffer to disk
       write(iscratch) ibuffer
       write(iscratch) cbuffer
       
       ! Reset the buffer
       nbuf=0
       ibuffer=0
       cbuffer=0_ib
       
    endif
    
    return
        
  end subroutine save_1h1I

!######################################################################
! load_1h1I: loads the stored 1H1I configurations into memory
!######################################################################
  subroutine load_1h1I(indx1h1I,conf1h1I,n1h1I,n_int_I,ibuffer,&
       cbuffer,nrec,nbuf,iscratch)

    use bitglobal
    use constants
    
    implicit none

    integer(is), intent(in)  :: n1h1I,nrec,nbuf,iscratch
    integer(is), intent(in)  :: n_int_I
    integer(is), intent(out) :: indx1h1I(2,n1h1I)
    integer(ib), intent(out) :: conf1h1I(n_int_I,2,n1h1I)
    integer(is), intent(in)  :: ibuffer(2,bufsize)
    integer(ib), intent(in)  :: cbuffer(n_int_I,2,bufsize)
    integer(is)              :: i,counter,lim1,lim2
    
    !
    ! Rewind to the start of the scratch file
    !
    rewind(iscratch)
    
    !
    ! Indices written to disk
    !
    counter=0
    do i=1,nrec
       counter=counter+1
       lim1=(counter-1)*bufsize+1
       lim2=counter*bufsize
       read(iscratch) indx1h1I(:,lim1:lim2)
       read(iscratch) conf1h1I(:,:,lim1:lim2)
    enddo

    !
    ! Remaining indices in the buffer
    !
    lim1=nrec*bufsize+1
    lim2=nrec*bufsize+nbuf
    indx1h1I(:,lim1:lim2)=ibuffer(:,1:nbuf)
    conf1h1I(:,:,lim1:lim2)=cbuffer(:,:,1:nbuf)
    
    return
    
  end subroutine load_1h1I

!######################################################################
  
end module holeconfs
