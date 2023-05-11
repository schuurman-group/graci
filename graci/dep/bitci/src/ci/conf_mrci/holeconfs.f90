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
  subroutine generate_hole_confs(cfgM,icvs,nroots)
    
    use constants
    use bitglobal
    use conftype

    implicit none
    
    ! MRCI configuration derived types for all irreps
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)
    
    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical                    :: lcvs

    ! Number of roots per irrep
    integer(is), intent(in)    :: nroots(0:nirrep-1)
    
    ! Everything else
    integer(is)                :: modus,i,istart

!----------------------------------------------------------------------
! Is this a CVS-MRCI calculation
!----------------------------------------------------------------------
    if (sum(icvs) > 0) then
       lcvs=.true.
    else
       lcvs=.false.
    endif

!----------------------------------------------------------------------
! First index of an irrep with a non-zero number of roots
!----------------------------------------------------------------------
    do i=0,nirrep-1
       if (nroots(i) > 0) then
          istart=i
          exit
       endif
    enddo

!----------------------------------------------------------------------
! 1-hole configurations
!----------------------------------------------------------------------
    ! First pass: determine the no. 1-hole configurations
    modus=0
    call builder_1hole(modus,cfgM(istart),icvs,lcvs)

    ! Allocate and initialise arrays
    allocate(cfgM(istart)%conf1h(cfgM(istart)%n_int_I,2,cfgM(istart)%n1h))
    allocate(cfgM(istart)%off1h(cfgM(istart)%nR+1))
    allocate(cfgM(istart)%a1h(cfgM(istart)%n1h))
    cfgM(istart)%conf1h=0_ib
    cfgM(istart)%off1h=0
    cfgM(istart)%a1h=0

    ! Second pass: fill in the 1-hole configuration and offset arrays
    modus=1
    call builder_1hole(modus,cfgM(istart),icvs,lcvs)

!----------------------------------------------------------------------
! 2-hole configurations
!----------------------------------------------------------------------
    ! First pass: determine the no. 2-hole configurations
    modus=0
    call builder_2hole(modus,cfgM(istart),icvs,lcvs)
    
    ! Allocate and initialise arrays
    allocate(cfgM(istart)%conf2h(cfgM(istart)%n_int_I,2,cfgM(istart)%n2h))
    allocate(cfgM(istart)%off2h(cfgM(istart)%nR+1))
    allocate(cfgM(istart)%a2h(2,cfgM(istart)%n2h))
    cfgM(istart)%conf2h=0_ib
    cfgM(istart)%off2h=0
    cfgM(istart)%a2h=0

    ! Second pass: fill in the 2-hole configuration and offset arrays
    modus=1
    call builder_2hole(modus,cfgM(istart),icvs,lcvs)

!----------------------------------------------------------------------
! Fill in the MRCI configuration derived types for the remaining irreps
!----------------------------------------------------------------------
    ! Loop over remaining irreps
    do i=istart+1,nirrep-1

       ! Cycle if there are no roots for this irrep
       if (nroots(i) == 0) cycle
       
       ! No. 1-hole and 2-hole configurations
       cfgM(i)%n1h=cfgM(istart)%n1h
       cfgM(i)%n2h=cfgM(istart)%n2h

       ! 1-hole configurations and offsets
       allocate(cfgM(i)%conf1h(cfgM(i)%n_int_I,2,cfgM(i)%n1h))
       allocate(cfgM(i)%off1h(cfgM(i)%nR+1))
       allocate(cfgM(i)%a1h(cfgM(i)%n1h))
       cfgM(i)%conf1h=cfgM(istart)%conf1h
       cfgM(i)%off1h=cfgM(istart)%off1h
       cfgM(i)%a1h=cfgM(istart)%a1h

       ! 2-hole configurations and offsets
       allocate(cfgM(i)%conf2h(cfgM(i)%n_int_I,2,cfgM(i)%n2h))
       allocate(cfgM(i)%off2h(cfgM(i)%nR+1))
       allocate(cfgM(i)%a2h(2,cfgM(i)%n2h))
       cfgM(i)%conf2h=cfgM(istart)%conf2h
       cfgM(i)%off2h=cfgM(istart)%off2h
       cfgM(i)%a2h=cfgM(istart)%a2h
       
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
  
end module holeconfs
