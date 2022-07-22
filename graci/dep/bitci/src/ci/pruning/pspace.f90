!**********************************************************************
! Routines for the contruction of the P-space of configurations
!**********************************************************************
module pspace

  implicit none

contains

!######################################################################
! pspace_conf_indices: Determination of the indices of the
!                      configurations that generate CSFs corresponding
!                      above threshold A-vector elements
!######################################################################
  subroutine pspace_conf_indices(cfg,Athrsh,Avec,csfdim,confdim,nroots,&
       nvec,vecmap,i1I,i2I,i1E,i2E,i1I1E,n1I,n2I,n1E,n2E,n1I1E,&
       nexplicit,iexplicit)

    use constants
    use bitglobal
    use conftype
    use utils
    
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: csfdim,confdim,nroots,nvec
    integer(is), intent(in)  :: n1I,n2I,n1E,n2E,n1I1E
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! Configuration selection threshold
    real(dp), intent(in)     :: Athrsh
    
    ! A-vectors
    real(dp), intent(in)     :: Avec(csfdim,nvec)

    ! Indices of the A-vectors of interest
    integer(is), intent(in)  :: vecmap(nroots)
    
    ! Surviving configuration flags
    integer(is), intent(out) :: i1I(n1I),i2I(n2I),i1E(n1E),i2E(n2E),&
                                i1I1E(n1I1E)

    ! Explicitly added CSFs
    integer(is), intent(in)  :: nexplicit
    integer(is), intent(in)  :: iexplicit(csfdim)
    
    ! Working array
    real(dp), allocatable    :: work(:)
    
    ! Everything else
    integer(is)              :: i,n,n1,ioff,csf,root
    integer(is), allocatable :: indx(:),iok(:)
    integer(is)              :: refdim
    real(dp)                 :: sumsq,goal

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    allocate(indx(csfdim), iok(csfdim))

    allocate(work(csfdim))
    
    i1I=0; i2I=0; i1E=0; i2E=0; i1I1E=0

!----------------------------------------------------------------------
! Determine the indices of the dominant CSFs
!----------------------------------------------------------------------
    ! Dimension of the ref space
    refdim=cfg%csfs0h(cfg%n0h+1)-1

    iok=0
    
    ! Loop over roots
    do n1=1,nroots
       n=vecmap(n1)

       ! Makesure that the first refdim elements are zeroed
       work=0.0d0
       work(refdim+1:csfdim)=Avec(refdim+1:csfdim,n)
       
       ! Sort the A-vector elements in order of decreasing absolute
       ! value
       call dsortindxa1('D',csfdim,abs(work),indx)

       ! Truncation threshold
       goal=Athrsh*dot_product(work,work)

       ! Determine which CSFs need to be included to reach the
       ! desired squared A-vector norm
       sumsq=0.0d0
       do i=1,csfdim
          
          sumsq=sumsq+work(indx(i))**2
          
          iok(indx(i))=1
          
          if (sumsq > goal) exit

       enddo
       
    enddo

!----------------------------------------------------------------------
! Add in the explicitly included CSFs
!----------------------------------------------------------------------
    do i=1,nexplicit
       if (iexplicit(i) == 1) iok(i)=1
    enddo
    
!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
    if (n1I > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1I configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1I(n),cfg%off1I(n+1)-1

             ! Loop over the CSFs generated by this configuration
             do csf=cfg%csfs1I(ioff),cfg%csfs1I(ioff+1)-1

                if (iok(csf) == 1) i1I(ioff)=1
                
             enddo

          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
    if (n2I > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2I configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2I(n),cfg%off2I(n+1)-1

             ! Loop over the CSFs generated by this configuration
             do csf=cfg%csfs2I(ioff),cfg%csfs2I(ioff+1)-1
                
                if (iok(csf) == 1) i2I(ioff)=1
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
    if (n1E > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1E configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1E(n),cfg%off1E(n+1)-1

             ! Loop over the CSFs generated by this configuration
             do csf=cfg%csfs1E(ioff),cfg%csfs1E(ioff+1)-1

                if (iok(csf) == 1) i1E(ioff)=1
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
    if (n2E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2E(n),cfg%off2E(n+1)-1

             ! Loop over the CSFs generated by this configuration
             do csf=cfg%csfs2E(ioff),cfg%csfs2E(ioff+1)-1

                if (iok(csf) == 1) i2E(ioff)=1
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
    if (n1I1E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 1I1E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1

             ! Loop over the CSFs generated by this configuration
             do csf=cfg%csfs1I1E(ioff),cfg%csfs1I1E(ioff+1)-1

                if (iok(csf) == 1) i1I1E(ioff)=1
                
             enddo
             
          enddo

       enddo

    endif

    return

  end subroutine pspace_conf_indices

!######################################################################
! pspace_set_new_confs: Fills in an MRCI configuration derived type
!                       with the surviving configuration information
!######################################################################
  subroutine pspace_set_new_confs(cfg,cfg_new,i1I,i2I,i1E,i2E,i1I1E,&
       n1I,n2I,n1E,n2E,n1I1E,confscr,nconf)

    use constants
    use bitglobal
    use conftype
    use filter_confs
    use iomod
    
    implicit none

    ! MRCI configuration derived types
    type(mrcfg), intent(in)  :: cfg
    type(mrcfg), intent(out) :: cfg_new

    ! Dimensions
    integer(is), intent(in) :: n1I,n2I,n1E,n2E,n1I1E
    
    ! Surviving configuration flags
    integer(is), intent(out) :: i1I(n1I),i2I(n2I),i1E(n1E),i2E(n2E),&
                                i1I1E(n1I1E)

    ! Configuration scratch file number
    integer(is), intent(in)  :: confscr

    ! Total number of surviving configurations
    integer(is), intent(out) :: nconf
    
    ! Numbers of surviving configurations
    integer(is)              :: n1I_new,n2I_new,n1E_new,n2E_new,&
                                n1I1E_new

    ! Everything else
    integer(is)              :: i,n_int_I,n0h,ntotal
    integer(is)              :: iscratch

!----------------------------------------------------------------------
! Numbers of surviving configurations
!----------------------------------------------------------------------
    n1I_new=sum(i1I)
    n2I_new=sum(i2I)
    n1E_new=sum(i1E)
    n2E_new=sum(i2E)
    n1I1E_new=sum(i1I1E)

!----------------------------------------------------------------------    
! Add the reference space configuration information
!----------------------------------------------------------------------
    ! Dimensions
    n_int_I=cfg%n_int_I
    n0h=cfg%n0h
    
    ! Total number of ref confs across all irreps
    cfg_new%nR=cfg%nR

    ! Allocate arrays
    allocate(cfg_new%conf0h(n_int_I,2,n0h))
    allocate(cfg_new%sop0h(n_int_I,2,n0h))
    allocate(cfg_new%confR(n_int_I,2,cfg_new%nR))
    allocate(cfg_new%sopR(n_int_I,2,cfg_new%nR))
    allocate(cfg_new%m2c(nmo))
    allocate(cfg_new%c2m(nmo))
    
    ! Irrep number
    cfg_new%irrep=cfg%irrep
    
    ! No. ref. confs of the current irrep
    cfg_new%n0h=n0h
    
    ! MO subspace dimensions
    cfg_new%n_int_I=n_int_I
    cfg_new%nmoI=cfg%nmoI
    cfg_new%nmoE=cfg%nmoE
    
    ! Ref confs and SOPs of the current irrep
    cfg_new%conf0h=cfg%conf0h
    cfg_new%sop0h=cfg%sop0h
    
    ! Ref confs and SOPs across all irreps
    cfg_new%confR=cfg%confR
    cfg_new%sopR=cfg%sopR
    
    ! MO mapping arrays
    cfg_new%m2c=cfg%m2c
    cfg_new%c2m=cfg%c2m

!----------------------------------------------------------------------
! Hole configurations
!----------------------------------------------------------------------
    ! 1-hole confs
    cfg_new%n1h=cfg%n1h
    
    allocate(cfg_new%conf1h(n_int_I,2,cfg_new%n1h))
    allocate(cfg_new%off1h(cfg_new%nR+1))
    allocate(cfg_new%a1h(cfg_new%n1h))
    cfg_new%conf1h=cfg%conf1h
    cfg_new%off1h=cfg%off1h
    cfg_new%a1h=cfg%a1h
    
    ! 2-hole confs
    cfg_new%n2h=cfg%n2h
    allocate(cfg_new%conf2h(n_int_I,2,cfg_new%n2h))
    allocate(cfg_new%off2h(cfg_new%nR+1))
    allocate(cfg_new%a2h(2,cfg_new%n2h))
    cfg_new%conf2h=cfg%conf2h
    cfg_new%off2h=cfg%off2h
    cfg_new%a2h=cfg%a2h
    
!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
    cfg_new%n1I=n1I_new
    
    allocate(cfg_new%a1I(cfg_new%n1I))
    allocate(cfg_new%off1I(cfg_new%n1h+1))

    if (n1I_new > 0) call set_new_confs_1I(cfg,cfg_new,i1I,n1I)
    
!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
    cfg_new%n2I=n2I_new

    allocate(cfg_new%a2I(2,cfg_new%n2I))
    allocate(cfg_new%off2I(cfg_new%n2h+1))

    if (n2I_new > 0) call set_new_confs_2I(cfg,cfg_new,i2I,n2I)
    
!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
    cfg_new%n1E=n1E_new

    allocate(cfg_new%a1E(cfg_new%n1E))
    allocate(cfg_new%off1E(cfg_new%n1h+1))

    if (n1E_new > 0) call set_new_confs_1E(cfg,cfg_new,i1E,n1E)
    
!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
    cfg_new%n2E=n2E_new

    allocate(cfg_new%a2E(2,cfg_new%n2E))
    allocate(cfg_new%off2E(cfg_new%n2h+1))

    if (n2E_new > 0) call set_new_confs_2E(cfg,cfg_new,i2E,n2E)
    
!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
    cfg_new%n1I1E=n1I1E_new
    
    allocate(cfg_new%a1I1E(2,cfg_new%n1I1E))
    allocate(cfg_new%off1I1E(cfg_new%n2h+1))

    if (n1I1E_new > 0) call set_new_confs_1I1E(cfg,cfg_new,i1I1E,n1I1E)

!----------------------------------------------------------------------
! Remove any hole configurations which do not generate any full
! configurations
!----------------------------------------------------------------------
    call filter_hole_confs(cfg_new)

!----------------------------------------------------------------------
! Set the total number of configurations
!----------------------------------------------------------------------
    ntotal=cfg_new%n0h+cfg_new%n1I+cfg_new%n1E+cfg_new%n2E&
         +cfg_new%n1I1E+cfg_new%n2I
    
    nconf=ntotal
    
!----------------------------------------------------------------------
! Write the configurations to disk
!----------------------------------------------------------------------
    ! Open the scratch file
    iscratch=scrunit(confscr)
    open(iscratch,file=scrname(confscr),form='unformatted',&
         status='unknown')

    ! Subspace dimensions
    write(iscratch) cfg_new%n_int_I
    write(iscratch) cfg_new%nmoI
    write(iscratch) cfg_new%nmoE
    write(iscratch) ntotal
    write(iscratch) cfg_new%nR
    write(iscratch) cfg_new%n1h
    write(iscratch) cfg_new%n2h
    write(iscratch) cfg_new%n0h
    write(iscratch) cfg_new%n1I
    write(iscratch) cfg_new%n2I
    write(iscratch) cfg_new%n1E
    write(iscratch) cfg_new%n2E
    write(iscratch) cfg_new%n1I1E

    ! Configuration information
    write(iscratch) cfg_new%confR
    write(iscratch) cfg_new%conf1h
    write(iscratch) cfg_new%a1h
    write(iscratch) cfg_new%off1h
    write(iscratch) cfg_new%conf2h
    write(iscratch) cfg_new%a2h
    write(iscratch) cfg_new%off2h
    write(iscratch) cfg_new%conf0h
    if (cfg_new%n1I > 0) then
       write(iscratch) cfg_new%a1I
       write(iscratch) cfg_new%off1I
    endif
    if (cfg_new%n2I > 0) then
       write(iscratch) cfg_new%a2I
       write(iscratch) cfg_new%off2I
    endif
    if (cfg_new%n1E > 0) then
       write(iscratch) cfg_new%a1E
       write(iscratch) cfg_new%off1E
    endif
    if (cfg_new%n2E > 0) then
       write(iscratch) cfg_new%a2E
       write(iscratch) cfg_new%off2E
    endif
    if (cfg_new%n1I1E > 0 ) then
       write(iscratch) cfg_new%a1I1E
       write(iscratch) cfg_new%off1I1E
    endif
    
    ! MO mapping arrays
    write(iscratch) cfg_new%m2c
    write(iscratch) cfg_new%c2m

    ! Number of CSFs as a function of the the number of open shells
    write(iscratch) ncsfs
    
    ! Close the scratch file
    close(iscratch)

    return
    
  end subroutine pspace_set_new_confs
    
!######################################################################
! set_new_confs_1I: Fills in the new 1I conf informations
!######################################################################
  subroutine set_new_confs_1I(cfg,cfg_new,i1I,n1I)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! MRCI configuration derived types
    type(mrcfg), intent(in)    :: cfg
    type(mrcfg), intent(inout) :: cfg_new

    ! Surviving configurations indices
    integer(is), intent(in)  :: n1I
    integer(is), intent(in)  :: i1I(n1I)
    
    ! Everything else
    integer(is)              :: n1h,n,ioff,counter
    integer(is), allocatable :: ngen(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n1h=cfg_new%n1h

    allocate(ngen(n1h))
    ngen=0

!----------------------------------------------------------------------
! Fill in the new 1I configuration information
!----------------------------------------------------------------------
    counter=0

    ! Loop over 1-hole configurations
    do n=1,n1h

       ! Loop over the 1I configurations generated by the 1-hole
       ! configuration
       do ioff=cfg%off1I(n),cfg%off1I(n+1)-1

          if (i1I(ioff) == 1) then

             ! Update the surviving configuration counter
             counter=counter+1

             ! Update the no. surviving confs generated by the
             ! current hole conf
             ngen(n)=ngen(n)+1

             ! Fill in the creation operator index
             cfg_new%a1I(counter)=cfg%a1I(ioff)
             
          endif

       enddo
       
    enddo

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    cfg_new%off1I(1)=1
    do n=2,n1h+1
       cfg_new%off1I(n)=cfg_new%off1I(n-1)+ngen(n-1)
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ngen)

    return
    
  end subroutine set_new_confs_1I
    
!######################################################################
! set_new_confs_2I: Fills in the new 2I conf informations
!######################################################################
  subroutine set_new_confs_2I(cfg,cfg_new,i2I,n2I)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! MRCI configuration derived types
    type(mrcfg), intent(in)    :: cfg
    type(mrcfg), intent(inout) :: cfg_new

    ! Surviving configurations indices
    integer(is), intent(in)  :: n2I
    integer(is), intent(in)  :: i2I(n2I)
    
    ! Everything else
    integer(is)              :: n2h,n,ioff,counter
    integer(is), allocatable :: ngen(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n2h=cfg_new%n2h

    allocate(ngen(n2h))
    ngen=0

!----------------------------------------------------------------------
! Fill in the new 2I configuration information
!----------------------------------------------------------------------
    counter=0

    ! Loop over 2-hole configurations
    do n=1,n2h

       ! Loop over the 2I configurations generated by the 2-hole
       ! configuration
       do ioff=cfg%off2I(n),cfg%off2I(n+1)-1

          if (i2I(ioff) == 1) then

             ! Update the surviving configuration counter
             counter=counter+1

             ! Update the no. surviving confs generated by the
             ! current hole conf
             ngen(n)=ngen(n)+1

             ! Fill in the creation operator index
             cfg_new%a2I(:,counter)=cfg%a2I(:,ioff)
             
          endif

       enddo
       
    enddo

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    cfg_new%off2I(1)=1
    do n=2,n2h+1
       cfg_new%off2I(n)=cfg_new%off2I(n-1)+ngen(n-1)
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ngen)

    return
    
  end subroutine set_new_confs_2I

!######################################################################
! set_new_confs_1E: Fills in the new 1E conf informations
!######################################################################
  subroutine set_new_confs_1E(cfg,cfg_new,i1E,n1E)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! MRCI configuration derived types
    type(mrcfg), intent(in)    :: cfg
    type(mrcfg), intent(inout) :: cfg_new

    ! Surviving configurations indices
    integer(is), intent(in)  :: n1E
    integer(is), intent(in)  :: i1E(n1E)
    
    ! Everything else
    integer(is)              :: n1h,n,ioff,counter
    integer(is), allocatable :: ngen(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n1h=cfg_new%n1h

    allocate(ngen(n1h))
    ngen=0

!----------------------------------------------------------------------
! Fill in the new 1E configuration information
!----------------------------------------------------------------------
    counter=0

    ! Loop over 1-hole configurations
    do n=1,n1h

       ! Loop over the 1E configurations generated by the 1-hole
       ! configuration
       do ioff=cfg%off1E(n),cfg%off1E(n+1)-1

          if (i1E(ioff) == 1) then

             ! Update the surviving configuration counter
             counter=counter+1

             ! Update the no. surviving confs generated by the
             ! current hole conf
             ngen(n)=ngen(n)+1

             ! Fill in the creation operator index
             cfg_new%a1E(counter)=cfg%a1E(ioff)
             
          endif

       enddo
       
    enddo

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    cfg_new%off1E(1)=1
    do n=2,n1h+1
       cfg_new%off1E(n)=cfg_new%off1E(n-1)+ngen(n-1)
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ngen)

    return
    
  end subroutine set_new_confs_1E

!######################################################################
! set_new_confs_2E: Fills in the new 2E conf informations
!######################################################################
  subroutine set_new_confs_2E(cfg,cfg_new,i2E,n2E)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! MRCI configuration derived types
    type(mrcfg), intent(in)    :: cfg
    type(mrcfg), intent(inout) :: cfg_new

    ! Surviving configurations indices
    integer(is), intent(in)  :: n2E
    integer(is), intent(in)  :: i2E(n2E)
    
    ! Everything else
    integer(is)              :: n2h,n,ioff,counter
    integer(is), allocatable :: ngen(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n2h=cfg_new%n2h

    allocate(ngen(n2h))
    ngen=0

!----------------------------------------------------------------------
! Fill in the new 2E configuration information
!----------------------------------------------------------------------
    counter=0

    ! Loop over 2-hole configurations
    do n=1,n2h

       ! Loop over the 2E configurations generated by the 2-hole
       ! configuration
       do ioff=cfg%off2E(n),cfg%off2E(n+1)-1

          if (i2E(ioff) == 1) then

             ! Update the surviving configuration counter
             counter=counter+1

             ! Update the no. surviving confs generated by the
             ! current hole conf
             ngen(n)=ngen(n)+1

             ! Fill in the creation operator index
             cfg_new%a2E(:,counter)=cfg%a2E(:,ioff)
             
          endif

       enddo
       
    enddo

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    cfg_new%off2E(1)=1
    do n=2,n2h+1
       cfg_new%off2E(n)=cfg_new%off2E(n-1)+ngen(n-1)
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ngen)

    return
    
  end subroutine set_new_confs_2E

!######################################################################
! set_new_confs_1I1E: Fills in the new 1I1E conf informations
!######################################################################
  subroutine set_new_confs_1I1E(cfg,cfg_new,i1I1E,n1I1E)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! MRCI configuration derived types
    type(mrcfg), intent(in)    :: cfg
    type(mrcfg), intent(inout) :: cfg_new

    ! Surviving configurations indices
    integer(is), intent(in)  :: n1I1E
    integer(is), intent(in)  :: i1I1E(n1I1E)
    
    ! Everything else
    integer(is)              :: n2h,n,ioff,counter
    integer(is), allocatable :: ngen(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n2h=cfg_new%n2h

    allocate(ngen(n2h))
    ngen=0

!----------------------------------------------------------------------
! Fill in the new 1I1E configuration information
!----------------------------------------------------------------------
    counter=0

    ! Loop over 2-hole configurations
    do n=1,n2h

       ! Loop over the 1I1E configurations generated by the 2-hole
       ! configuration
       do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1

          if (i1I1E(ioff) == 1) then

             ! Update the surviving configuration counter
             counter=counter+1

             ! Update the no. surviving confs generated by the
             ! current hole conf
             ngen(n)=ngen(n)+1

             ! Fill in the creation operator index
             cfg_new%a1I1E(:,counter)=cfg%a1I1E(:,ioff)
             
          endif

       enddo
       
    enddo

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    cfg_new%off1I1E(1)=1
    do n=2,n2h+1
       cfg_new%off1I1E(n)=cfg_new%off1I1E(n-1)+ngen(n-1)
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ngen)

    return
    
  end subroutine set_new_confs_1I1E
  
!######################################################################
  
end module pspace