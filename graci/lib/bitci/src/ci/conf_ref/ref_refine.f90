!**********************************************************************
! Refinement of the MRCI or MR-ENPT2 reference space
!**********************************************************************
module ref_refine

  public refine_ref_space, refine_ref_space_pt2
  
contains
  
!######################################################################
! refine_ref_space: refinement of the reference space based on the
!                   dominant configurations of the MRCI eigenvectors
!######################################################################
#ifdef CBINDING
  subroutine refine_ref_space(confscrM,confscrR,vecscr,nroots,nextra,&
       cthrsh,minrnorm,ndconf) bind(c,name="refine_ref_space")
#else
  subroutine refine_ref_space(confscrM,confscrR,vecscr,nroots,nextra,&
       cthrsh,minrnorm,ndconf)
#endif
    
    use constants
    use bitglobal
    use conftype
    use refconf
    use mrciutils

    implicit none
    
    ! MRCI configuration scratch file numbers
    integer(is), intent(in)  :: confscrM(0:nirrep-1)

    ! Reference space configuration scratch file numbers
    integer(is), intent(out) :: confscrR(0:nirrep-1)
    
    ! MRCI eigenpair scratch file numbers
    integer(is), intent(in)  :: vecscr(0:nirrep-1)
    
    ! Number of roots per irrep
    integer(is), intent(in)  :: nroots(0:nirrep-1)

    ! Number of extra states per irrep
    integer(is), intent(in)  :: nextra(0:nirrep-1)
    
    ! Configuration selection threshold
    real(dp), intent(in)     :: cthrsh

    ! Minimum reference space norm
    real(dp), intent(out)    :: minrnorm

    ! New number of reference configurations per irrep
    integer(is), intent(out) :: ndconf(0:nirrep-1)
    
    ! MRCI configuration derived type
    type(mrcfg), allocatable :: cfg(:)

    ! Dominant configurations
    integer(is), allocatable :: id(:)
    integer(ib), allocatable :: dconf(:,:,:)
    integer(ib), allocatable :: dsop(:,:,:)

    ! Updated internal MO space information
    integer(is)              :: nmoI,nmoE
    integer(is)              :: Ilist(nmo),Elist(nmo)

    ! Old canonical-to-MRCI MO mapping
    integer(is)              :: m2c_old(nmo)
    
    ! Updated canonical-to-MRCI MO mapping
    integer(is)              :: m2c_new(nmo),c2m_new(nmo)

    ! Selection threshold damping
    integer(is)              :: nsurvive
    real(dp)                 :: scale,thrsh
    real(dp), parameter      :: damp=0.95d0
    
    ! Everything else
    integer(is)              :: irrep,nconf,ndconf_tot,i,start
    integer(is)              :: old,sumd,istart,iend
    integer(is)              :: rdim
    real(dp)                 :: rnorm

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Configuration derived types
    allocate(cfg(0:nirrep-1))

!----------------------------------------------------------------------
! Set up the MRCI configuration derived types
!----------------------------------------------------------------------
    do irrep=0,nirrep-1
       call cfg(irrep)%initialise(irrep,confscrM(irrep))
    enddo

!----------------------------------------------------------------------
! Total number of MRCI configurations
!----------------------------------------------------------------------
    nconf=0
    do irrep=0,nirrep-1
       nconf=nconf+cfg(irrep)%confdim
    enddo

!----------------------------------------------------------------------
! Minimum reference space norm
!----------------------------------------------------------------------
    minrnorm=1.0d0

    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Loop over roots for the current irrep
       do i=1,nroots(irrep)

          ! Norm of the wavefunction projected onto the reference
          ! space
          call refnorm(i,cfg(irrep),vecscr(irrep),rnorm)

          ! Update the minimum reference space norm
          if (rnorm < minrnorm) minrnorm=rnorm
          
       enddo
       
    enddo
    
!----------------------------------------------------------------------
! Determine the indices of the above-threshold configurations
!----------------------------------------------------------------------
    allocate(id(nconf))
    id=0

    ! Initialise the starting point in the id array
    start=1

    ! Initialise the above-threshold configuration counter
    old=0

    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Initialise the no. surviving CSFs
       nsurvive=0

       ! Adjust the selection threshold until enough CSFs
       ! have been retained
       scale=1.0d0
       thrsh=cthrsh
       do while(nsurvive < nroots(irrep)+nextra(irrep))
              
          ! Loop over roots for the current irrep
          do i=1,nroots(irrep)

             ! Selection threshold
             thrsh=thrsh*scale
             
             ! Fill in the indices of the above-threshold configurations
             ! for this root
             call fill_above_threshold(id,nconf,cfg(irrep),&
                  start,vecscr(irrep),i,thrsh,nsurvive)

             ! Update the scaling factor
             scale=scale*damp
             
          enddo

       enddo
          
       ! Update the no. above-threshold configurations
       sumd=sum(id)
       ndconf(irrep)=sumd-old       
       old=sumd
              
       ! Update the starting point in the id array
       start=start+cfg(irrep)%confdim
       
    enddo

!----------------------------------------------------------------------
! Get the above threshold configuration bit strings
!----------------------------------------------------------------------
    ! Total no. dominant/above-threshold configurations
    ndconf_tot=sum(id)

    ! Allocate the dconf array
    allocate(dconf(n_int,2,ndconf_tot))
    dconf=0_ib

    ! Initialise the starting point in the id array
    start=1
    
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Start and end points in the dconf array for this irrep
       istart=sum(ndconf(0:irrep-1))+1
       iend=istart+ndconf(irrep)-1

       ! Get the above threshold configurations for this irrep
       call get_dominant_confs(cfg(irrep),dconf(:,:,istart:iend),&
            ndconf(irrep),id,nconf,start)

       ! Update the starting point in the id array
       start=start+cfg(irrep)%confdim
       
    enddo

!----------------------------------------------------------------------
! Rearrange the new reference space configurations to correspond to
! canonical MO ordering
!----------------------------------------------------------------------
    ! Old MO mapping array
    m2c_old=cfg(0)%m2c

    ! Put the configurations into canonical ordering
    call reorder_confs(m2c_old,dconf,ndconf_tot)
    
!----------------------------------------------------------------------
! Update the internal-external MO spaces
!----------------------------------------------------------------------
    ! Get the lists of internal and external MOs
    call get_internal_external_mos(dconf,ndconf_tot,nmoI,nmoE,Ilist,&
         Elist)

    ! Construct the new canonical-to-MRCI MO index mapping array
    call get_mo_mapping(nmoI,nmoE,Ilist,Elist,m2c_new,c2m_new)

!----------------------------------------------------------------------
! Get the SOPs corresponding to the above threshold configurations
!----------------------------------------------------------------------
    ! Allocate the SOP array
    allocate(dsop(n_int,2,ndconf_tot))
    dsop=0_ib

    ! Compute the SOPs
    do i=1,ndconf_tot
       dsop(:,:,i)=conf_to_sop(dconf(:,:,i),n_int)
    enddo

!----------------------------------------------------------------------
! Re-arrange the new reference configurations such that the internal
! MOs come before the external MOs
!----------------------------------------------------------------------
    call rearrange_ref_confs(c2m_new,dconf,dsop,ndconf_tot)
    
!----------------------------------------------------------------------
! Write the new reference configuration files
!----------------------------------------------------------------------
    call write_ref_confs(dconf,dsop,ndconf_tot,m2c_new,c2m_new,&
         nmoI,nmoE,ndconf,confscrR)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(cfg)
    deallocate(id)
    deallocate(dconf)
    deallocate(dsop)
        
    return
  
  end subroutine refine_ref_space

!######################################################################
! refine_ref_space_pt2: refinement of the reference space based on the
!                       dominant configurations of the MR-ENPT2
!                       eigenvectors
!######################################################################
#ifdef CBINDING
  subroutine refine_ref_space_pt2(confscrM,confscrR,vecscr,Qscr,&
       nroots,nextra,cmin,alpha,beta,minrnorm,ndconf) &
       bind(c,name="refine_ref_space_pt2")
#else
  subroutine refine_ref_space_pt2(confscrM,confscrR,vecscr,Qscr,&
       nroots,nextra,cmin,alpha,beta,minrnorm,ndconf)
#endif

    use constants
    use bitglobal
    use conftype
    use refconf
    use mrciutils
    use iomod
    
    ! MRCI configuration scratch file numbers
    integer(is), intent(in)  :: confscrM(0:nirrep-1)
    
    ! Reference space configuration scratch file numbers
    integer(is), intent(out) :: confscrR(0:nirrep-1)
    
    ! MRCI eigenpair scratch file numbers
    integer(is), intent(in)  :: vecscr(0:nirrep-1)

    ! Q-space info scratch file numbers
    integer(is), intent(in)  :: Qscr(0:nirrep-1)
    
    ! Number of roots per irrep
    integer(is), intent(in)  :: nroots(0:nirrep-1)

    ! Number of extra states per irrep
    integer(is), intent(in)  :: nextra(0:nirrep-1)
    
    ! Configuration selection parameters
    real(dp), intent(in)     :: cmin,alpha,beta
  
    ! Minimum reference space norm
    real(dp), intent(out)    :: minrnorm
  
    ! New number of reference configurations per irrep
    integer(is), intent(out) :: ndconf(0:nirrep-1)
    
    ! MRCI configuration derived type
    type(mrcfg), allocatable :: cfg(:)
    
    ! Dominant configurations
    integer(is), allocatable :: id(:)
    integer(ib), allocatable :: dconf(:,:,:)
    integer(ib), allocatable :: dsop(:,:,:)
    
    ! Updated internal MO space information
    integer(is)              :: nmoI,nmoE
    integer(is)              :: Ilist(nmo),Elist(nmo)

    ! Old canonical-to-MRCI MO mapping
    integer(is)              :: m2c_old(nmo)
    
    ! Updated canonical-to-MRCI MO mapping
    integer(is)              :: m2c_new(nmo),c2m_new(nmo)

    ! Selection threshold damping
    integer(is)              :: nsurvive
    real(dp)                 :: scale
    real(dp), parameter      :: damp=0.950
    
    ! Everything else
    integer(is)              :: iscratch,maxroots
    integer(is)              :: irrep,nconf,ndconf_tot,i,n,start
    integer(is)              :: old,sumd,istart,iend
    integer(is)              :: rdim
    real(dp)                 :: rnorm,thrsh
    real(dp), allocatable    :: Qnorm(:,:)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Configuration derived types
    allocate(cfg(0:nirrep-1))

    ! Q-space norms
    maxroots=maxval(nroots)    
    allocate(Qnorm(maxroots,0:nirrep-1))
    Qnorm=0.0d0
    
!----------------------------------------------------------------------
! Set up the MRCI configuration derived types
!----------------------------------------------------------------------
    do irrep=0,nirrep-1
       call cfg(irrep)%initialise(irrep,confscrM(irrep))
    enddo

!----------------------------------------------------------------------
! Read the norms of the 1st-order wave function corrections from disk
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Open scratch file
       iscratch=scrunit(Qscr(irrep))
       open(iscratch,file=scrname(Qscr(irrep)),form='unformatted',&
            status='old')

       ! No. roots
       read(iscratch) n

       ! Sanity check
       if (n /= nroots(irrep)) then
          errmsg='Error reading Qinfo file: wrong no. roots'
          call error_control
       endif

       ! Q-space norms
       read(iscratch) Qnorm(1:nroots(irrep),irrep)
       
       ! Close scratch file
       close(iscratch)
       
    enddo
    
!----------------------------------------------------------------------
! Total number of MRCI configurations
!----------------------------------------------------------------------
    nconf=0
    do irrep=0,nirrep-1
       nconf=nconf+cfg(irrep)%confdim
    enddo

!----------------------------------------------------------------------
! Minimum reference space norm
!----------------------------------------------------------------------
    minrnorm=1.0d0

    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Loop over roots for the current irrep
       do i=1,nroots(irrep)

          ! Norm of the wavefunction projected onto the reference
          ! space
          call refnorm(i,cfg(irrep),vecscr(irrep),rnorm)

          ! Update the minimum reference space norm
          if (rnorm < minrnorm) minrnorm=rnorm
          
       enddo
       
    enddo
    
!----------------------------------------------------------------------
! Determine the indices of the above-threshold configurations
!----------------------------------------------------------------------
    allocate(id(nconf))
    id=0

    ! Initialise the starting point in the id array
    start=1

    ! Initialise the above-threshold configuration counter
    old=0

    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Initialise the no. surviving CSFs
       nsurvive=0

       ! Adjust the selection threshold until enough CSFs
       ! have been retained
       scale=1.0d0
       do while(nsurvive < nroots(irrep)+nextra(irrep))
       
          ! Loop over roots for the current irrep
          do i=1,nroots(irrep)

             ! Selection threshold for this root
             thrsh=max(cmin, alpha/(cosh(beta*Qnorm(i,irrep))**2))
             thrsh=thrsh*scale
             
             ! Fill in the indices of the above-threshold configurations
             ! for this root
             call fill_above_threshold(id,nconf,cfg(irrep),start,&
                  vecscr(irrep),i,thrsh,nsurvive)

             ! Update the scaling factor
             scale=scale*damp
             
          enddo

       enddo
          
       ! Update the no. above-threshold configurations
       sumd=sum(id)
       ndconf(irrep)=sumd-old       
       old=sumd
              
       ! Update the starting point in the id array
       start=start+cfg(irrep)%confdim
       
    enddo

!----------------------------------------------------------------------
! Get the above threshold configuration bit strings
!----------------------------------------------------------------------
    ! Total no. dominant/above-threshold configurations
    ndconf_tot=sum(id)

    ! Allocate the dconf array
    allocate(dconf(n_int,2,ndconf_tot))
    dconf=0_ib

    ! Initialise the starting point in the id array
    start=1
    
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Start and end points in the dconf array for this irrep
       istart=sum(ndconf(0:irrep-1))+1
       iend=istart+ndconf(irrep)-1

       ! Get the above threshold configurations for this irrep
       call get_dominant_confs(cfg(irrep),dconf(:,:,istart:iend),&
            ndconf(irrep),id,nconf,start)

       ! Update the starting point in the id array
       start=start+cfg(irrep)%confdim
       
    enddo

!----------------------------------------------------------------------
! Rearrange the new reference space configurations to correspond to
! canonical MO ordering
!----------------------------------------------------------------------
    ! Old MO mapping array
    m2c_old=cfg(0)%m2c

    ! Put the configurations into canonical ordering
    call reorder_confs(m2c_old,dconf,ndconf_tot)
    
!----------------------------------------------------------------------
! Update the internal-external MO spaces
!----------------------------------------------------------------------
    ! Get the lists of internal and external MOs
    call get_internal_external_mos(dconf,ndconf_tot,nmoI,nmoE,Ilist,&
         Elist)

    ! Construct the new canonical-to-MRCI MO index mapping array
    call get_mo_mapping(nmoI,nmoE,Ilist,Elist,m2c_new,c2m_new)

!----------------------------------------------------------------------
! Get the SOPs corresponding to the above threshold configurations
!----------------------------------------------------------------------
    ! Allocate the SOP array
    allocate(dsop(n_int,2,ndconf_tot))
    dsop=0_ib

    ! Compute the SOPs
    do i=1,ndconf_tot
       dsop(:,:,i)=conf_to_sop(dconf(:,:,i),n_int)
    enddo

!----------------------------------------------------------------------
! Re-arrange the new reference configurations such that the internal
! MOs come before the external MOs
!----------------------------------------------------------------------
    call rearrange_ref_confs(c2m_new,dconf,dsop,ndconf_tot)
    
!----------------------------------------------------------------------
! Write the new reference configuration files
!----------------------------------------------------------------------
    call write_ref_confs(dconf,dsop,ndconf_tot,m2c_new,c2m_new,&
         nmoI,nmoE,ndconf,confscrR)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(cfg)
    deallocate(id)
    deallocate(dconf)
    deallocate(dsop)
    
    return
  
  end subroutine refine_ref_space_pt2
  
!######################################################################
! refnorm: calculation of the norm of a single wavefunction projected
!          onto the reference space
!######################################################################
  subroutine refnorm(sindx,cfg,vecscr,rnorm)

    use constants
    use bitglobal
    use conftype
    use iomod
    
    implicit none

    ! State index
    integer(is), intent(in) :: sindx
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! MRCI eigenvector scratch file number
    integer(is), intent(in) :: vecscr

    ! Reference space norm
    real(dp), intent(out)   :: rnorm

    ! Eigenvector
    real(dp), allocatable   :: vec(:)

    ! Everything else
    integer(is)             :: csfdim,rdim

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    csfdim=cfg%csfdim
    allocate(vec(csfdim))
    
!----------------------------------------------------------------------
! Read in the eigenvector
!----------------------------------------------------------------------
    call read_single_vector(vecscr,vec,csfdim,sindx)

!----------------------------------------------------------------------
! Reference space norm
!----------------------------------------------------------------------
    rdim=cfg%csfs0h(cfg%n0h+1)-1
    rnorm=sqrt(dot_product(vec(1:rdim),vec(1:rdim)))
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(vec)
    
    return
    
  end subroutine refnorm
  
!######################################################################
! fill_above_threshold: determines the configurations with absolute
!                       coefficient values above the threshold cthrsh
!######################################################################
  subroutine fill_above_threshold(id,nconf,cfg,start,vecscr,sindx,&
       cthrsh,nsurvive)

    use constants
    use bitglobal
    use conftype
    use utils
    use iomod
    
    implicit none

    ! Total number of configurations
    integer(is), intent(in)    :: nconf

    ! Above-threshold configuration flags
    integer(is), intent(inout) :: id(nconf)
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in)    :: cfg

    ! Starting point in the id array
    integer(is), intent(in)    :: start

    ! MRCI eigenvector scratch file number
    integer(is), intent(in)    :: vecscr

    ! State index
    integer(is), intent(in)    :: sindx

    ! Configuration threshold
    real(dp), intent(in)       :: cthrsh

    ! Number of surviving CSFs
    integer(is), intent(out)   :: nsurvive
    
    ! Eigenvector
    real(dp), allocatable      :: vec(:)
    
    ! Everything else
    integer(is)                :: csfdim,iconf,n,csf,ioff
    integer(is)                :: counter
    integer(is)                :: n_int_I
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    csfdim=cfg%csfdim
    allocate(vec(csfdim))
    
!----------------------------------------------------------------------
! Read in the eigenvector
!----------------------------------------------------------------------
    call read_single_vector(vecscr,vec,csfdim,sindx)
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ! Configuration counter
    iconf=start-1

    ! No. surviving CSFs
    nsurvive=0
    
!----------------------------------------------------------------------
! Reference space configurations
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I

    ! Loop over ref space configurations
    do n=1,cfg%n0h

       ! Increment the configuration counter
       iconf=iconf+1
       
       ! Loop over the CSFs generated by this configuration
       do csf=cfg%csfs0h(n),cfg%csfs0h(n+1)-1
          
          ! Cycle if this is not a dominant CSF
          if (abs(vec(csf)) < cthrsh) cycle
          
          ! Flag the configuration as being dominant
          id(iconf)=1
          
       enddo

       ! Update the no. surviving CSFs
       if (id(iconf) == 1) nsurvive=nsurvive+ &
            cfg%csfs0h(n+1)-cfg%csfs0h(n)
       
    enddo

!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
    if (cfg%n1I > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1I configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1I(n),cfg%off1I(n+1)-1

             ! Increment the configuration counter
             iconf=iconf+1
             
             ! Loop over the CSFs generated by this configuration
             do csf=cfg%csfs1I(ioff),cfg%csfs1I(ioff+1)-1

                ! Cycle if this is not a dominant CSF
                if (abs(vec(csf)) < cthrsh) cycle
                
                ! Flag the configuration as being dominant
                id(iconf)=1

             enddo

             ! Update the no. surviving CSFs
             if (id(iconf) == 1) nsurvive=nsurvive+ &
                  cfg%csfs1I(ioff+1)-cfg%csfs1I(ioff)
             
          enddo
          
       enddo

    endif

!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
    if (cfg%n2I > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2I configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2I(n),cfg%off2I(n+1)-1

             ! Increment the configuration counter
             iconf=iconf+1

             ! Loop over the CSFs generated by this configuration
             do csf=cfg%csfs2I(ioff),cfg%csfs2I(ioff+1)-1

                ! Cycle if this is not a dominant CSF
                if (abs(vec(csf)) < cthrsh) cycle
                
                ! Flag the configuration as being dominant
                id(iconf)=1

             enddo

             ! Update the no. surviving CSFs
             if (id(iconf) == 1) nsurvive=nsurvive+ &
                  cfg%csfs2I(ioff+1)-cfg%csfs2I(ioff)
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
    if (cfg%n1E > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1E configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1E(n),cfg%off1E(n+1)-1

             ! Increment the configuration counter
             iconf=iconf+1

             ! Loop over the CSFs generated by this configuration
             do csf=cfg%csfs1E(ioff),cfg%csfs1E(ioff+1)-1

                ! Cycle if this is not a dominant CSF
                if (abs(vec(csf)) < cthrsh) cycle
                
                ! Flag the configuration as being dominant
                id(iconf)=1
                
             enddo

             ! Update the no. surviving CSFs
             if (id(iconf) == 1) nsurvive=nsurvive+ &
                  cfg%csfs1E(ioff+1)-cfg%csfs1E(ioff)
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
    if (cfg%n2E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2E(n),cfg%off2E(n+1)-1

             ! Increment the configuration counter
             iconf=iconf+1

             ! Loop over the CSFs generated by this configuration
             do csf=cfg%csfs2E(ioff),cfg%csfs2E(ioff+1)-1

                ! Cycle if this is not a dominant CSF
                if (abs(vec(csf)) < cthrsh) cycle
                
                ! Flag the configuration as being dominant
                id(iconf)=1
                
             enddo

             ! Update the no. surviving CSFs
             if (id(iconf) == 1) nsurvive=nsurvive+ &
                  cfg%csfs2E(ioff+1)-cfg%csfs2E(ioff)
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
    if (cfg%n1I1E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 1I1E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1

             ! Increment the configuration counter
             iconf=iconf+1

             ! Loop over the CSFs generated by this configuration
             do csf=cfg%csfs1I1E(ioff),cfg%csfs1I1E(ioff+1)-1

                ! Cycle if this is not a dominant CSF
                if (abs(vec(csf)) < cthrsh) cycle
                
                ! Flag the configuration as being dominant
                id(iconf)=1
                
             enddo

             ! Update the no. surviving CSFs
             if (id(iconf) == 1) nsurvive=nsurvive+ &
                  cfg%csfs1I1E(ioff+1)-cfg%csfs1I1E(ioff)
             
          enddo

       enddo

    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(vec)
    
    return
    
  end subroutine fill_above_threshold
  
!######################################################################
! get_dominant_confs: returns the configurations with above-threshold
!                     coefficients for a single irrep
!######################################################################
  subroutine get_dominant_confs(cfg,dconf,ndconf,id,nconf,start)

    use constants
    use bitglobal
    use conftype
    
    implicit none
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in)   :: cfg

    ! Above-threshold configurations
    integer(is), intent(in)    :: ndconf
    integer(ib), intent(out)   :: dconf(n_int,2,ndconf)

    ! Total number of configurations
    integer(is), intent(in)    :: nconf

    ! Above-threshold configuration flags
    integer(is), intent(inout) :: id(nconf)

    ! Starting point in the id array
    integer(is), intent(in)    :: start

    ! Everything else
    integer(is)                :: iconf,n,ioff,counter
    integer(is)                :: n_int_I
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ! Configuration counter
    iconf=start-1

    ! Above threshold configuration counter
    counter=0
    
!----------------------------------------------------------------------
! Reference space configurations
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I

    ! Loop over ref space configurations
    do n=1,cfg%n0h

       ! Increment the configuration counter
       iconf=iconf+1

       ! Save the configuration if it is above threshold
       if (id(iconf) == 1) then
          counter=counter+1
          dconf(1:n_int_I,:,counter)=cfg%conf0h(:,:,n)
       endif
       
    enddo

!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
    if (cfg%n1I > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1I configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1I(n),cfg%off1I(n+1)-1

             ! Increment the configuration counter
             iconf=iconf+1

             ! Save the configuration if it is above threshold
             if (id(iconf) == 1) then
                counter=counter+1
                dconf(:,:,counter)=cfg%conf1I(:,:,ioff)
             endif
             
          enddo
          
       enddo

    endif
    
!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
    if (cfg%n2I > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2I configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2I(n),cfg%off2I(n+1)-1

             ! Increment the configuration counter
             iconf=iconf+1

             ! Save the configuration if it is above threshold
             if (id(iconf) == 1) then
                counter=counter+1
                dconf(:,:,counter)=cfg%conf2I(:,:,ioff)
             endif
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
    if (cfg%n1E > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1E configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1E(n),cfg%off1E(n+1)-1

             ! Increment the configuration counter
             iconf=iconf+1
             
             ! Save the configuration if it is above threshold
             if (id(iconf) == 1) then
                counter=counter+1
                dconf(:,:,counter)=cfg%conf1E(:,:,ioff)
             endif
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
    if (cfg%n2E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2E(n),cfg%off2E(n+1)-1

             ! Increment the configuration counter
             iconf=iconf+1

             ! Save the configuration if it is above threshold
             if (id(iconf) == 1) then
                counter=counter+1
                dconf(:,:,counter)=cfg%conf2E(:,:,ioff)
             endif
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
    if (cfg%n1I1E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 1I1E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1

             ! Increment the configuration counter
             iconf=iconf+1

             ! Save the configuration if it is above threshold
             if (id(iconf) == 1) then
                counter=counter+1
                dconf(:,:,counter)=cfg%conf1I1E(:,:,ioff)
             endif
             
          enddo

       enddo

    endif

    return
    
  end subroutine get_dominant_confs

!######################################################################
  
end module ref_refine
