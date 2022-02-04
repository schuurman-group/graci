!**********************************************************************
! Derived data type for the storage of the MRCI configurations for a
! single irrep. Uses the following tree-type data structure, similar
! to that described in Chem. Phys., 225, 197 (1997):
!
!         R            -> Reference confs
!         |
!        / \
!       /   \
!      /     \
!     1H      2H       -> 1-Hole and 2-Hole confs
!    /|      /|\
!   / |     / | \
! 1I 2I   2I  2E 1I1E  -> Confs with 1 or 2 created particles in the
!                         internal/external MO spaces
!**********************************************************************

module conftype

  use constants
  
  private

  type, public :: mrcfg
     
     ! Irrep number
     integer(is)              :: irrep

     ! Total numbers of confs and CSFs
     integer(is)              :: confdim,csfdim

     ! Reference space configurations of the
     ! symmetry 'irrep'
     integer(is)              :: n0h
     integer(is)              :: n_int_I,nmoI,nmoE
     integer(ib), allocatable :: conf0h(:,:,:)
     integer(ib), allocatable :: sop0h(:,:,:)

     ! Reference space configurations accross
     ! all symmetries
     integer(is)              :: nR
     integer(ib), allocatable :: confR(:,:,:)
     integer(ib), allocatable :: sopR(:,:,:)
     
     ! 1-hole and 2-hole configurations
     integer(is)              :: n1h,n2h
     integer(is), allocatable :: a1h(:),a2h(:,:)
     integer(is), allocatable :: off1h(:),off2h(:)
     integer(ib), allocatable :: conf1h(:,:,:),conf2h(:,:,:)
          
     ! 1E configurations
     integer(is)              :: n1E
     integer(is), allocatable :: a1E(:)
     integer(is), allocatable :: off1E(:)
     integer(ib), allocatable :: conf1E(:,:,:)
     integer(ib), allocatable :: sop1E(:,:,:)
     
     ! 2E configurations
     integer(is)              :: n2E
     integer(is), allocatable :: a2E(:,:)
     integer(is), allocatable :: off2E(:)
     integer(ib), allocatable :: conf2E(:,:,:)
     integer(ib), allocatable :: sop2E(:,:,:)
     
     ! 1I configurations
     integer(is)              :: n1I
     integer(is), allocatable :: a1I(:)
     integer(is), allocatable :: off1I(:)
     integer(ib), allocatable :: conf1I(:,:,:)
     integer(ib), allocatable :: sop1I(:,:,:)
     
     ! 2I configurations
     integer(is)              :: n2I
     integer(is), allocatable :: a2I(:,:)
     integer(is), allocatable :: off2I(:)
     integer(ib), allocatable :: conf2I(:,:,:)
     integer(ib), allocatable :: sop2I(:,:,:)
     
     ! 1I-1E configurations
     integer(is)              :: n1I1E
     integer(is), allocatable :: a1I1E(:,:)
     integer(is), allocatable :: off1I1E(:)
     integer(ib), allocatable :: conf1I1E(:,:,:)
     integer(ib), allocatable :: sop1I1E(:,:,:)

     ! Starting points of the CSFs generated by each configuration
     integer(is), allocatable :: csfs0h(:)
     integer(is), allocatable :: csfs1I(:)
     integer(is), allocatable :: csfs2I(:)
     integer(is), allocatable :: csfs1E(:)
     integer(is), allocatable :: csfs2E(:)
     integer(is), allocatable :: csfs1I1E(:)

     ! Complete sets of confs, SOPs, CSF offsets and
     ! creation operator indices
     ! Not filled in by default, but very useful for ref space
     ! merging
     integer(ib), allocatable :: confall(:,:,:)
     integer(ib), allocatable :: sopall(:,:,:)
     integer(is), allocatable :: csfsall(:)
     
     ! MO mapping arrays
     integer(is), allocatable :: m2c(:),c2m(:)

   contains

     ! Read in the configuration information from disk and initialise
     ! the derived data type
     procedure, non_overridable :: initialise
     
     ! Set the CFS offsets
     procedure, non_overridable :: set_csf_offsets

     ! Compute the conf and SOP bit strings
     procedure, non_overridable :: compute_confs_sops

     ! Fill in the complete arrays of confs, SOPs and CSF offsers
     procedure, non_overridable :: concatenate_arrays
     
     ! Finalisation/deallocation of arrays
     procedure, non_overridable :: finalise
     
  end type mrcfg

contains

!######################################################################
! initialise: Reads the configuration information from disk and fills
!             in the conf and SOP bit string arrays
!######################################################################
  subroutine initialise(cfg,irrep,scrnum)

    use constants
    use bitglobal

    implicit none

    class(mrcfg), intent(inout) :: cfg
    integer(is), intent(in)     :: irrep
    integer(is), intent(in)     :: scrnum
    integer(is)                 :: iscratch
    
!----------------------------------------------------------------------
! Set the irrep number
!----------------------------------------------------------------------
    cfg%irrep=irrep

!----------------------------------------------------------------------
! Read the configuration information from disk
!----------------------------------------------------------------------
    ! Open the scratch file
    iscratch=scrunit(scrnum)
    open(iscratch,file=scrname(scrnum),form='unformatted',status='old')
    
    ! Subspace dimensions
    read(iscratch) cfg%n_int_I
    read(iscratch) cfg%nmoI
    read(iscratch) cfg%nmoE
    read(iscratch) cfg%confdim
    read(iscratch) cfg%nR
    read(iscratch) cfg%n1h
    read(iscratch) cfg%n2h
    read(iscratch) cfg%n0h
    read(iscratch) cfg%n1I
    read(iscratch) cfg%n2I
    read(iscratch) cfg%n1E
    read(iscratch) cfg%n2E
    read(iscratch) cfg%n1I1E

    ! Allocate arrays
    allocate(cfg%confR(cfg%n_int_I,2,cfg%nR))
    allocate(cfg%conf0h(cfg%n_int_I,2,cfg%n0h))
    allocate(cfg%conf1h(cfg%n_int_I,2,cfg%n1h))
    allocate(cfg%a1h(cfg%n1h))
    allocate(cfg%off1h(cfg%nR+1))
    allocate(cfg%conf2h(cfg%n_int_I,2,cfg%n2h))
    allocate(cfg%a2h(2,cfg%n2h))
    allocate(cfg%off2h(cfg%nR+1))
    allocate(cfg%a1E(cfg%n1E))
    allocate(cfg%off1E(cfg%n1h+1))
    allocate(cfg%a2E(2,cfg%n2E))
    allocate(cfg%off2E(cfg%n2h+1))
    allocate(cfg%a1I(cfg%n1I))
    allocate(cfg%off1I(cfg%n1h+1))
    allocate(cfg%a2I(2,cfg%n2I))
    allocate(cfg%off2I(cfg%n2h+1))
    allocate(cfg%a1I1E(2,cfg%n1I1E))
    allocate(cfg%off1I1E(cfg%n2h+1))
    allocate(cfg%m2c(nmo))
    allocate(cfg%c2m(nmo))
    
    ! Configuration information
    read(iscratch) cfg%confR
    read(iscratch) cfg%conf1h
    read(iscratch) cfg%a1h
    read(iscratch) cfg%off1h
    read(iscratch) cfg%conf2h
    read(iscratch) cfg%a2h
    read(iscratch) cfg%off2h
    read(iscratch) cfg%conf0h
    if (cfg%n1I > 0) then
       read(iscratch) cfg%a1I
       read(iscratch) cfg%off1I
    endif
    if (cfg%n2I > 0) then
       read(iscratch) cfg%a2I
       read(iscratch) cfg%off2I
    endif
    if (cfg%n1E > 0) then
       read(iscratch) cfg%a1E
       read(iscratch) cfg%off1E
    endif
    if (cfg%n2E > 0) then
       read(iscratch) cfg%a2E
       read(iscratch) cfg%off2E
    endif
    if (cfg%n1I1E > 0) then
       read(iscratch) cfg%a1I1E
       read(iscratch) cfg%off1I1E
    endif
    
    ! MO mapping arrays
    read(iscratch) cfg%m2c
    read(iscratch) cfg%c2m
    
    ! Close the scratch file
    close(iscratch)

!----------------------------------------------------------------------
! Set up the CSF offsets
!----------------------------------------------------------------------
    call cfg%set_csf_offsets

!----------------------------------------------------------------------
! Compute up the conf and SOP bit strings
!----------------------------------------------------------------------
    call cfg%compute_confs_sops
    
    return
    
  end subroutine initialise

!######################################################################
! set_csf_offsets: Sets up the CSF offsets
!######################################################################
  subroutine set_csf_offsets(cfg)
    
    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    class(mrcfg), intent(inout) :: cfg

    ! Full configurations and SOPs
    integer(ib)                 :: conf(n_int,2)
    integer(ib)                 :: conf1(n_int,2)
    integer(ib)                 :: sop(n_int,2)
    
    ! Everything else
    integer(is)                 :: sum,nopen,counter,n_int_I
    integer(is)                 :: iconf,ioff,imo,imo1,imo2,n

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(cfg%csfs0h(cfg%n0h+1))
    allocate(cfg%csfs1I(cfg%n1I+1))
    allocate(cfg%csfs2I(cfg%n2I+1))
    allocate(cfg%csfs1E(cfg%n1E+1))
    allocate(cfg%csfs2E(cfg%n2E+1))
    allocate(cfg%csfs1I1E(cfg%n1I1E+1))
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I
    cfg%csfdim=0
    sum=1

!----------------------------------------------------------------------
! Reference space configurations
!----------------------------------------------------------------------
    ! Loop over reference configurations
    do iconf=1,cfg%n0h

       ! Construct the corresponding SOP
       conf=0_ib
       conf(1:n_int_I,:)=cfg%conf0h(:,:,iconf)
       sop=conf_to_sop(conf,n_int)

       ! Number of open shells
       nopen=sop_nopen(sop,n_int)

       ! Sum the number of CSFs generated by this configuration
       cfg%csfdim=cfg%csfdim+ncsfs(nopen)

       ! Offset
       cfg%csfs0h(iconf)=sum
       sum=sum+ncsfs(nopen)
       
    enddo

    ! Final offset
    cfg%csfs0h(cfg%n0h+1)=sum

!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
    if (cfg%n1I > 0) then

       ! Initialise the configuration counter
       counter=0

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
       
          ! Loop over the 1I configurations generated by the
          ! current 1-hole configuration
          do ioff=cfg%off1I(n),cfg%off1I(n+1)-1

             ! Increment the configuration counter
             counter=counter+1

             ! Construct the configuration
             imo=cfg%a1I(counter)
             conf1=0_ib
             conf1(1:n_int_I,:)=cfg%conf1h(:,:,n)
             conf=create_electron(conf1,n_int,imo)
             
             ! Construct the SOP
             sop=conf_to_sop(conf,n_int)
             
             ! Number of open shells
             nopen=sop_nopen(sop,n_int)

             ! Sum the number of CSFs generated by this configuration
             cfg%csfdim=cfg%csfdim+ncsfs(nopen)

             ! Offset
             cfg%csfs1I(counter)=sum
             sum=sum+ncsfs(nopen)

          enddo
          
       enddo

       ! Final offset
       cfg%csfs1I(cfg%n1I+1)=sum
       
    endif

!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
    if (cfg%n2I > 0) then
       
       ! Initialise the configuration counter
       counter=0

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
       
          ! Loop over the 2I configurations generated by the
          ! current 2-hole configuration
          do ioff=cfg%off2I(n),cfg%off2I(n+1)-1
          
             ! Increment the configuration counter
             counter=counter+1
          
             ! Construct the configuration
             imo1=cfg%a2I(1,counter)
             imo2=cfg%a2I(2,counter)
             conf1=0_ib
             conf1(1:n_int_I,:)=cfg%conf2h(:,:,n)
             conf=create_electron(conf1,n_int,imo1)
             conf1=conf
             conf=create_electron(conf1,n_int,imo2)
             
             ! Construct the SOP
             sop=conf_to_sop(conf,n_int)
             
             ! Number of open shells
             nopen=sop_nopen(sop,n_int)
          
             ! Sum the number of CSFs generated by this configuration
             cfg%csfdim=cfg%csfdim+ncsfs(nopen)

             ! Offset
             cfg%csfs2I(counter)=sum
             sum=sum+ncsfs(nopen)
             
          enddo
          
       enddo

       ! Final offset
       cfg%csfs2I(cfg%n2I+1)=sum

    endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
    if (cfg%n1E > 0) then

       ! Initialise the configuration counter
       counter=0

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1E configurations generated by the
          ! current 1-hole configuration
          do ioff=cfg%off1E(n),cfg%off1E(n+1)-1
             
             ! Increment the configuration counter
             counter=counter+1
             
             ! Construct the configuration
             imo=cfg%a1E(counter)
             conf1=0_ib
             conf1(1:n_int_I,:)=cfg%conf1h(:,:,n)
             conf=create_electron(conf1,n_int,imo)

             ! Construct the SOP
             sop=conf_to_sop(conf,n_int)
             
             ! Number of open shells
             nopen=sop_nopen(sop,n_int)
          
             ! Sum the number of CSFs generated by this configuration
             cfg%csfdim=cfg%csfdim+ncsfs(nopen)

             ! Offset
             cfg%csfs1E(counter)=sum
             sum=sum+ncsfs(nopen)
             
          enddo
          
       enddo

       ! Final offset
       cfg%csfs1E(cfg%n1E+1)=sum
       
    endif

!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
    if (cfg%n2E > 0) then
    
       ! Initialise the configuration counter
       counter=0

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
     
          ! Loop over the 2E configurations generated by the
          ! current 2-hole configuration
          do ioff=cfg%off2E(n),cfg%off2E(n+1)-1
             
             ! Increment the configuration counter
             counter=counter+1

             ! Construct the configuration
             imo1=cfg%a2E(1,counter)
             imo2=cfg%a2E(2,counter)
             conf1=0_ib
             conf1(1:n_int_I,:)=cfg%conf2h(:,:,n)
             conf=create_electron(conf1,n_int,imo1)
             conf1=conf
             conf=create_electron(conf1,n_int,imo2)

             ! Construct the SOP
             sop=conf_to_sop(conf,n_int)
             
             ! Number of open shells
             nopen=sop_nopen(sop,n_int)
          
             ! Sum the number of CSFs generated by this configuration
             cfg%csfdim=cfg%csfdim+ncsfs(nopen)

             ! Offset
             cfg%csfs2E(counter)=sum
             sum=sum+ncsfs(nopen)
             
          enddo

       enddo

       ! Final offset
       cfg%csfs2E(cfg%n2E+1)=sum
       
    endif

!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
    if (cfg%n1I1E > 0) then
    
       ! Initialise the configuration counter
       counter=0

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
     
          ! Loop over the 1I1E configurations generated by
          ! the current 2-hole configuration
          do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1
        
             ! Increment the configuration counter
             counter=counter+1

             ! Construct the configuration
             imo1=cfg%a1I1E(1,counter)
             imo2=cfg%a1I1E(2,counter)
             conf1=0_ib
             conf1(1:n_int_I,:)=cfg%conf2h(:,:,n)
             conf=create_electron(conf1,n_int,imo1)
             conf1=conf
             conf=create_electron(conf1,n_int,imo2)

             ! Construct the SOP
             sop=conf_to_sop(conf,n_int)
             
             ! Number of open shells
             nopen=sop_nopen(sop,n_int)
          
             ! Sum the number of CSFs generated by this configuration
             cfg%csfdim=cfg%csfdim+ncsfs(nopen)

             ! Offset
             cfg%csfs1I1E(counter)=sum
             sum=sum+ncsfs(nopen)
             
          enddo
          
       enddo

       ! Final offset
       cfg%csfs1I1E(cfg%n1I1E+1)=sum

    endif
    
    return
    
  end subroutine set_csf_offsets

!######################################################################
! compute_confs_sops: Computes the conf and SOP bit strings
!######################################################################
  subroutine compute_confs_sops(cfg)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    class(mrcfg), intent(inout) :: cfg

    ! Working arrays
    integer(ib)                 :: conf1(n_int,2)
    
    ! Everything else
    integer(is)                 :: n_int_I
    integer(is)                 :: iconf,n,ioff,imo,imo1,imo2,counter
    real(dp)                    :: mem
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I
    
    allocate(cfg%sop0h(n_int_I,2,cfg%n0h))
    allocate(cfg%sopR(n_int_I,2,cfg%nR))
    allocate(cfg%sop1I(n_int,2,cfg%n1I))
    allocate(cfg%sop2I(n_int,2,cfg%n2I))
    allocate(cfg%sop1E(n_int,2,cfg%n1E))
    allocate(cfg%sop2E(n_int,2,cfg%n2E))
    allocate(cfg%sop1I1E(n_int,2,cfg%n1I1E))
    allocate(cfg%conf1I(n_int,2,cfg%n1I))
    allocate(cfg%conf2I(n_int,2,cfg%n2I))
    allocate(cfg%conf1E(n_int,2,cfg%n1E))
    allocate(cfg%conf2E(n_int,2,cfg%n2E))
    allocate(cfg%conf1I1E(n_int,2,cfg%n1I1E))

    cfg%conf1I=0_ib; cfg%conf2I=0_ib; cfg%conf1E=0_ib
    cfg%conf2E=0_ib; cfg%conf1I1E=0_ib
    cfg%sopR=0_ib; cfg%sop0h=0_ib; cfg%sop1I=0_ib; cfg%sop2I=0_ib;
    cfg%sop1E=0_ib; cfg%sop2E=0_ib; cfg%sop1I1E=0_ib

!----------------------------------------------------------------------
! Reference space SOPs across all irreps
!----------------------------------------------------------------------
    ! Loop over reference configurations
    do iconf=1,cfg%nR
       
       ! Construct the SOP
       cfg%sopR(:,:,iconf)=conf_to_sop(cfg%confR(:,:,iconf),n_int_I)

    enddo
    
!----------------------------------------------------------------------
! Reference space SOPs for the single irrep under consideration
!----------------------------------------------------------------------
    ! Loop over reference configurations
    do iconf=1,cfg%n0h
       
       ! Construct the SOP
       cfg%sop0h(:,:,iconf)=conf_to_sop(cfg%conf0h(:,:,iconf),n_int_I)

    enddo

!----------------------------------------------------------------------    
! 1I configurations
!----------------------------------------------------------------------
    if (cfg%n1I > 0) then
       
       ! Initialise the configuration counter
       counter=0

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1I configurations generated by the
          ! current 1-hole configuration
          do ioff=cfg%off1I(n),cfg%off1I(n+1)-1

             ! Increment the configuration counter
             counter=counter+1

             ! Construct the configuration
             imo=cfg%a1I(counter)
             conf1=0_ib
             conf1(1:n_int_I,:)=cfg%conf1h(:,:,n)
             cfg%conf1I(:,:,counter)=create_electron(conf1,n_int,imo)
             
             ! Construct the SOP
             cfg%sop1I(:,:,counter)=&
                  conf_to_sop(cfg%conf1I(:,:,counter),n_int)
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
    if (cfg%n2I > 0) then

       ! Initialise the configuration counter
       counter=0

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
       
          ! Loop over the 2I configurations generated by the
          ! current 2-hole configuration
          do ioff=cfg%off2I(n),cfg%off2I(n+1)-1
          
             ! Increment the configuration counter
             counter=counter+1
          
             ! Construct the configuration
             imo1=cfg%a2I(1,counter)
             imo2=cfg%a2I(2,counter)
             conf1=0_ib
             conf1(1:n_int_I,:)=cfg%conf2h(:,:,n)
             cfg%conf2I(:,:,counter)=create_electron(conf1,n_int,imo1)
             conf1=cfg%conf2I(:,:,counter)
             cfg%conf2I(:,:,counter)=create_electron(conf1,n_int,imo2)
             
             ! Construct the SOP
             cfg%sop2I(:,:,counter)=&
                  conf_to_sop(cfg%conf2I(:,:,counter),n_int)
             
          enddo
          
       enddo
       
    endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
    if (cfg%n1E > 0) then

       ! Initialise the configuration counter
       counter=0

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1E configurations generated by the
          ! current 1-hole configuration
          do ioff=cfg%off1E(n),cfg%off1E(n+1)-1
             
             ! Increment the configuration counter
             counter=counter+1
             
             ! Construct the configuration
             imo=cfg%a1E(counter)
             conf1=0_ib
             conf1(1:n_int_I,:)=cfg%conf1h(:,:,n)
             cfg%conf1E(:,:,counter)=create_electron(conf1,n_int,imo)

             ! Construct the SOP
             cfg%sop1E(:,:,counter)=&
                  conf_to_sop(cfg%conf1E(:,:,counter),n_int)
             
          enddo
          
       enddo

    endif

!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
    if (cfg%n2E > 0) then
    
       ! Initialise the configuration counter
       counter=0

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
     
          ! Loop over the 2E configurations generated by the
          ! current 2-hole configuration
          do ioff=cfg%off2E(n),cfg%off2E(n+1)-1
             
             ! Increment the configuration counter
             counter=counter+1

             ! Construct the configuration
             imo1=cfg%a2E(1,counter)
             imo2=cfg%a2E(2,counter)
             conf1=0_ib
             conf1(1:n_int_I,:)=cfg%conf2h(:,:,n)
             cfg%conf2E(:,:,counter)=create_electron(conf1,n_int,imo1)
             conf1=cfg%conf2E(:,:,counter)
             cfg%conf2E(:,:,counter)=create_electron(conf1,n_int,imo2)

             ! Construct the SOP
             cfg%sop2E(:,:,counter)=&
                  conf_to_sop(cfg%conf2E(:,:,counter),n_int)
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
    if (cfg%n1I1E > 0) then
    
       ! Initialise the configuration counter
       counter=0

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
     
          ! Loop over the 1I1E configurations generated by
          ! the current 2-hole configuration
          do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1
        
             ! Increment the configuration counter
             counter=counter+1

             ! Construct the configuration
             imo1=cfg%a1I1E(1,counter)
             imo2=cfg%a1I1E(2,counter)
             conf1=0_ib
             conf1(1:n_int_I,:)=cfg%conf2h(:,:,n)
             cfg%conf1I1E(:,:,counter)=create_electron(conf1,n_int,imo1)
             conf1=cfg%conf1I1E(:,:,counter)
             cfg%conf1I1E(:,:,counter)=create_electron(conf1,n_int,imo2)

             ! Construct the SOP
             cfg%sop1I1E(:,:,counter)=&
                  conf_to_sop(cfg%conf1I1E(:,:,counter),n_int)
             
          enddo
          
       enddo

    endif
    
    return
    
  end subroutine compute_confs_sops

!######################################################################
! concatenate_arrays: Concatenates all the confs, SOPs, CSF offsets
!                     and creation operator indices into the confall,
!                     sopall and csfsall
!######################################################################
  subroutine concatenate_arrays(cfg)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    class(mrcfg), intent(inout) :: cfg

    ! Everything else
    integer(is)              :: n,iconf,n_int_I,counter
    integer(is)              :: sum,nopen

!----------------------------------------------------------------------
! Allocation and initialisation
!----------------------------------------------------------------------
    ! Confs
    if (allocated(cfg%confall)) deallocate(cfg%confall)
    allocate(cfg%confall(n_int,2,cfg%confdim))
    cfg%confall=0_ib

    ! SOPs
    if (allocated(cfg%sopall)) deallocate(cfg%sopall)
    allocate(cfg%sopall(n_int,2,cfg%confdim))
    cfg%sopall=0_ib

    ! CSF offsets
    if (allocated(cfg%csfsall)) deallocate(cfg%csfsall)
    allocate(cfg%csfsall(cfg%confdim+1))
    cfg%csfsall=0

    ! Conf counter
    counter=0

    ! CSF sum
    sum=1

!----------------------------------------------------------------------
! Reference configurations
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I

    ! Loop over reference configurations
    do iconf=1,cfg%n0h

       ! Increment the counter
       counter=counter+1

       ! Save the conf and SOP
       cfg%confall(1:n_int_I,:,counter)=cfg%conf0h(:,:,iconf)
       cfg%sopall(1:n_int_I,:,counter)=cfg%sop0h(:,:,iconf)
     
       ! Number of open shells
       nopen=sop_nopen(cfg%sopall(:,:,counter),n_int)
       
       ! CSF offset
       cfg%csfsall(counter)=sum
       sum=sum+ncsfs(nopen)

    enddo

!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
    if (cfg%n1I > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
        
          ! Loop over the 1I configurations generated by the 1-hole
          ! configuration
          do iconf=cfg%off1I(n),cfg%off1I(n+1)-1

             ! Increment the counter
             counter=counter+1
           
             ! Save the conf and SOP
             cfg%confall(:,:,counter)=cfg%conf1I(:,:,iconf)
             cfg%sopall(:,:,counter)=cfg%sop1I(:,:,iconf)
           
             ! Number of open shells
             nopen=sop_nopen(cfg%sopall(:,:,counter),n_int)
     
             ! CSF offset
             cfg%csfsall(counter)=sum
             sum=sum+ncsfs(nopen)

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
          do iconf=cfg%off2I(n),cfg%off2I(n+1)-1

             ! Increment the counter
             counter=counter+1
           
             ! Save the conf and SOP
             cfg%confall(:,:,counter)=cfg%conf2I(:,:,iconf)
             cfg%sopall(:,:,counter)=cfg%sop2I(:,:,iconf)
             
             ! Number of open shells
             nopen=sop_nopen(cfg%sopall(:,:,counter),n_int)
           
             ! CSF offset
             cfg%csfsall(counter)=sum
             sum=sum+ncsfs(nopen)

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
          do iconf=cfg%off1E(n),cfg%off1E(n+1)-1

             ! Increment the counter
             counter=counter+1
           
             ! Save the conf and SOP
             cfg%confall(:,:,counter)=cfg%conf1E(:,:,iconf)
             cfg%sopall(:,:,counter)=cfg%sop1E(:,:,iconf)
           
             ! Number of open shells
             nopen=sop_nopen(cfg%sopall(:,:,counter),n_int)
             
             ! CSF offset
             cfg%csfsall(counter)=sum
             sum=sum+ncsfs(nopen)

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
          do iconf=cfg%off2E(n),cfg%off2E(n+1)-1

             ! Increment the counter
             counter=counter+1
           
             ! Save the conf and SOP
             cfg%confall(:,:,counter)=cfg%conf2E(:,:,iconf)
             cfg%sopall(:,:,counter)=cfg%sop2E(:,:,iconf)
             
             ! Number of open shells
             nopen=sop_nopen(cfg%sopall(:,:,counter),n_int)
             
             ! CSF offset
             cfg%csfsall(counter)=sum
             sum=sum+ncsfs(nopen)

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
          do iconf=cfg%off1I1E(n),cfg%off1I1E(n+1)-1
             
             ! Increment the counter
             counter=counter+1
             
             ! Save the conf and SOP
             cfg%confall(:,:,counter)=cfg%conf1I1E(:,:,iconf)
             cfg%sopall(:,:,counter)=cfg%sop1I1E(:,:,iconf)
           
             ! Number of open shells
             nopen=sop_nopen(cfg%sopall(:,:,counter),n_int)
           
             ! CSF offset
             cfg%csfsall(counter)=sum
             sum=sum+ncsfs(nopen)

          enddo
           
       enddo

    endif

!----------------------------------------------------------------------
! Final CSF offset element
!----------------------------------------------------------------------
    cfg%csfsall(cfg%confdim+1)=sum
    
    return
    
  end subroutine concatenate_arrays
    
!######################################################################
! finalise: Deallocates all arrays and zeros all scalar variables
!######################################################################
  subroutine finalise(cfg)

    use constants

    implicit none

    class(mrcfg), intent(inout) :: cfg

    !
    ! Zero all scalar variables
    !
    cfg%irrep=0
    cfg%n_int_I=0
    cfg%nmoI=0
    cfg%nmoE=0
    cfg%confdim=0
    cfg%nR=0
    cfg%n1h=0
    cfg%n2h=0
    cfg%n0h=0
    cfg%n1I=0
    cfg%n2I=0
    cfg%n1E=0
    cfg%n2E=0
    cfg%n1I1E=0

    !
    ! Deallocate all arrays
    !
    if (allocated(cfg%confR)) deallocate(cfg%confR)
    if (allocated(cfg%a1h)) deallocate(cfg%a1h)
    if (allocated(cfg%off1h)) deallocate(cfg%off1h)
    if (allocated(cfg%a2h)) deallocate(cfg%a2h)
    if (allocated(cfg%off2h)) deallocate(cfg%off2h)
    if (allocated(cfg%conf0h)) deallocate(cfg%conf0h)
    if (allocated(cfg%conf1h)) deallocate(cfg%conf1h)
    if (allocated(cfg%conf2h)) deallocate(cfg%conf2h)
    if (allocated(cfg%a1E)) deallocate(cfg%a1E)
    if (allocated(cfg%off1E)) deallocate(cfg%off1E)
    if (allocated(cfg%a2E)) deallocate(cfg%a2E)
    if (allocated(cfg%off2E)) deallocate(cfg%off2E)
    if (allocated(cfg%a1I)) deallocate(cfg%a1I)
    if (allocated(cfg%off1I)) deallocate(cfg%off1I)
    if (allocated(cfg%a2I)) deallocate(cfg%a2I)
    if (allocated(cfg%off2I)) deallocate(cfg%off2I)
    if (allocated(cfg%a1I1E)) deallocate(cfg%a1I1E)
    if (allocated(cfg%off1I1E)) deallocate(cfg%off1I1E)
    if (allocated(cfg%m2c)) deallocate(cfg%m2c)
    if (allocated(cfg%c2m)) deallocate(cfg%c2m)
    if (allocated(cfg%csfs0h)) deallocate(cfg%csfs0h)
    if (allocated(cfg%csfs1I)) deallocate(cfg%csfs1I)
    if (allocated(cfg%csfs2I)) deallocate(cfg%csfs2I)
    if (allocated(cfg%csfs1E)) deallocate(cfg%csfs1E)
    if (allocated(cfg%csfs2E)) deallocate(cfg%csfs2E)
    if (allocated(cfg%csfs1I1E)) deallocate(cfg%csfs1I1E)
    if (allocated(cfg%sopR)) deallocate(cfg%sopR)
    if (allocated(cfg%sop0h)) deallocate(cfg%sop0h)
    if (allocated(cfg%sop1I)) deallocate(cfg%sop1I)
    if (allocated(cfg%sop2I)) deallocate(cfg%sop2I)
    if (allocated(cfg%sop1E)) deallocate(cfg%sop1E)
    if (allocated(cfg%sop2E)) deallocate(cfg%sop2E)
    if (allocated(cfg%sop1I1E)) deallocate(cfg%sop1I1E)
    if (allocated(cfg%conf1I)) deallocate(cfg%conf1I)
    if (allocated(cfg%conf2I)) deallocate(cfg%conf2I)
    if (allocated(cfg%conf1E)) deallocate(cfg%conf1E)
    if (allocated(cfg%conf2E)) deallocate(cfg%conf2E)
    if (allocated(cfg%conf1I1E)) deallocate(cfg%conf1I1E)
    if (allocated(cfg%confall)) deallocate(cfg%confall)
    if (allocated(cfg%sopall)) deallocate(cfg%sopall)
    if (allocated(cfg%csfsall)) deallocate(cfg%csfsall)
    
    return
    
  end subroutine finalise
    
!######################################################################
  
end module conftype
