!**********************************************************************
! Routines for the removal of dead wood from the reference space
!**********************************************************************

!######################################################################
! prune_ref_space: Removal of ref space configurations that do not
!                  contribute appreciably to the ref space eigenvectors
!######################################################################
#ifdef CBINDING
subroutine prune_ref_space(irrep,nroots,confscr,nconf,vecscr) &
     bind(c,name="prune_ref_space")
#else
subroutine prune_ref_space(irrep,nroots,confscr,nconf,vecscr)
#endif

  use bitglobal
  use constants
  use refconf
  use mrciutils
  use iomod

  ! Irrep number number of roots
  integer(is), intent(in)    :: irrep
  integer(is), intent(in)    :: nroots
  
  ! Array of reference configuration scratch file number
  integer(is), intent(in)    :: confscr

  ! Array of numbers of referecne configurations
  integer(is), intent(inout) :: nconf

  ! Eigenvector scratch file index
  integer(is), intent(in)    :: vecscr

  ! Reference space eigenpairs
  integer(is)                :: refdim
  real(dp), allocatable      :: e0(:),vec0(:,:)

  ! Number of configurations found in the scratch file
  integer(is)                :: nconf1

  ! Number of 64-bit integers required to represent the configurations
  integer(is)                :: n_int_I

  ! Internal and external MOs
  integer(is)                :: nmoI,nmoE
  integer(is)                :: Ilist(nmo),Elist(nmo)
  
  ! Reference configurations and SOPs
  integer(ib), allocatable   :: conf(:,:,:),sop(:,:,:)
  integer(ib), allocatable   :: conf_new(:,:,:),sop_new(:,:,:)
  
  ! MO mapping arrays
  integer(is)                :: m2c(nmo),c2m(nmo)

  ! CSF offsets
  integer(is), allocatable   :: offset(:)

  ! Deadwood flags
  integer(is), allocatable   :: idead(:)
  
  ! Everything else
  integer(is)                :: iscratch,i,n,iconf
  integer(is), allocatable   :: iroots(:)
  
!----------------------------------------------------------------------
! Determine the no. ref space CSFs
!----------------------------------------------------------------------
  ! Open the ref conf scratch file
  iscratch=scrunit(vecscr)
  open(iscratch,file=scrname(vecscr),form='unformatted',status='old')

  ! Read the no. CSFs
  read(iscratch) refdim

  ! Close the ref conf scratch file
  close(iscratch)
  
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  ! Reference space eigenvalues
  allocate(e0(nroots))
  e0=0.0d0

  ! Reference space eigenvectors
  allocate(vec0(refdim,nroots))
  vec0=0.0d0

  ! Required eigenvectors
  allocate(iroots(nroots))
  iroots=0
    
!----------------------------------------------------------------------
! Read in the reference space eigenvectors
!----------------------------------------------------------------------
! Note that we call read_some_eigenpairs here in case extra roots were
! computed in the ref space diagonalisation, as is done, e.g., when
! pruned MRCI is being used
!----------------------------------------------------------------------
  ! For now we will just use the first nroots eigenvectors
  do i=1,nroots
     iroots(i)=i
  enddo
  
  ! Read in the eigenpairs
  call read_some_eigenpairs(vecscr,vec0,e0,refdim,nroots,iroots)

!----------------------------------------------------------------------
! Read in the reference space configurations
!----------------------------------------------------------------------
  call read_ref_confs(confscr,nconf1,n_int_I,nmoI,nmoE,conf,sop,m2c,&
       c2m)

!----------------------------------------------------------------------
! Sanity check on the number of configurations
!----------------------------------------------------------------------
  if (nconf /= nconf1) then
     errmsg='Error in prune_ref_space: '&
          //'inconsistent configuration numbers'
     call error_control
  endif

!----------------------------------------------------------------------
! Get the CSF offsets
!----------------------------------------------------------------------
  allocate(offset(nconf1+1))
  offset=0

  call get_csf_offsets(nconf1,refdim,offset,sop,n_int_I)

!----------------------------------------------------------------------
! Determine the indices of the deadwood configurations
!----------------------------------------------------------------------
  allocate(idead(nconf1))
  idead=0
  
  call get_deadwood(nconf1,refdim,nroots,n_int_I,vec0,offset,idead)

!----------------------------------------------------------------------
! Construct the array of surviving configurations
!----------------------------------------------------------------------
  ! No. surviving confs
  nconf=nconf1-sum(idead)

  ! Allocate arrays
  allocate(conf_new(n_int,2,nconf))
  conf_new=0_ib
  
  ! Fill in the conf array
  n=0
  do iconf=1,nconf1
     if (idead(iconf) == 0) then
        n=n+1
        conf_new(1:n_int_I,:,n)=conf(:,:,iconf)
     endif
  enddo

!----------------------------------------------------------------------
! Put the surviving configurations into the canonical MO ordering
!----------------------------------------------------------------------
  call reorder_confs(m2c,conf_new,nconf)

!----------------------------------------------------------------------
! Update the internal-external MO spaces
!----------------------------------------------------------------------
  ! Get the lists of internal and external MOs
  call get_internal_external_mos(conf_new,nconf,nmoI,nmoE,Ilist,Elist)

  ! Construct the new canonical-to-MRCI MO index mapping array
  call get_mo_mapping(nmoI,nmoE,Ilist,Elist,m2c,c2m)

!----------------------------------------------------------------------
! Generate the surviving SOPs
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(sop_new(n_int,2,nconf))
  sop_new=0_ib

  ! Compute the SOPs
  do iconf=1,nconf
     sop_new(:,:,iconf)=conf_to_sop(conf_new(:,:,iconf),n_int)
  enddo
  
!----------------------------------------------------------------------
! Re-arrange the surviving configurations s.t. the internal MOs come
! before the external MOs
!----------------------------------------------------------------------
  call rearrange_ref_confs(c2m,conf_new,sop_new,nconf)

!----------------------------------------------------------------------
! Write the surviving configurations to file
!----------------------------------------------------------------------
  call rewrite_ref_confs(conf_new,sop_new,nconf,m2c,c2m,nmoI,nmoE,&
       confscr)
  
  return
  
end subroutine prune_ref_space

!######################################################################
! get_csf_offsets: determines the starting points for the CSFs
!                  generated by each configuration
!######################################################################
subroutine get_csf_offsets(nconf,refdim,offset,sop,n_int_I)

  use constants
  use bitglobal
  
  implicit none

  ! Dimensions
  integer(is), intent(in)  :: nconf,refdim,n_int_I

  ! CSF offsets
  integer(is), intent(out) :: offset(nconf+1)

  ! SOPs
  integer(ib), intent(in)  :: sop(n_int_I,2,nconf)

  ! Everything else
  integer(is), allocatable :: nopen(:)
  integer(is)              :: i,k,sum

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(nopen(nconf))
  
!----------------------------------------------------------------------
! Number of open shells for each configuration
!----------------------------------------------------------------------
  ! Initialisation
  nopen=0

  ! Loop over configurations
  do i=1,nconf

     ! Number of open shells
     do k=1,n_int_I
        nopen(i)=nopen(i)+popcnt(sop(k,1,i))
     enddo
     
  enddo

!----------------------------------------------------------------------
! Offsets
!----------------------------------------------------------------------
  sum=1
  offset=0

  do i=1,nconf
     offset(i)=sum
     sum=sum+ncsfs(nopen(i))
  enddo

  offset(nconf+1)=refdim+1
  
  return
  
end subroutine get_csf_offsets

!######################################################################
! get_deadwood: determines the indices of the deadwood configurations
!######################################################################
subroutine get_deadwood(nconf,refdim,nroots,n_int_I,vec0,offset,idead)

  use constants
  use bitglobal
  
  implicit none

  ! Dimensions
  integer(is), intent(in)  :: nconf,refdim,nroots,n_int_I

  ! Eigenvectors
  real(dp), intent(in)     :: vec0(refdim,nroots)

  ! CSF offsets
  integer(is), intent(in)  :: offset(nconf+1)

  ! Deadwood flags
  integer(is), intent(out) :: idead(nconf)

  ! Deadwood threshold (hard-coded for now)
  real(dp), parameter      :: cthrsh=0.005d0
  
  ! Everything else
  integer(is)              :: n,iconf,icsf

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
  idead=1
  
!----------------------------------------------------------------------
! Determine the deadwood configurations
!----------------------------------------------------------------------
  ! Loop over states
  do n=1,nroots

     ! Loop over configurations
     do iconf=1,nconf

        ! Loop over CSFs generated by this CSF
        do icsf=offset(iconf),offset(iconf+1)-1

           ! Above threshold coefficient?
           if (abs(vec0(icsf,n)) > cthrsh) then
              idead(iconf)=0
           endif
           
        enddo
        
     enddo
     
  enddo

  return
  
end subroutine get_deadwood

!######################################################################
! rewrite_ref_confs: re-writes a single ref conf scratch file with
!                    the non-deadwood configuration information
!######################################################################
subroutine rewrite_ref_confs(conf_new,sop_new,nconf,m2c,c2m,nmoI,&
     nmoE,confscr)

  use constants
  use bitglobal
  use iomod
  
  implicit none

  ! Dimensions
  integer(is), intent(in) :: nconf

  ! Configurations and SOPs
  integer(ib), intent(in) :: conf_new(n_int,2,nconf)
  integer(ib), intent(in) :: sop_new(n_int,2,nconf)

  ! No. internal and external MOs
  integer(is), intent(in) :: nmoI,nmoE

  ! Ref conf scratch file number
  integer(is), intent(in) :: confscr

  ! MO mapping arrays
  integer(is), intent(in) :: m2c(nmo),c2m(nmo)

  ! Everything else
  integer(is)             :: n_int_I
  integer(is)             :: iscratch

!----------------------------------------------------------------------
! Number of 64-bit integers needed to represent each reference space
! configuration and SOP bit string
!----------------------------------------------------------------------
  n_int_I=(nmoI-1)/64+1

!----------------------------------------------------------------------
! Overwrite the ref conf scratch file
!----------------------------------------------------------------------
  ! Open the scratch file
  iscratch=scrunit(confscr)
  open(iscratch,file=scrname(confscr),form='unformatted',&
       status='unknown')

  ! Number of reference space configurations
  write(iscratch) nconf
       
  ! Number of 64-bit integers needed to represent each reference
  ! space configuration and SOP bit string
  write(iscratch) n_int_I
  
  ! Number of internal and external MOs
  write(iscratch) nmoI
  write(iscratch) nmoE
  
  ! Configurations
  write(iscratch) conf_new(1:n_int_I,:,:)
       
  ! SOPs
  write(iscratch) sop_new(1:n_int_I,:,:)
  
  ! MO mapping arrays
  write(iscratch) m2c
  write(iscratch) c2m
  
  ! Close the scratch file
  close(iscratch)
  
  return
  
end subroutine rewrite_ref_confs
