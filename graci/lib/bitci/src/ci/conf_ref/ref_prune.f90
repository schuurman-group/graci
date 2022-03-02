!**********************************************************************
! Routines for the removal of dead wood from the reference space
!**********************************************************************

!######################################################################
! prune_ref_space: Removal of ref space configurations that do not
!                  contribute appreciably to the ref space eigenvectors
!######################################################################
#ifdef CBINDING
subroutine prune_ref_space(nroots,confscr,nconf,vecscr) &
     bind(c,name="prune_ref_space")
#else
subroutine prune_ref_space(nroots,confscr,nconf,vecscr)
#endif
  
  use bitglobal
  use constants
  use refconf
  use mrciutils
  use iomod

  implicit none
  
  ! Number of roots per irrep
  integer(is), intent(in)    :: nroots(0:nirrep-1)

  ! Array of reference configuration scratch file numbers
  integer(is), intent(in)    :: confscr(0:nirrep-1)

  ! Array of numbers of referecne configurations
  integer(is), intent(inout) :: nconf(0:nirrep-1)

  ! Eigenvector scratch file indices
  integer(is), intent(in)    :: vecscr(0:nirrep-1)

  ! Reference space eigenpairs
  integer(is)                :: refdim(0:nirrep-1)
  real(dp), allocatable      :: e0(:),vec0(:,:)

  ! Number of configurations found in the scratch file
  integer(is)                :: nconf1

  ! Number of 64-bit integers required to represent the internal MOs
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
  integer(is)                :: iscratch,i,n,iconf,irrep
  integer(is), allocatable   :: iroots(:)
  integer(is)                :: nroots_tot
  integer(is)                :: refdim_tot
  integer(is)                :: nconf_tot,nconf_tot_old
  integer(ib), allocatable   :: conf1(:,:,:),sop1(:,:,:)
  integer(is)                :: istart,iend
  
!----------------------------------------------------------------------
! Determine the no. ref space CSFs per irrep
!----------------------------------------------------------------------
  ! Loop over irreps
  do irrep=0,nirrep-1
  
     ! Open the ref conf scratch file
     iscratch=scrunit(vecscr(irrep))
     open(iscratch,file=scrname(vecscr(irrep)),form='unformatted',&
          status='old')

     ! No. CSFs
     read(iscratch) refdim(irrep)
     
     ! Close the ref conf scratch file
     close(iscratch)

  enddo
  
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  ! Total number of confs
  nconf_tot=sum(nconf)

  ! Total no. CSFs
  refdim_tot=sum(refdim)

  ! Total no. roots
  nroots_tot=sum(nroots)

  ! Reference space confs and SOPs across all irreps
  allocate(conf(n_int,2,nconf_tot))
  allocate(sop(n_int,2,nconf_tot))
  conf=0_ib
  sop=0_ib

  ! Deadwood conf indices
  allocate(idead(nconf_tot))
  idead=0
  
!----------------------------------------------------------------------
! Read in the reference space configurations
!----------------------------------------------------------------------
! Note that n_int_I, nmoI, nmoE, m2c, and c2m will be the same
! for all irreps
!----------------------------------------------------------------------
  ! Initialisation
  istart=1
  iend=0
  
  ! Loop over irreps
  do irrep=0,nirrep-1

     ! Start and end points in the total conf & SOP arrays
     if (irrep > 0) istart=istart+nconf(irrep-1)
     iend=iend+nconf(irrep)

     ! Read in the configurations
     call read_ref_confs(confscr(irrep),nconf1,n_int_I,nmoI,nmoE,&
          conf1,sop1,m2c,c2m)

     ! Sanity check on the number of configurations
     if (nconf(irrep) /= nconf1) then
        errmsg='Error in prune_ref_space: '&
             //'inconsistent configuration numbers'
        call error_control
     endif

     ! Save the confs & SOPs
     conf(1:n_int_I,:,istart:iend)=conf1
     sop(1:n_int_I,:,istart:iend)=sop1
     
     ! Deallocate arrays
     deallocate(conf1,sop1)
     
  enddo

!----------------------------------------------------------------------
! Determine the deadwood configurations
!----------------------------------------------------------------------
! Note that we call read_some_eigenpairs here in case extra roots were
! computed in the ref space diagonalisation, as is done, e.g., when
! pruned MRCI is being used
!----------------------------------------------------------------------
  ! Initialisation
  istart=1
  iend=0

  ! Loop over irreps
  do irrep=0,nirrep-1
     
     ! Allocate arrays
     allocate(iroots(nroots(irrep)))
     iroots=0
     allocate(vec0(refdim(irrep),nroots(irrep)))
     vec0=0.0d0
     allocate(e0(nroots(irrep)))
     e0=0.0d0
     allocate(offset(nconf(irrep)+1))
     offset=0
     
     ! Start and end points in the total conf & SOP arrays
     if (irrep > 0) istart=istart+nconf(irrep-1)
     iend=iend+nconf(irrep)
     
     ! Indices of the roots to be read in
     do i=1,nroots(irrep)
        iroots(i)=i
     enddo
     
     ! Read in the reference space eigenvectors for this irrep
     call read_some_eigenpairs(vecscr(irrep),vec0,e0,refdim(irrep),&
          nroots(irrep),iroots)

     ! Get the CSF offsets
     call get_csf_offsets(nconf(irrep),refdim(irrep),offset,&
          sop(:,:,istart:iend))

     ! Determine the indices of the deadwood configurations for
     ! this irrep
     call get_deadwood(nconf(irrep),refdim(irrep),nroots(irrep),&
          vec0,offset,idead(istart:iend))

     !print*,irrep,nconf(irrep)-sum(idead(istart:iend))
     
     ! Deallocate arrays
     deallocate(iroots,vec0,e0,offset)
     
  enddo
  
!----------------------------------------------------------------------
! Construct the array of surviving configurations
!----------------------------------------------------------------------
  ! No. surviving confs
  nconf_tot_old=nconf_tot
  nconf_tot=nconf_tot_old-sum(idead)

  ! Allocate arrays
  allocate(conf_new(n_int,2,nconf_tot))
  conf_new=0_ib
  
  ! Fill in the new conf array
  n=0
  do iconf=1,nconf_tot_old
     if (idead(iconf) == 0) then
        n=n+1
        conf_new(:,:,n)=conf(:,:,iconf)
     endif
  enddo
  
!----------------------------------------------------------------------
! Put the surviving configurations into the canonical MO ordering
!----------------------------------------------------------------------
  call reorder_confs(m2c,conf_new,nconf_tot)
  
!----------------------------------------------------------------------
! Update the internal-external MO spaces
!----------------------------------------------------------------------
  ! Get the lists of internal and external MOs
  call get_internal_external_mos(conf_new,nconf_tot,nmoI,nmoE,Ilist,&
       Elist)
  
  ! Construct the new canonical-to-MRCI MO index mapping array
  call get_mo_mapping(nmoI,nmoE,Ilist,Elist,m2c,c2m)
  
!----------------------------------------------------------------------
! Generate the surviving SOPs in the canonical MO ordering
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(sop_new(n_int,2,nconf_tot))
  sop_new=0_ib

  ! Compute the SOPs
  do iconf=1,nconf_tot
     sop_new(:,:,iconf)=conf_to_sop(conf_new(:,:,iconf),n_int)
  enddo

!----------------------------------------------------------------------
! Re-arrange the surviving configurations s.t. the internal MOs come
! before the external MOs
!----------------------------------------------------------------------
  call rearrange_ref_confs(c2m,conf_new,sop_new,nconf_tot)

!----------------------------------------------------------------------
! Determine the new number of configurations per irrep
!----------------------------------------------------------------------
  nconf=0
  
  do iconf=1,nconf_tot
     irrep=sop_sym_mrci(sop_new(:,:,iconf),m2c)
     nconf(irrep)=nconf(irrep)+1
  enddo

!----------------------------------------------------------------------
! Write the surviving configurations to file
!----------------------------------------------------------------------
  call rewrite_ref_confs(conf_new,sop_new,nconf_tot,nconf,m2c,c2m,&
       nmoI,nmoE,confscr)

  return
  
end subroutine prune_ref_space

!######################################################################
! get_csf_offsets: determines the starting points for the CSFs
!                  generated by each configuration
!######################################################################
subroutine get_csf_offsets(nconf,refdim,offset,sop)

  use constants
  use bitglobal
  
  implicit none

  ! Dimensions
  integer(is), intent(in)  :: nconf,refdim

  ! CSF offsets
  integer(is), intent(out) :: offset(nconf+1)

  ! SOPs
  integer(ib), intent(in)  :: sop(n_int,2,nconf)

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
     do k=1,n_int
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
subroutine get_deadwood(nconf,refdim,nroots,vec0,offset,idead)

  use constants
  use bitglobal
  use utils
  
  implicit none

  ! Dimensions
  integer(is), intent(in)  :: nconf,refdim,nroots

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
subroutine rewrite_ref_confs(conf_new,sop_new,nconf_tot,nconf,&
     m2c,c2m,nmoI,nmoE,confscr)

  use constants
  use bitglobal
  use iomod
  
  implicit none

  ! Dimensions
  integer(is), intent(in) :: nconf(0:nirrep-1)
  integer(is), intent(in) :: nconf_tot

  ! Configurations and SOPs
  integer(ib), intent(in) :: conf_new(n_int,2,nconf_tot)
  integer(ib), intent(in) :: sop_new(n_int,2,nconf_tot)

  ! No. internal and external MOs
  integer(is), intent(in) :: nmoI,nmoE

  ! Ref conf scratch file numbers
  integer(is), intent(in) :: confscr(0:nirrep-1)

  ! MO mapping arrays
  integer(is), intent(in) :: m2c(nmo),c2m(nmo)

  ! Everything else
  integer(is)             :: n_int_I
  integer(is)             :: irrep,istart,iend
  integer(is)             :: iscratch

  
!----------------------------------------------------------------------
! Number of 64-bit integers needed to represent each reference space
! configuration and SOP bit string
!----------------------------------------------------------------------
  n_int_I=(nmoI-1)/64+1

!----------------------------------------------------------------------
! Overwrite the ref conf scratch files
!----------------------------------------------------------------------
  istart=1
  iend=0

  ! Loop over irreps
  do irrep=0,nirrep-1

     ! Start and end points in the total conf & SOP arrays
     if (irrep > 0) istart=istart+nconf(irrep-1)
     iend=iend+nconf(irrep)

     ! Open the scratch file
     iscratch=scrunit(confscr(irrep))
     open(iscratch,file=scrname(confscr(irrep)),form='unformatted',&
          status='unknown')
  
     ! Number of reference space configurations
     write(iscratch) nconf(irrep)
     
     ! Number of 64-bit integers needed to represent each reference
     ! space configuration and SOP bit string
     write(iscratch) n_int_I
  
     ! Number of internal and external MOs
     write(iscratch) nmoI
     write(iscratch) nmoE
  
     ! Configurations
     write(iscratch) conf_new(1:n_int_I,:,istart:iend)
       
     ! SOPs
     write(iscratch) sop_new(1:n_int_I,:,istart:iend)
  
     ! MO mapping arrays
     write(iscratch) m2c
     write(iscratch) c2m
  
     ! Close the scratch file
     close(iscratch)

  enddo
     
  return
  
end subroutine rewrite_ref_confs
