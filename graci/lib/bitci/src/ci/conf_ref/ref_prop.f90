!**********************************************************************
! Routines for the generation of reference space configurations via
! the approximate re-expression of those of another geometry in terms
! of the MOs of the current one
!**********************************************************************

!######################################################################
! 
!######################################################################
#ifdef CBINDING
subroutine ref_space_propagate(nmo0,nmo1,smat,conffile0_in,nconf,&
     confscr) bind(c,name="ref_space_propagate")
#else
subroutine ref_space_propagate(nmo0,nmo1,smat,conffile0_in,nconf,&
     confscr)
#endif
  
  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use dethash
  use mrciutils
  use iomod

  ! MO basis dimensions
  integer(is), intent(in)            :: nmo0,nmo1

  ! MO overlaps
  real(dp), intent(in)               :: smat(nmo0,nmo1)

  ! Scratch file name for the ref space configurations to propagate
  ! forwards
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: conffile0_in
  character(len=255)                 :: conffile0
  integer(is)                        :: length
#else
  character(len=*)                   :: conffile0_in
  character(len=255)                 :: conffile0
#endif

  ! Number of generated reference configurations
  integer(is), intent(out)           :: nconf

  ! Generated reference configuration scratch file number
  integer(is), intent(out)           :: confscr

  ! Reference confs, SOPs, etc. at the previous geometry
  integer(is)                        :: scrnum0
  integer(is)                        :: nconf0,n_int_I,nmoI,nmoE
  integer(ib), allocatable           :: conf_r0(:,:,:),sop_r0(:,:,:)
  integer(ib), allocatable           :: conf(:,:,:),sop(:,:,:)     
  integer(is), allocatable           :: m2c0(:),c2m0(:)

  ! Hash table
  type(dhtbl)                        :: h
  integer(is)                        :: initial_size
  integer(ib)                        :: key(n_int,2)
  
  ! Everything else
  integer(is)                        :: i,j,k,n,imo0,imo
  
!----------------------------------------------------------------------
! If C bindings are on, then convert the previous geometry
! ref conf file name from the C char type to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(conffile0_in)
  call c2fstr(conffile0_in,conffile0,length)
#else
  conffile0=adjustl(trim(conffile0_in))
#endif

!----------------------------------------------------------------------
! Register the previous geometry ref conf scratch file
!----------------------------------------------------------------------
  call register_scratch_file(scrnum0,conffile0)

!----------------------------------------------------------------------
! Read in the previous geometry ref confs
!----------------------------------------------------------------------
  allocate(m2c0(nmo0),c2m0(nmo0))

  call read_ref_confs(scrnum0,nconf0,n_int_I,nmoI,nmoE,conf_r0,&
       sop_r0,m2c0,c2m0)

!----------------------------------------------------------------------
! Put the ref confs into the 'canonical' MO order
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(conf(n_int,2,nconf0))
  conf=0_ib

  ! Resize the conf array
  conf(1:n_int_I,:,:)=conf_r0

  ! Reorder the conf array
  call reorder_confs(m2c0,conf,nconf0)
  
!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------
  ! Initialisation
  initial_size=nconf0
  call h%initialise_table(initial_size)
  
!----------------------------------------------------------------------
! Construct the current geometry ref confs
!----------------------------------------------------------------------
  ! Loop over previous geometry ref confs
  do i=1,nconf0

     ! Set the new conf assuming a one-to-one MO mapping

     print*,''

     key=0_ib
     do j=1,2
        do k=1,n_int
           do n=0,n_bits-1
              if (btest(conf(k,j,i),n)) then
                 imo0=n+1+(k-1)*n_bits
                 print*,imo0
              endif
           enddo
        enddo
     enddo
     
  enddo
  
  !print*,''
  !do i=1,nmo0
  !   print*,i,maxval(smat(i,:))
  !enddo
  !stop

 
  !
  ! To begin with, we can just assume a one-to-one MO mapping
  ! i.e., MO rotation angles of zero or pi/2
  !
  
  !
  ! We should stick the generated confs into a hash table
  ! to avoid duplicates
  !

  !
  ! We should check the irreps of the generated confs...
  !
  
  return
  
end subroutine ref_space_propagate
