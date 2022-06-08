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
  !integer(is)                        :: m2c(nmo0),c2m(nmo0)
  integer(is), allocatable           :: m2c0(:),c2m0(:)
  
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

  call read_ref_confs(scrnum0,nconf0,n_int_I,nmoI,nmoE,conf_r0,sop_r0,&
       m2c0,c2m0)

  !
  ! Note that we are going to need to put the r0 confs into the
  ! 'canonical' MO ordering before continuing
  ! Also, the conf_r0 array is going to have to be re-dimensioned
  ! to (***n_int***,2,nconf0) before doing this, as its leading
  ! dimension is currently n_int_I...
  !


  print*,''
  print*,'Note that we are going to need to put the r0 confs into the'
  print*,'canonical MO ordering before continuing'
  print*,'Also, the conf_r0 array is going to have to be re-dimensioned'
  print*,'to (***n_int***,2,nconf0) before doing this, as its leading'
  print*,'dimension is currently n_int_I...'
  STOP
  
  return
  
end subroutine ref_space_propagate
