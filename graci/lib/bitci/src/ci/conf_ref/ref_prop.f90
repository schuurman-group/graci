!**********************************************************************
! Routines for the generation of reference space configurations via
! the approximate re-expression of those of another geometry in terms
! of the MOs of the current one
!**********************************************************************

!######################################################################
! 
!######################################################################
#ifdef CBINDING
subroutine ref_space_propagate(irrep,nmo0,nmo1,smat,conffile0_in,&
     nconf,confscr) bind(c,name="ref_space_propagate")
#else
subroutine ref_space_propagate(irrep,nmo0,nmo1,smat,conffile0_in,&
     nconf,confscr)
#endif
  
  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use dethash
  use mrciutils
  use iomod

  implicit none

  ! Irrep
  integer(is), intent(in)            :: irrep
  
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

  ! Reference confs and SOPs
  integer(is)                        :: scrnum0
  integer(is)                        :: nconf0
  integer(ib), allocatable           :: conf_r0(:,:,:),sop_r0(:,:,:)
  integer(ib), allocatable           :: conf(:,:,:),sop(:,:,:)     
  integer(is), allocatable           :: m2c(:),c2m(:)

  ! Internal MO space information
  integer(is)                        :: nmoI,nmoE,n_int_I
  integer(is)                        :: Ilist(nmo),Elist(nmo)
  
  ! Hash table
  type(dhtbl)                        :: h
  integer(is)                        :: initial_size
  integer(ib)                        :: key(n_int,2),key1(n_int,2)

  ! Working arrays
  integer(is)                        :: nsocc,ndocc
  integer(is), allocatable           :: socc(:),docc(:)
  real(dp), allocatable              :: abs_smatT(:,:)

  ! IO variables
  integer(is)                        :: iscratch
  character(len=2)                   :: amult,asym
  character(len=60)                  :: filename
  
  ! Everything else
  integer(is)                        :: i,n,imo0,imo(1)
  
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
  allocate(m2c(nmo0),c2m(nmo0))

  call read_ref_confs(scrnum0,nconf0,n_int_I,nmoI,nmoE,conf_r0,&
       sop_r0,m2c,c2m)

!----------------------------------------------------------------------
! Put the ref confs into the 'canonical' MO order
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(conf(n_int,2,nconf0))
  allocate(sop(n_int,2,nconf0))
  conf=0_ib
  sop=0_ib  

  ! Resize the conf and sop arrays
  conf(1:n_int_I,:,:)=conf_r0
  sop(1:n_int_I,:,:)=sop_r0
  
  ! Reorder the conf and sop arrays
  call reorder_confs(m2c,conf,nconf0)
  call reorder_confs(m2c,sop,nconf0)
  
!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------
  ! Initialisation
  initial_size=nconf0
  call h%initialise_table(initial_size)
  
!----------------------------------------------------------------------
! Construct the current geometry ref confs
!----------------------------------------------------------------------
! Important: as a first pass, we will assume a one-to-one mapping
!            between the previous and current MOs
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(socc(nmo0))
  allocate(docc(nmo0))
  allocate(abs_smatT(nmo,nmo0))
  socc=0
  docc=0
  abs_smatT=0.0d0

  ! Absolute value of the transpose of the MO overlap matrix
  abs_smatT=abs(transpose(smat))
  
  ! Loop over previous geometry ref confs
  do i=1,nconf0

     ! Get the lists of singly- and doubly-occupied MOs
     call sop_socc_list(sop(:,:,i),n_int,socc,nmo0,nsocc)
     call sop_docc_list(sop(:,:,i),n_int,docc,nmo0,ndocc)

     ! Initialisation
     key=0_ib
     
     ! Doubly-occupied MOs
     do n=1,ndocc

        ! Previous geometry MO index
        imo0=docc(n)

        ! Current geometry MO with the greatest overlap
        ! with imo0
        imo=maxloc(abs_smatT(:,imo0))
        
        ! Create the electrons
        key1=create_electron(key,n_int,imo(1))
        key=create_electron(key1,n_int,imo(1))
        
     enddo

     ! Singly-occupied MOs
     do n=1,nsocc

        ! Previous geometry MO index
        imo0=socc(n)

        ! Current geometry MO with the greatest overlap
        ! with imo0
        imo=maxloc(abs_smatT(:,imo0))
        
        ! Create the electron
        key1=create_electron(key,n_int,imo(1))
        key=key1
        
     enddo

     ! Insert the conf into the hash table
     call h%insert_key(key)
     
  enddo

!----------------------------------------------------------------------
! Deallocate the old conf and sop arrays
!----------------------------------------------------------------------
  deallocate(conf)
  deallocate(sop)
  deallocate(m2c)
  deallocate(c2m)
  
!----------------------------------------------------------------------
! Retrieve the confs from the hash table
!----------------------------------------------------------------------
  ! Number of current geometry ref confs generated
  nconf=h%n_keys_stored

  ! Allocate arrays
  allocate(conf(n_int,2,nconf))
  allocate(sop(n_int,2,nconf))
  conf=0_ib
  sop=0_ib

  ! Fetch the confs
  call h%retrieve_keys(conf)
  
!----------------------------------------------------------------------
! Generate the SOPs
!----------------------------------------------------------------------
  do i=1,nconf
     sop(:,:,i)=conf_to_sop(conf(:,:,i),n_int)
  enddo

!----------------------------------------------------------------------
! Get the list of internal MOs
! Here, Ilist is the list of internal MO indices, Elist is the list of
! external MO indices, nmoI is the number of internal MOs, and nmoE
! is the number of external MOs
!----------------------------------------------------------------------
  ! Get the lists of internal and external MOs
  call get_internal_external_mos(conf,nconf,nmoI,nmoE,Ilist,Elist)

  ! No. (n_bits)-bit integers needed to represent each component
  ! of each conf/SOP
  n_int_I=(nmoI-1)/n_bits+1
  
!----------------------------------------------------------------------
! Construct the canonical-to-MRCI MO index mapping array
!----------------------------------------------------------------------
  allocate(m2c(nmo))
  allocate(c2m(nmo))
  m2c=0
  c2m=0
  
  call get_mo_mapping(nmoI,nmoE,Ilist,Elist,m2c,c2m)

!----------------------------------------------------------------------
! Re-arrange the reference space confs and SOPs such that the internal
! MOs come before the external MOs
!----------------------------------------------------------------------
  call reorder_confs(c2m,conf,nconf)
  call reorder_confs(c2m,sop,nconf)
  
!----------------------------------------------------------------------
! Check on the irreps generated by the confs
!----------------------------------------------------------------------
  do i=1,nconf
     if (sop_sym_mrci(sop(:,:,i),m2c) /= irrep) then
        errmsg='Error in ref_space_propagate: ' &
             //'incorrect irrep encountered'
        call error_control
     endif
  enddo

!----------------------------------------------------------------------
! Write the reference space configuration information to disk
!----------------------------------------------------------------------
  ! Scratch file name
  write(amult,'(i0)') imult
  write(asym,'(i0)') irrep
  call scratch_name('refconf.mult'//trim(amult)//'.sym'//trim(asym),&
       filename)

  ! Register the scratch file
  call register_scratch_file(confscr,filename)

  ! Open the scratch file
  iscratch=scrunit(confscr)
  open(iscratch,file=scrname(confscr),form='unformatted',&
       status='unknown')

  ! Number of reference space configurations
  write(iscratch) nconf
       
  ! Number of (n_bits)-bit integers needed to represent each
  ! reference space configuration and SOP bit string
  write(iscratch) n_int_I
  
  ! Number of internal and external MOs
  write(iscratch) nmoI
  write(iscratch) nmoE
  
  ! Configurations
  write(iscratch) conf(1:n_int_I,:,:)
  
  ! SOPs
  write(iscratch) sop(1:n_int_I,:,:)
  
  ! MO mapping arrays
  write(iscratch) m2c
  write(iscratch) c2m
  
  ! Close the scratch file
  close(iscratch)
  
  return
  
end subroutine ref_space_propagate
