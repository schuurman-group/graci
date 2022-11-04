!**********************************************************************
! Routines for the generation of reference space configurations via
! the approximate re-expression of those of another geometry in terms
! of the MOs of the current one
!**********************************************************************

#ifdef CBINDING
subroutine ref_space_propagate(nroots,nmo0,nmo1,smat,conffile0_in,&
     nconf,confscr) bind(c,name="ref_space_propagate")
#else
subroutine ref_space_propagate(nroots,nmo0,nmo1,smat,conffile0_in,&
     nconf,confscr)
#endif
  
  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use dethash
  use mrciutils
  use iomod

  implicit none

  ! No. roots
  integer(is), intent(in)            :: nroots(0:nirrep-1)
  
  ! MO basis dimensions
  integer(is), intent(in)            :: nmo0,nmo1

  ! MO overlaps
  real(dp), intent(in)               :: smat(nmo0,nmo1)

  ! Scratch file name for the ref space configurations to propagate
  ! forwards
#ifdef CBINDING
  character(len=1,kind=C_CHAR), intent(in) :: conffile0_in(255,0:nirrep-1)
  character(len=255)                       :: conffile0(0:nirrep-1)
  integer(is)                              :: length
#else
  character(len=255), intent(in)           :: conffile0_in(0:nirrep-1)
  character(len=255)                       :: conffile0(0:nirrep-1
#endif

  ! Number of generated reference configurations
  integer(is), intent(out)           :: nconf(0:nirrep-1)

  ! Generated reference configuration scratch file numbers
  integer(is), intent(out)           :: confscr(0:nirrep-1)

  ! Reference confs and SOPs
  integer(is)                        :: scrnum0(0:nirrep-1)
  integer(is)                        :: nconf0(0:nirrep-1)
  integer(is)                        :: maxconf0,maxconf
  integer(ib), allocatable           :: conf_r0(:,:,:,:),sop_r0(:,:,:,:)
  integer(ib), allocatable           :: conf(:,:,:,:),sop(:,:,:,:)
  integer(is), allocatable           :: m2c(:),c2m(:)

  ! Internal MO space information
  integer(is)                        :: nmoI,nmoE,n_int_I
  integer(is)                        :: Ilist(nmo),Elist(nmo)
  
  ! Hash table
  type(dhtbl)                        :: h(0:nirrep-1)
  integer(is)                        :: initial_size
  integer(ib)                        :: key(n_int,2),key1(n_int,2)

  ! Working arrays
  integer(is)                        :: nsocc,ndocc
  integer(is), allocatable           :: socc(:),docc(:)
  integer(is), allocatable           :: imap(:)
  real(dp), allocatable              :: abs_smatT(:,:)

  ! IO variables
  integer(is)                        :: iscratch
  character(len=2)                   :: amult,asym
  character(len=250)                 :: filename
  
  ! Everything else
  integer(is)                        :: irrep,i,n,imo0,imo(1),ntot
  integer(is)                        :: istart,iend
  integer(ib), allocatable           :: work(:,:,:)
  
!----------------------------------------------------------------------
! If C bindings are on, then convert the previous geometry ref conf
! file names from the C char type to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  do irrep=0,nirrep-1
     length=cstrlen(conffile0_in(:,irrep))
     call c2fstr(conffile0_in(:,irrep),conffile0(irrep),length)
  enddo
#else
  do irrep=0,nirrep-1
     conffile0(irrep)=adjustl(trim(conffile0_in(irrep)))
  enddo
#endif

!----------------------------------------------------------------------
! Register the previous geometry ref conf scratch files
!----------------------------------------------------------------------
  do irrep=0,nirrep-1
     call register_scratch_file(scrnum0(irrep),conffile0(irrep))
  enddo

!----------------------------------------------------------------------
! Read in the previous geometry ref confs
!----------------------------------------------------------------------
  allocate(m2c(nmo0),c2m(nmo0))

  call read_ref_confs_all(scrnum0,nconf0,maxconf0,conf_r0,sop_r0,&
       n_int_I,nmoI,nmoE,m2c,c2m,nroots)

!----------------------------------------------------------------------
! Put the ref confs into the 'canonical' MO order
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(conf(n_int,2,maxconf0,0:nirrep-1))
  allocate(sop(n_int,2,maxconf0,0:nirrep-1))
  conf=0_ib
  sop=0_ib  

  ! Resize the conf and sop arrays
  conf(1:n_int_I,:,:,:)=conf_r0
  sop(1:n_int_I,:,:,:)=sop_r0
  
  ! Reorder the conf and sop arrays
  do irrep=0,nirrep-1
     n=nconf0(irrep)
     call reorder_confs(m2c,conf(:,:,1:n,irrep),n)
     call reorder_confs(m2c,sop(:,:,1:n,irrep),n)
  enddo

!----------------------------------------------------------------------
! Initialise the hash tables
!----------------------------------------------------------------------
  do irrep=0,nirrep-1
     initial_size=nconf0(irrep)
     call h(irrep)%initialise_table(initial_size)
  enddo

!----------------------------------------------------------------------
! Construct the MO mapping array
!----------------------------------------------------------------------
  allocate(imap(nmo0))
  allocate(abs_smatT(nmo,nmo0))
  imap=0
  abs_smatT=0.0d0

  ! Absolute value of the transpose of the MO overlap matrix
  abs_smatT=abs(transpose(smat))

  ! Loop over ref geometry MOs
  do i=1,nmo0

     ! Current MO with the greatest overlap with the ref MO
     imo=maxloc(abs_smatT(:,i))
     imap(i)=imo(1)

     ! Remove the selected current MO from the list
     abs_smatT(imo(1),:)=0.0d0
     
  enddo
  
!----------------------------------------------------------------------
! Construct the current geometry ref confs
!----------------------------------------------------------------------
! Important: as a first pass, we will assume a one-to-one mapping
!            between the previous and current MOs
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(socc(nmo0))
  allocate(docc(nmo0))
  socc=0
  docc=0
  
  ! Loop over irreps
  do irrep=0,nirrep-1

     ! Loop over the previous goemetry ref confs for this irrep
     do i=1,nconf0(irrep)

        ! Get the lists of singly- and doubly-occupied MOs
        call sop_socc_list(sop(:,:,i,irrep),n_int,socc,nmo0,nsocc)
        call sop_docc_list(sop(:,:,i,irrep),n_int,docc,nmo0,ndocc)

        ! Initialisation
        key=0_ib
     
        ! Doubly-occupied MOs
        do n=1,ndocc

           ! Previous geometry MO index
           imo0=docc(n)

           ! Current geometry MO with the greatest overlap
           ! with imo0
           imo(1)=imap(imo0)
           
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
           imo(1)=imap(imo0)
           
           ! Create the electron
           key1=create_electron(key,n_int,imo(1))
           key=key1
        
        enddo
        
        ! Insert the conf into the hash table
        call h(irrep)%insert_key(key)
        
     enddo
     
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
  do irrep=0,nirrep-1
     nconf(irrep)=h(irrep)%n_keys_stored
  enddo

  ! Allocate arrays
  maxconf=maxval(nconf)
  
  ! Allocate arrays
  allocate(conf(n_int,2,maxconf,0:nirrep-1))
  allocate(sop(n_int,2,maxconf,0:nirrep-1))
  conf=0_ib
  sop=0_ib

  ! Fetch the confs
  do irrep=0,nirrep-1
     call h(irrep)%retrieve_keys(conf(:,:,1:nconf(irrep),irrep))
  enddo

!----------------------------------------------------------------------
! Generate the SOPs
!----------------------------------------------------------------------
  do irrep=0,nirrep-1
     do i=1,nconf(irrep)
        sop(:,:,i,irrep)=conf_to_sop(conf(:,:,i,irrep),n_int)
     enddo
  enddo

!----------------------------------------------------------------------
! Get the list of internal MOs
! Here, Ilist is the list of internal MO indices, Elist is the list of
! external MO indices, nmoI is the number of internal MOs, and nmoE
! is the number of external MOs
!----------------------------------------------------------------------
  ! Work array of all confs
  ntot=sum(nconf)
  allocate(work(n_int,2,ntot))
  work=0_ib

  istart=1
  do irrep=0,nirrep-1
     iend=istart+nconf(irrep)-1
     work(:,:,istart:iend)=conf(:,:,1:nconf(irrep),irrep)
     istart=iend+1
  enddo
  
  ! Get the lists of internal and external MOs
  call get_internal_external_mos(work,ntot,nmoI,nmoE,Ilist,Elist)

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
  do irrep=0,nirrep-1
     n=nconf(irrep)
     call reorder_confs(c2m,conf(:,:,1:n,irrep),n)
     call reorder_confs(c2m,sop(:,:,1:n,irrep),n)
  enddo

!----------------------------------------------------------------------
! Check on the irreps generated by the confs
!----------------------------------------------------------------------
  do irrep=0,nirrep-1
     do i=1,nconf(irrep)
        if (sop_sym_mrci(sop(:,:,i,irrep),m2c) /= irrep) then
           errmsg='Error in ref_space_propagate: ' &
                //'incorrect irrep encountered'
           call error_control
        endif
     enddo
  enddo

!----------------------------------------------------------------------
! Write the reference space configuration information to disk
!----------------------------------------------------------------------
  do irrep=0,nirrep-1
     
     ! Scratch file name
     write(amult,'(i0)') imult
     write(asym,'(i0)') irrep
     call scratch_name('refconf.mult'//trim(amult)//'.sym' &
          //trim(asym),filename)

     ! Register the scratch file
     call register_scratch_file(confscr(irrep),filename)

     ! Open the scratch file
     iscratch=scrunit(confscr(irrep))
     open(iscratch,file=scrname(confscr(irrep)),form='unformatted',&
          status='unknown')

     ! Number of reference space configurations
     write(iscratch) nconf(irrep)
       
     ! Number of (n_bits)-bit integers needed to represent each
     ! reference space configuration and SOP bit string
     write(iscratch) n_int_I
  
     ! Number of internal and external MOs
     write(iscratch) nmoI
     write(iscratch) nmoE
  
     ! Configurations
     write(iscratch) conf(1:n_int_I,:,1:nconf(irrep),irrep)
  
     ! SOPs
     write(iscratch) sop(1:n_int_I,:,1:nconf(irrep),irrep)
     
     ! MO mapping arrays
     write(iscratch) m2c
     write(iscratch) c2m
  
     ! Close the scratch file
     close(iscratch)

  enddo
     
  return

end subroutine ref_space_propagate
