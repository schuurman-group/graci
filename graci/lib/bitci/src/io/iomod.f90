module iomod

  save

  !
  ! Error message to be written when terminating the program
  !
  character(len=200) :: errmsg
  
contains

!#######################################################################
! cstrlen: returns the length of a C character array
!#######################################################################
  function cstrlen(cstring)

    use constants
    use iso_c_binding, only: C_CHAR, C_NULL_CHAR
    
    implicit none

    integer(is)                        :: cstrlen
    character(kind=C_CHAR), intent(in) :: cstring(*)

    cstrlen=0
    do
       cstrlen=cstrlen+1
       if (cstring(cstrlen) == C_NULL_CHAR) exit
    enddo
    
    return
    
  end function cstrlen
    
!#######################################################################
! c2fstr: converts a C character array to a fortran fixed-length string
!#######################################################################
  subroutine c2fstr(cstring,fstring,length)

    use constants 
    use iso_c_binding, only: C_CHAR, C_NULL_CHAR

    implicit none


    character(kind=C_CHAR), intent(in) :: cstring(length)
    character(len=255), intent(out)    :: fstring
    integer(is), intent(in)            :: length
    integer(is)                        :: i

    fstring=''

    do i = 1,255
       if (cstring(i) == C_NULL_CHAR) exit
       fstring(i:i) = cstring(i) 
    enddo
  
    fstring = adjustl(fstring)
    
    return
    
  end subroutine c2fstr
  
!#######################################################################
! f2Cstr: converts a fortran fixed-length string to a C character array
!#######################################################################
  subroutine f2cstr(fstring,cstring,length)

    use constants 
    use iso_c_binding, only: C_CHAR, C_NULL_CHAR

    implicit none

    character(len=255), intent(in)      :: fstring
    character(kind=C_CHAR), intent(out) :: cstring(length)
    integer(is), intent(in)             :: length
    integer(is)                         :: i
    
    do i=1,len_trim(fstring)
       cstring(i)=fstring(i:i)
    enddo

    i=len_trim(fstring)+1
    cstring(i)=C_NULL_CHAR
    
    return
    
  end subroutine f2cstr
    
!#######################################################################
! freeunit: determines the first free unit number
!#######################################################################
  subroutine freeunit(unit)

    use constants
      
    implicit none

    integer(is) :: unit,i
    logical     :: lopen
      
    !
    ! Find the first free unit
    ! N.B. Save the first 20 io units for standard files
    !
    do i=20,1000
       inquire(unit=i,opened=lopen)
       if (.not.lopen) then
          unit=i
          exit
       endif
    enddo

    return

  end subroutine freeunit

!#######################################################################
! error_control: writes the passed string to the screen and, if open,
!                the log file, then terminates the program
!#######################################################################
  subroutine error_control
      
    implicit none
      
    logical :: lopen

    !
    ! Write the error message to the screen and terminate
    ! the program
    !
    write(6,'(/,2x,a,/)') trim(errmsg)
    stop
    
  end subroutine error_control
  
!#######################################################################
! scratch_name: for a given string filename, returns the string
!               scratchdir/filename
!#######################################################################
  subroutine scratch_name(filename,path)

    use constants
    use bitglobal
    
    implicit none

    character(len=*) :: filename,path

    path=trim(scratchdir)//'/'//trim(filename)
    
    return
    
  end subroutine scratch_name
    
!#######################################################################
! register_scratch_file: registers a new scratch file unit and filename
!#######################################################################
  subroutine register_scratch_file(unit,filename)

    use constants
    use bitglobal
    
    implicit none

    integer(is), intent(inout)      :: unit
    character(len=*), intent(in)    :: filename

    integer(is)                     :: olddim
    integer(is), allocatable        :: iswap(:)
    character(len=255), allocatable :: aswap(:)
    
    !
    ! Increment the number of scratch files in use
    !
    nscratch=nscratch+1

    !
    ! If we are past the dimensions of the scratch file arrays,
    ! then re-allocate them
    !
    if (nscratch > maxunits) then
       olddim=maxunits
       allocate(iswap(olddim),aswap(olddim))
       iswap=scrunit
       aswap=scrname
       maxunits = maxunits+100
       deallocate(scrunit,scrname)
       allocate(scrunit(maxunits),scrname(maxunits))
       scrunit=0
       scrname=''
       scrunit(1:olddim)=iswap
       scrname(1:olddim)=aswap
       deallocate(iswap,aswap)
    endif
    
    !
    ! Register the new scratch file
    !
    unit=nscratch
    scrunit(nscratch)=1000+nscratch
    scrname(nscratch)=filename
    
    return
    
  end subroutine register_scratch_file

!#######################################################################
! read_det_file: Reads the determinants and offset information from
!                the scratch file numbered scrnum
!#######################################################################
  subroutine read_det_file(scrnum,ndet,offdim_a,offdim_b,nsym,da,db,&
       nunique_a,nunique_b,offset_a,offset_b,mapab)

    use constants
    use bitglobal
    
    implicit none

    integer(is), intent(in)  :: scrnum
    integer(is), intent(out) :: ndet,offdim_a,offdim_b
    integer(is), intent(out) :: nsym(0:nirrep-1)
    integer(ib), allocatable :: da(:,:,:),db(:,:,:)
    integer(is), intent(out) :: nunique_a(0:nirrep-1),&
                                nunique_b(0:nirrep-1)
    integer(is), allocatable :: offset_a(:,:),offset_b(:,:)
    integer(is), allocatable :: mapab(:)
    integer(is)              :: iscratch
    
    !
    ! Open the scratch file
    !
    iscratch=scrunit(scrnum)
    open(iscratch,file=scrname(scrnum),form='unformatted',&
         status='old')
    
    !
    ! Number of determinants
    !
    read(iscratch) ndet
    
    !
    ! Leading dimensions of the offset arrays
    !
    read(iscratch) offdim_a
    read(iscratch) offdim_b

    !
    ! Number of determinants per irrep
    !
    read(iscratch) nsym

    !
    ! Determinants in alpha- and beta-major order
    !
    allocate(da(n_int,2,ndet))
    allocate(db(n_int,2,ndet))
    read(iscratch) da
    read(iscratch) db
    
    !
    ! Number of unique alpha and beta strings per irrep
    !
    read(iscratch) nunique_a
    read(iscratch) nunique_b
    
    !
    ! Offsets for the alpha- and beta-major ordered determinants
    !
    allocate(offset_a(offdim_a,0:nirrep-1))
    allocate(offset_b(offdim_b,0:nirrep-1))
    read(iscratch) offset_a
    read(iscratch) offset_b
    
    !
    ! Mapping between the alpha- and beta-major ordered determinants
    !
    allocate(mapab(ndet))
    read(iscratch) mapab
    
    !
    ! Close the scratch file
    !
    close(iscratch)
    
    return
    
  end subroutine read_det_file

!#######################################################################
! read_dets_alpha_major: Reads only the determinants in alpha-major
!                        order from the scratch file numbered scrnum
!#######################################################################
  subroutine read_dets_alpha_major(scrnum,ndet,da)

    use constants
    use bitglobal
    
    implicit none

    integer(is), intent(in)  :: scrnum
    integer(is), intent(out) :: ndet
    integer(ib), allocatable :: da(:,:,:)
    integer(is)              :: iscratch
    integer(is)              :: idum
    integer(is)              :: nsym(0:nirrep-1)
    
    !
    ! Open the scratch file
    !
    iscratch=scrunit(scrnum)
    open(iscratch,file=scrname(scrnum),form='unformatted',&
         status='old')

    !
    ! Number of determinants
    !
    read(iscratch) ndet

    !
    ! Read paset the leading dimensions of the offset arrays and
    ! the number of determinants per irrep
    !
    read(iscratch) idum
    read(iscratch) idum
    read(iscratch) nsym

    !
    ! Determinants in alpha- and beta-major order
    !
    allocate(da(n_int,2,ndet))
    read(iscratch) da

    !
    ! Close the scratch file
    !
    close(iscratch)
    
    return
    
  end subroutine read_dets_alpha_major
    
!#######################################################################
! read_ref_confs: Reads the reference configuration information from
!                 the scratch file numbered scrnum
!#######################################################################
  subroutine read_ref_confs(scrnum,nconf,n_int_I,nmoI,nmoE,conf,sop,&
       m2c,c2m)

    use constants
    use bitglobal
    
    implicit none
    
    integer(is), intent(in)  :: scrnum
    integer(is), intent(out) :: nconf,n_int_I,nmoI,nmoE
    integer(ib), allocatable :: conf(:,:,:),sop(:,:,:)
    integer(is), intent(out) :: m2c(nmo),c2m(nmo)
    integer(is)              :: iscratch

    !
    ! Open the scratch file
    !
    iscratch=scrunit(scrnum)
    open(iscratch,file=scrname(scrnum),form='unformatted',&
         status='old')
    
    !
    ! Number of reference space configurations
    !
    read(iscratch) nconf
    
    !
    ! Number of 64-bit integers needed to represent each reference
    ! space configuration and SOP bit string
    !
    read(iscratch) n_int_I
    
    !
    ! Number of internal and external MOs
    !
    read(iscratch) nmoI
    read(iscratch) nmoE
    
    !
    ! Configurations
    !
    allocate(conf(n_int_I,2,nconf))
    read(iscratch) conf
    
    !
    ! SOP for the singly-occupied orbitals
    !
    allocate(sop(n_int_I,2,nconf))
    read(iscratch) sop

    ! MO mapping arrays
    read(iscratch) m2c
    read(iscratch) c2m
    
    !
    ! Close the scratch file
    !
    close(iscratch)
    
    return
    
  end subroutine read_ref_confs

!#######################################################################
! read_nconf0: Reads the number of reference configurations from the
!              the scratch file numbered scrnum
!#######################################################################
  subroutine read_nconf0(scrnum,nconf0)

    use constants
    use bitglobal
    
    implicit none

    integer(is), intent(in)  :: scrnum
    integer(is), intent(out) :: nconf0
    integer(is)              :: iscratch    
    
    !
    ! Open the scratch file
    !
    iscratch=scrunit(scrnum)
    open(iscratch,file=scrname(scrnum),form='unformatted',&
         status='old')
    
    !
    ! Number of reference space configurations
    !
    read(iscratch) nconf0

    !
    ! Close the scratch file
    !
    close(iscratch)
    
    return
    
  end subroutine read_nconf0

!#######################################################################
! read_intext_dims: Reads the internal and external space dimensions
!                   from the the reference space configuration scratch
!                   file numbered scrnum
!#######################################################################
  subroutine read_n_int_I(scrnum,n_int_I,nmoI,nmoE)

    use constants
    use bitglobal
    
    implicit none

    integer(is), intent(in)  :: scrnum
    integer(is), intent(out) :: n_int_I,nmoI,nmoE
    integer(is)              :: iscratch,nconf

    !
    ! Open the scratch file
    !
    iscratch=scrunit(scrnum)
    open(iscratch,file=scrname(scrnum),form='unformatted',&
         status='old')
    
    !
    ! Number of reference space configurations
    !
    read(iscratch) nconf
    
    !
    ! Number of 64-bit integers needed to represent each reference
    ! space configuration and SOP bit string
    !
    read(iscratch) n_int_I
    
    !
    ! Number of internal and external MOs
    !
    read(iscratch) nmoI
    read(iscratch) nmoE

    !
    ! Close the scratch file
    !
    close(iscratch)
    
    return
    
  end subroutine read_n_int_I
    
!#######################################################################
! read_conf0_all: reads the reference space configurations for all
!                 irreps from disk
!######################################################################
  subroutine read_ref_confs_all(scrnum,nconf0,maxconf0,conf0h,sop0h,&
       n_int_I,nmoI,nmoE,m2c,c2m,nroots)
    
    use constants
    use bitglobal
    
    implicit none
    
    ! Reference space configuration scratch file numbers
    integer(is), intent(in)  :: scrnum(0:nirrep-1)
    
    ! Number of reference space configurations for each irrep
    integer(is), intent(out) :: nconf0(0:nirrep-1)
    integer(is), intent(out) :: maxconf0

    ! Number of roots for each irrep
    integer(is), intent(in)  :: nroots(0:nirrep-1)
  
    ! Reference space configurations
    integer(ib), allocatable :: conf0h(:,:,:,:)
    integer(ib), allocatable :: sop0h(:,:,:,:)
    integer(is), intent(out) :: n_int_I,nmoI,nmoE

    ! MO mappings
    integer(is), intent(out) :: m2c(nmo),c2m(nmo)
    
    ! Everything else
    integer(is)              :: i,k,iscratch,idum
    
!----------------------------------------------------------------------
! Read the number of reference space configurations of each symmetry
!----------------------------------------------------------------------
    nconf0=0
    do i=0,nirrep-1
       if (nroots(i) > 0) call read_nconf0(scrnum(i),nconf0(i))
    enddo

!----------------------------------------------------------------------
! Read the internal and external space dimensions from disk
!----------------------------------------------------------------------
    do i=0,nirrep-1
       if (nroots(i) > 0) then
          call read_n_int_I(scrnum(i),n_int_I,nmoI,nmoE)
          exit
       endif
    enddo

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    maxconf0=maxval(nconf0)

    allocate(conf0h(n_int_I,2,maxconf0,0:nirrep-1))
    conf0h=0_ib

    allocate(sop0h(n_int_I,2,maxconf0,0:nirrep-1))
    sop0h=0_ib
    
!----------------------------------------------------------------------
! Read the reference space configurations from disk
!----------------------------------------------------------------------
    ! Loop over irreps
    do i=0,nirrep-1

       ! Cycle if there are no roots for this irrep
       if (nroots(i) == 0) cycle

       ! Open the scratch file
       iscratch=scrunit(scrnum(i))
       open(iscratch,file=scrname(scrnum(i)),form='unformatted',&
            status='old')

       ! Skip past the dimensions that we have already read
       read(iscratch) idum ! nconf
       read(iscratch) idum ! n_int_I
       read(iscratch) idum ! nmoI
       read(iscratch) idum ! nmoE

       ! Configurations
       read(iscratch) conf0h(:,:,1:nconf0(i),i)

       ! SOPs
       read(iscratch) sop0h(:,:,1:nconf0(i),i)

       ! MO mapping arrays
       read(iscratch) m2c
       read(iscratch) c2m
       
       ! Close the scratch file
       close(iscratch)
       
    enddo
    
    return

  end subroutine read_ref_confs_all

!######################################################################
! read_energies: reads the eigenvalues for a single irrep from the
!                scratch file numbered scrnum
!######################################################################
  subroutine read_energies(scrnum,nroots,ener)
  
    use constants
    use bitglobal
    
    implicit none
    
    integer(is), intent(in) :: scrnum,nroots
    integer(is)             :: iscratch
    integer(is)             :: ndum,nroots1
    real(dp), intent(out)   :: ener(nroots)

    !
    ! Open the scratch file
    !
    iscratch=scrunit(scrnum)
    open(iscratch,file=scrname(scrnum),form='unformatted',&
         status='old')
  
    !
    ! Read past the N-electron basis dimension
    !
    read(iscratch) ndum
  
    !
    ! Check on the number of roots
    !
    read(iscratch) nroots1
    if (nroots /= nroots1) then
       write(6,'(/,2x,a)') 'Error in read_energies: inconsistent '&
            //'number of roots'
       stop
    endif
    
    !
    ! Read in the eigenvalues
    !
    read(iscratch) ener

    !
    ! Close the scratch file
    !
    close(iscratch)
    
    return
    
  end subroutine read_energies
  
!######################################################################
! read_single_vector: reads a single eigenvector indexed n from the
!                           scratch file numbered scrnum
!######################################################################
  subroutine read_single_vector(scrnum,vec,dim,n)

    use constants
    use bitglobal
    
    implicit none

    integer(is), intent(in) :: scrnum,dim,n
    integer(is)             :: iscratch
    integer(is)             :: ndum,nroots,i
    real(dp), intent(out)   :: vec(dim)
    real(dp), allocatable   :: ener(:)
    
    !
    ! Open the scratch file
    !
    iscratch=scrunit(scrnum)
    open(iscratch,file=scrname(scrnum),form='unformatted',&
         status='old')
  
    !
    ! Consistency check on the CSF basis dimension
    !
    read(iscratch) ndum
    if (ndum /= dim) then
       errmsg='Error in read_single_vector: '&
            //'wrong CSF basis dimension encountered'
       call error_control
    endif

    !
    ! Number of roots
    !
    read(iscratch) nroots

    !
    ! Consistency check on the number of roots
    !
    if (n > nroots) then
       errmsg='Error in read_single_vector: '&
            //'inconsistent root number requested'
       call error_control       
    endif

    !
    ! Read past the energies
    !
    allocate(ener(nroots))
    read(iscratch) ener

    !
    ! Read in the requested eigenvector
    !
    do i=1,n
       read(iscratch) vec       
    enddo
    
    !
    ! Close the scratch file
    !
    close(iscratch)
    
    return
    
  end subroutine read_single_vector
  
!######################################################################
! read_all_eigenpairs: reads all the eigenpairs from the scratch file
!                      numbered scrnum
!######################################################################
  subroutine read_all_eigenpairs(scrnum,vec,ener,dim,nroots)

    use constants
    use bitglobal
    
    implicit none

    integer(is), intent(in) :: scrnum,dim,nroots
    integer(is)             :: iscratch
    integer(is)             :: ndum,i
    real(dp), intent(out)   :: vec(dim,nroots),ener(nroots)
    
    !
    ! Open the scratch file
    !
    iscratch=scrunit(scrnum)
    open(iscratch,file=scrname(scrnum),form='unformatted',&
         status='old')

    !
    ! Consistency check on the CSF basis dimension
    !
    read(iscratch) ndum
    if (ndum /= dim) then
       errmsg='Error in read_all_eigenpairs: '&
            //'wrong CSF basis dimension encountered'
       call error_control
    endif

    !
    ! Consistency check on the number of roots
    !
    read(iscratch) ndum
    if (ndum /= nroots) then
       errmsg='Error in read_all_eigenpairs: '&
            //'wrong no. roots encountered'
       call error_control
    endif

    !
    ! Read in the energies
    !
    read(iscratch) ener

    !
    ! Read in the eigenvectors
    !
    do i=1,nroots
       read(iscratch) vec(:,i)
    enddo

    !
    ! Close the scratch file
    !
    close(iscratch)
    
    return
    
  end subroutine read_all_eigenpairs
    
!######################################################################
! read_some_eigenpairs: reads a subset of the eigenpairs (as specified
!                       in the array iroots) from the scratch file
!                       numbered scrnum
!######################################################################
  subroutine read_some_eigenpairs(scrnum,vec,ener,dim,nroots,iroots)

    use constants
    use bitglobal
    
    implicit none

    integer(is), intent(in) :: scrnum,dim,nroots
    integer(is)             :: iscratch
    integer(is)             :: ndum,nentries,i,n
    integer(is)             :: ipos(1)
    integer(is), intent(in) :: iroots(nroots)
    real(dp), intent(out)   :: vec(dim,nroots),ener(nroots)
    real(dp), allocatable   :: ener1(:),vec1(:)
    
    !
    ! Open the scratch file
    !
    iscratch=scrunit(scrnum)
    open(iscratch,file=scrname(scrnum),form='unformatted',&
         status='old')

    !
    ! Consistency check on the CSF basis dimension
    !
    read(iscratch) ndum
    if (ndum /= dim) then
       errmsg='Error in read_all_eigenpairs: '&
            //'wrong CSF basis dimension encountered'
       call error_control
    endif

    !
    ! Number of roots on file
    !
    read(iscratch) nentries

    !
    ! Allocate working arrays
    !
    allocate(ener1(nentries))
    allocate(vec1(dim))
    
    !
    ! Read in the energies
    !
    read(iscratch) ener1

    !
    ! Read in the eigenvectors
    !
    n=0
    do i=1,nentries
       ! Read in the next eigenvector
       read(iscratch) vec1(:)
       ! Save the eigenvector (and energy) if it is in the subset
       ! of requested roots
       ipos=findloc(iroots,value=i)
       if (ipos(1) /= 0) then
          n=n+1
          vec(:,n)=vec1(:)
          ener(n)=ener1(i)
       endif
    enddo

    !
    ! Close the scratch file
    !
    close(iscratch)
    
    return
    
  end subroutine read_some_eigenpairs
    
!######################################################################
  
end module iomod
