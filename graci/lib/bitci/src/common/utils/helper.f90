!**********************************************************************
! Helper routines that return system information and determinants
! written to disk
!**********************************************************************

!######################################################################
! get_n_int: simply returns the value of n_int
!######################################################################
#ifdef CBINDING
subroutine get_n_int(n) bind(c,name="get_n_int")
#else
subroutine get_n_int(n)
#endif
  
  use constants
  use bitglobal
  
  implicit none
  
  integer(is) :: n

  n=n_int
  
  return

end subroutine get_n_int

!######################################################################
! retrieve_base_det: returns the base determinant
!######################################################################
#ifdef CBINDING
subroutine retrieve_base_det(d) bind(c,name="retrieve_base_det")
#else
subroutine retrieve_base_det(d)
#endif
  
  use constants
  use bitglobal
  
  integer(ib) :: d(n_int,2)
  
  d=det0
  
  return
  
end subroutine retrieve_base_det

!######################################################################
! retrieve_mrci_dets: returns the MRCI determinants (Soon to be
!                     deprecated!)
!######################################################################
#ifdef CBINDING
subroutine retrieve_mrci_dets(d,nmrci) &
     bind(c,name="retrieve_mrci_dets")
#else
subroutine retrieve_mrci_dets(d,nmrci)
#endif
  
  use constants
  use bitglobal
  use iomod
  
  integer(is), intent(in)  :: nmrci
  integer(ib), intent(out) :: d(n_int,2,nmrci)
  integer(is)              :: bufdim,nrecords
  integer(ib), allocatable :: buffer(:,:,:)
  integer(is)              :: nchk
  integer(is)              :: detcnt,ndet,i
  integer(is)              :: iscratch
  character(len=60)        :: mrcifile

  !
  ! Open the scratch file
  !
  call freeunit(iscratch)
  call scratch_name('detmrci',mrcifile)
  open(iscratch,file=mrcifile,form='unformatted',status='unknown')

  !
  ! Read the number of MRCI determinants
  !
  read(iscratch) nchk

  !
  ! Exit if the number of determinants written to disk is different
  ! to the size of the determinant array that was passed to us
  !
  if (nchk.ne.nmrci) then
     errmsg='Inconsistent numbers of determinants in retrieve_mrci_dets'
     call error_control
  endif
  
  !
  ! Read the buffer size and the number of records
  !
  read(iscratch) bufdim
  read(iscratch) nrecords

  !
  ! Allocate the buffer
  !
  allocate(buffer(n_int,2,bufdim))

  !
  ! Read in the determinants
  !
  detcnt=0
  do i=1,nrecords
     read(iscratch) ndet,buffer
     d(:,:,detcnt+1:detcnt+ndet)=buffer(:,:,1:ndet)
     detcnt=detcnt+ndet
  enddo

  !
  ! Deallocate the buffer
  !
  deallocate(buffer)
  
  !
  ! Close the scratch file
  !
  close(iscratch)
    
  return
    
end subroutine retrieve_mrci_dets

!#######################################################################
! retrieve_energies: Reads the eigenvalues from the scratch file
!                    numbered scrnum
!#######################################################################
#ifdef CBINDING
subroutine retrieve_energies(scrnum,nroots,ener) &
     bind(c,name="retrieve_energies")
#else
subroutine retrieve_energies(scrnum,nroots,ener)
#endif
  
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
     write(6,'(/,2x,a)') 'Error in retrieve_energies: inconsistent '&
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
  
end subroutine retrieve_energies

!######################################################################
! retieve_some_energies: Returns a subset of the eigenvalues from the
!                        scratch file numbered scrnum
!######################################################################
#ifdef CBINDING
subroutine retrieve_some_energies(scrnum,nroots,ener,iroots) &
     bind(c,name="retrieve_some_energies")
#else
subroutine retrieve_some_energies(scrnum,nroots,ener,iroots)
#endif
  
  use constants
  use bitglobal
  
  implicit none
  
  integer(is), intent(in) :: scrnum,nroots
  integer(is), intent(in) :: iroots(nroots)
  real(dp), intent(out)   :: ener(nroots)
  integer(is)             :: iscratch,ndum,nentries,n,i
  integer(is)             :: ipos(1)
  real(dp), allocatable   :: ener1(:)
  
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
  ! Number of roots on file
  !
  read(iscratch) nentries

  !
  ! Allocate working array
  !
  allocate(ener1(nentries))
  
  !
  ! Read in the eigenvalues
  !
  read(iscratch) ener1

  !
  ! Save the requested eigenvalues
  !
  n=0
  ! Loop over all eigenvalues
  do i=1,nentries
     ! Save the eigenvalue if it is in the subset of requested roots
     ipos=findloc(iroots,value=i)
     if (ipos(1) /= 0) then
        n=n+1
        ener(n)=ener1(i)
     endif
  enddo
  
  !
  ! Close the scratch file
  !
  close(iscratch)
  
  
  return
  
end subroutine retrieve_some_energies
  
!######################################################################
! retrieve_filename: given a scratch file number, returns the
!                    associated filename
!######################################################################
#ifdef CBINDING
subroutine retrieve_filename(scrnum,filename1) &
     bind(c,name="retrieve_filename")
#else
subroutine retrieve_filename(scrnum,filename1)
#endif

  use constants
  use bitglobal
  use iomod
  use iso_c_binding, only: C_CHAR, C_NULL_CHAR
  
  implicit none

  ! Scratch file number
  integer(is), intent(in)             :: scrnum

  ! Filename
#ifdef CBINDING
  character(kind=C_CHAR), intent(out) :: filename1(*)
  character(len=255)                  :: filename
  integer(is)                         :: length
#else
  character(len=*), intent(out)       :: filename1
  character(len=255)                  :: filename
#endif
  
  !
  ! Filename
  !
  filename=scrname(scrnum)

  !
  ! If C bindings are on, then convert to the C char type
  !
#ifdef CBINDING
  ! Length of the C char array
  length=cstrlen(filename1)
  
  ! Convert to a fortran fixed-length string
  call f2cstr(filename,filename1,length)
#else
  filename1=filename
#endif
  
  return
  
end subroutine retrieve_filename
