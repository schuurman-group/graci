!**********************************************************************
! Helper routines that return information about a bitci calculation
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

  implicit none
  
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

  implicit none
  
  integer(is), intent(in)  :: nmrci
  integer(ib), intent(out) :: d(n_int,2,nmrci)
  integer(is)              :: bufdim,nrecords
  integer(ib), allocatable :: buffer(:,:,:)
  integer(is)              :: nchk
  integer(is)              :: detcnt,ndet,i
  integer(is)              :: iscratch
  character(len=250)       :: mrcifile

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

!######################################################################
! retrieve_nhpar: given a Hamiltonian name, returns the corresponding
!                 number of parameter values (including desel)
!######################################################################
#ifdef CBINDING
  subroutine retrieve_nhpar(ham1,npar) &
       bind(c,name="retrieve_nhpar")
#else
    subroutine retrieve_nhpar(ham1,npar)
#endif

  use constants
  use bitglobal
  use hparam
  use iomod
  use iso_c_binding, only: C_CHAR, C_NULL_CHAR
  
  implicit none

  ! Input: Hamiltonian name
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: ham1(*)
  character(len=255)                 :: ham
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: ham1
  character(len=255)                 :: ham
#endif

  ! Ouput: no. Hamiltonian parameters
  integer(is), intent(out)           :: npar

!----------------------------------------------------------------------
! If C bindings are on, then convert the Hamiltonian and calculation
! labels from the C char type to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(ham1)
  call c2fstr(ham1,ham,length)
#else
  ham=adjustl(trim(ham1))
#endif

!----------------------------------------------------------------------
! Number of Hamiltonian parameters
!----------------------------------------------------------------------
  select case(trim(ham))
     
  case('grimme')
     npar=size(grimme1)

  case('grimme_short')
     npar=size(grimme1_short)

  case('r2016')
     npar=size(r2016)

  case('r2016_short')
     npar=size(r2016_short)

  case('r2017')
     npar=size(r2017)

  case('r2017_short')
     npar=size(r2017_short)

  case('r2018')
     npar=size(r2018)

  case('r2018_short')
     npar=size(r2018_short)

  case('r2022')
     npar=size(r2022)

  case('qe8')
     npar=size(qe8)

  case('qe8_qhort')
     npar=size(qe8_short)

  case('cvs-qe8')
     npar=size(cvs_qe8)

  case default
     errmsg='Error in retrieve_nhpar: unrecognised Hamiltonian name'
     call error_control

  end select

!----------------------------------------------------------------------
! Account for desel
!----------------------------------------------------------------------
  npar=npar+1
  
  return
  
end subroutine retrieve_nhpar

!######################################################################
! retrieve_hparams: given a Hamiltonian name, returns the corresponding
!                   parameter values
!######################################################################
#ifdef CBINDING
subroutine retrieve_hpar(ham1,dim,params) &
     bind(c,name="retrieve_hpar")
#else
subroutine retrieve_hpar(ham1,dim,params)
#endif

  use constants
  use bitglobal
  use hparam
  use iomod
  use iso_c_binding, only: C_CHAR, C_NULL_CHAR
  
  implicit none

  ! Input: Hamiltonian name
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: ham1(*)
  character(len=255)                 :: ham
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: ham1
  character(len=255)                 :: ham
#endif

  ! Ouput: Hamiltonian parameters
  integer(is), intent(in)            :: dim
  real(dp), intent(out)              :: params(dim)

  ! Everything else
  integer(is)                        :: npar
  
!----------------------------------------------------------------------
! If C bindings are on, then convert the Hamiltonian and calculation
! labels from the C char type to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(ham1)
  call c2fstr(ham1,ham,length)
#else
  ham=adjustl(trim(ham1))
#endif

!----------------------------------------------------------------------
! Package up the Hamiltonian parameters excluding desel
!----------------------------------------------------------------------
  select case(trim(ham))
     
  case('grimme')
     npar=size(grimme1)+1
     if (npar > dim) goto 999
     params(1:npar-1)=grimme1
     params(npar)=1.d0

  case('grimme_short')
     npar=size(grimme1_short)+1
     if (npar > dim) goto 999
     params(1:npar-1)=grimme1_short
     params(npar)=0.8d0

  case('r2016')
     npar=size(r2016)+1
     if (npar > dim) goto 999
     params(1:npar-1)=r2016
     params(npar)=1.d0

  case('r2016_short')
     npar=size(r2016_short)+1
     if (npar > dim) goto 999
     params(1:npar-1)=r2016_short
     params(npar)=0.8d0

  case('r2017')
     npar=size(r2017)+1
     if (npar > dim) goto 999
     params(1:npar-1)=r2017
     params(npar)=1.d0

  case('heil17_short')
     npar=size(r2017_short)+1
     if (npar > dim) goto 999
     params(1:npar-1)=r2017_short
     params(npar)=0.8d0

  case('r2018')
     npar=size(r2018)+1
     if (npar > dim) goto 999
     params(1:npar-1)=r2018
     params(npar)=1.d0

  case('r2018_short')
     npar=size(r2018_short)+1
     if (npar > dim) goto 999
     params(1:npar-1)=r2018_short
     params(npar)=0.8d0

  case('r2022')
     npar=size(r2022)+1
     if (npar > dim) goto 999
     params(1:npar-1)=r2022
     params(npar)=1.d0

  case('qe8')
     npar=size(qe8)+1
     if (npar > dim) goto 999
     params(1:npar-1)=qe8
     params(npar)=1.d0

  case('qe8_short')
     npar=size(qe8_short)+1
     if (npar > dim) goto 999
     params(1:npar-1)=qe8_short
     params(npar)=0.8d0

  case('cvs-qe8')
     npar=size(cvs_qe8)+1
     if (npar > dim) goto 999
     params(1:npar-1)=cvs_qe8
     params(npar)=1.d0

  case default
     errmsg='Error in retrieve_hpar: unrecognised Hamiltonian name'
     call error_control

  end select

  return

999 continue
  errmsg='Error in retrieve_hpar: size(hpar) < npar'
  call error_control
  
end subroutine retrieve_hpar

