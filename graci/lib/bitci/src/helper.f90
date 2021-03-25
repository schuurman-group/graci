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
! retrieve_cas_dets: returns the CAS determinants
!######################################################################
#ifdef CBINDING
subroutine retrieve_cas_dets(d,ncas) bind(c,name="retrieve_cas_dets")
#else
subroutine retrieve_cas_dets(d,ncas)
#endif
  
  use constants
  use bitglobal
  use iomod
  
  integer(is), intent(in) :: ncas
  integer(ib)             :: d(n_int,2,ncas)
  integer(is)             :: nchk
  integer(is)             :: iscratch
  character(len=60)       :: casfile
    
  !
  ! Open the scratch file
  !
  call freeunit(iscratch)
  call scratch_name('detcas',casfile)
  open(iscratch,file=casfile,form='unformatted',status='unknown')

  !
  ! Read the number of CAS determinants
  !
  read(iscratch) nchk

  !
  ! Exit if the number of determinants written to disk is different
  ! to the size of the determinant array that was passed to us
  !
  if (nchk.ne.ncas) then
     errmsg='Inconsistent numbers of determinants in retrieve_cas_dets'
     call error_control
  endif

  !
  ! Read the CAS determinants
  !
  read(iscratch) d
  
  !
  ! Close the scratch file
  !
  close(iscratch)
    
  return

end subroutine retrieve_cas_dets

!######################################################################
! retrieve_ras_dets: returns the RAS determinants
!######################################################################
#ifdef CBINDING
subroutine retrieve_ras_dets(d,nras) bind(c,name="retrieve_ras_dets")
#else
subroutine retrieve_ras_dets(d,nras)
#endif
  
  use constants
  use bitglobal
  use iomod
  
  integer(is), intent(in) :: nras
  integer(ib)             :: d(n_int,2,nras)
  integer(is)             :: nchk
  integer(is)             :: iscratch
  character(len=60)       :: rasfile
    
  !
  ! Open the scratch file
  !
  call freeunit(iscratch)
  call scratch_name('detras',rasfile)
  open(iscratch,file=rasfile,form='unformatted',status='unknown')

  !
  ! Read the number of RAS determinants
  !
  read(iscratch) nchk

  !
  ! Exit if the number of determinants written to disk is different
  ! to the size of the determinant array that was passed to us
  !
  if (nchk.ne.nras) then
     errmsg='Inconsistent numbers of determinants in retrieve_ras_dets'
     call error_control
  endif

  !
  ! Read the RAS determinants
  !
  read(iscratch) d
  
  !
  ! Close the scratch file
  !
  close(iscratch)
    
  return

end subroutine retrieve_ras_dets

!######################################################################
! retrieve_mrci_dets: returns the MRCI determinants
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
  integer(is)              :: bufsize,nrecords
  integer(ib), allocatable :: buffer(:,:,:)
  integer(is)              :: nchk
  integer(is)              :: detcnt,ndets,i
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
  read(iscratch) bufsize
  read(iscratch) nrecords

  !
  ! Allocate the buffer
  !
  allocate(buffer(n_int,2,bufsize))

  !
  ! Read in the determinants
  !
  detcnt=0
  do i=1,nrecords
     read(iscratch) ndets,buffer
     d(:,:,detcnt+1:detcnt+ndets)=buffer(:,:,1:ndets)
     detcnt=detcnt+ndets
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

!######################################################################
! write_det_string: writes the determinant d to the character string
!                   string in the format 2...0...u...d...
!######################################################################
#ifdef CBINDING
subroutine write_det_string(d,string) bind(c,name="write_det_string")
#else
subroutine write_det_string(d,string)
#endif
  
  use constants
  use bitglobal
  use bitutils
  
  implicit none

  integer(ib), intent(in)       :: d(n_int,2)
  character(len=*), intent(out) :: string
  integer(is)                   :: ihomo,imo,i,k,occa,occb,occ
  
  !
  ! Initialisation
  !
  string=''

  !
  ! Index of the highest occupied orbital
  !
  ihomo=homo_index(d)

  !
  ! Write the determinant character string
  !
  do imo=1,ihomo

     ! Block index
     k=(imo-1)/64+1

     ! Orbital position with the block
     i=imo-1-(k-1)*64

     ! Alpha-spin occupancy
     if (btest(d(k,1),i)) then
        occa=1
     else
        occa=0
     endif

     ! Beta-spin occupancy
     if (btest(d(k,2),i)) then
        occb=1
     else
        occb=0
     endif

     ! Total occupancy
     occ=occa+occb

     ! Write the total occupancy to the character string
     if (occ.eq.2) then
        write(string(imo:imo),'(a1)') '2'
     else if (occ.eq.0) then
        write(string(imo:imo),'(a1)') '0'
     else
        if (occa.eq.1) then
           write(string(imo:imo),'(a1)') 'u'
        else
           write(string(imo:imo),'(a1)') 'd'
        endif
     endif
        
  enddo

  return
  
end subroutine write_det_string
