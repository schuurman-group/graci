!**********************************************************************
! Interface to the liboverlap library
!**********************************************************************

!######################################################################
! detoverlap: Calculation of wave function overlaps using the external
!             liboverlap library.
!             Takes as input a list of bra-ket state pairs and bitwf
!             wave function scratch file numbers, reads in the
!             required determinant bit strings, etc. and makes the
!             call to liboverlap.
!######################################################################
#ifdef CBINDING
subroutine detoverlap(irrep,nrootsB,nrootsK,npairs,iroots,wfscrB,&
     wfscrK,Sij) bind(c,name='detoverlap')
#else
subroutine detoverlap(irrep,nrootsB,nrootsK,npairs,iroots,wfscrB,&
     wfscrK,Sij)
#endif

  use constants
  use bitglobal
  use iowf
  use iomod

  implicit none

  ! Irrep and no. roots
  integer(is), intent(in)  :: irrep,nrootsB,nrootsK

  ! Indices of the pairs of states for which overlaps are requested
  integer(is), intent(in)  :: npairs
  integer(is), intent(in)  :: iroots(npairs,2)

  ! Scratch file numbers
  integer(is), intent(in)  :: wfscrB,wfscrK

  ! Wave function overlaps
  real(dp), intent(out)    :: Sij(npairs,2)

  ! Determinant bit strings
  integer(ib), allocatable :: detB(:,:,:),detK(:,:,:)

  ! Eigenvectors
  real(dp), allocatable    :: vecB(:,:),vecK(:,:)
  
  ! Everything else
  integer(is)              :: i,k
  integer(is)              :: ndetB,ndetK
  integer(is)              :: nvecB,nvecK
  integer(is), allocatable :: iBra(:),iKet(:)
  integer(is), allocatable :: Bmap(:),Kmap(:)
  integer(is), allocatable :: ireadB(:),ireadK(:)
  
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  write(6,'(/,52a)') ('-',i=1,52)
  write(6,'(2(x,a))') &
       'Wave function overlap calculation for the',&
       trim(irreplbl(irrep,ipg)),'subspace'
  write(6,'(52a)') ('-',i=1,52)

!----------------------------------------------------------------------
! Get the number of bra and ket determinants
!----------------------------------------------------------------------
  call read_ndet(wfscrB,ndetB)
  call read_ndet(wfscrK,ndetK)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(detB(n_intB,2,ndetB))
  detB=0_ib
  
  allocate(detK(n_intK,2,ndetB))
  detK=0_ib

  allocate(iBra(nrootsB), iKet(nrootsK))
  iBra=0; iKet=0

  allocate(Bmap(npairs), Kmap(npairs))
  Bmap=0; Kmap=0
  
!----------------------------------------------------------------------
! Which eigenvectors are needed?
!----------------------------------------------------------------------
  !
  ! Bra and ket states appearing in the requested 1-TDMs
  !
  do i=1,npairs
     iBra(iroots(i,1))=1
     iKet(iroots(i,2))=1
  enddo

  !
  ! Number of bra and ket eigenvectors
  !
  nvecB=sum(iBra)
  nvecK=sum(iKet)

  !
  ! Bra-ket pair to eigenvector mapping
  !
  ! Bmap(n) <-> index of the bra eigenvector needed to evaluate the
  !             n'th 1-TDM
  ! Kmap(n) <-> index of the Ket eigenvector needed to evaluate the
  !             n'th 1-TDM
  !
  do i=1,npairs
     Bmap(i)=sum(iBra(1:iroots(i,1)))
     Kmap(i)=sum(iKet(1:iroots(i,2)))
  enddo

!----------------------------------------------------------------------
! Read in the bra eigenvectors and determinant bit strings
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(vecB(ndetB,nvecB))
  vecB=0.0d0

  allocate(ireadB(nvecB))
  ireadB=0
  
  ! List of needed eigenvectors
  k=0
  do i=1,nrootsB
     if (iBra(i) == 1) then
        k=k+1
        ireadB(k)=i
     endif
  enddo

  ! Read in the eigenvectors  
  call read_detwf(wfscrB,ndetB,nvecB,n_intB,detB,vecB,ireadB)

!----------------------------------------------------------------------
! Read in the ket eigenvectors and determinant bit strings
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(vecK(ndetK,nvecK))
  vecK=0.0d0

  allocate(ireadK(nvecK))
  ireadK=0
  
  ! List of needed eigenvectors
  k=0
  do i=1,nrootsK
     if (iKet(i) == 1) then
        k=k+1
        ireadK(k)=i
     endif
  enddo

  ! Read in the eigenvectors  
  call read_detwf(wfscrK,ndetK,nvecK,n_intK,detK,vecK,ireadK)

!----------------------------------------------------------------------
! Call to liboverlap
!----------------------------------------------------------------------
  call overlap()

  return
  
end subroutine detoverlap
