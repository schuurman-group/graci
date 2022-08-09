!**********************************************************************
! Calculation of 1-TDMs for MRCI wavefunctions
!**********************************************************************

!######################################################################
! transition_density_mrci: Controls the calculation of MRCI 1-TDMs
!######################################################################
#ifdef CBINDING
subroutine transition_density_mrci(irrepB,irrepK,nrootsB,nrootsK,&
     npairs,iroots,rhoij,conffileB_in,vecfileB_in,conffileK_in,&
     vecfileK_in) bind(c,name='transition_density_mrci')
#else
subroutine transition_density_mrci(irrepB,irrepK,nrootsB,nrootsK,&
     npairs,iroots,rhoij,conffileB_in,vecfileB_in,conffileK_in,&
     vecfileK_in)
#endif

  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use iomod
  use conftype
  use merge
  use tdm
    
  implicit none

  ! MRCI configuration and eigenvector file names
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: conffileB_in(*),vecfileB_in(*),&
                                        conffileK_in(*),vecfileK_in(*)
  character(len=255)                 :: conffileB,vecfileB,&
                                        conffileK,vecfileK
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: conffileB_in,vecfileB_in,&
                                        conffileK_in,vecfileK_in
  character(len=255)                 :: conffileB,vecfileB,&
                                        conffileK,vecfileK
#endif

  ! Irreps and no. roots
  integer(is), intent(in)  :: irrepB,irrepK,nrootsB,nrootsK

  ! Indices of the pairs of states for which 1-TDMs are requested
  integer(is), intent(in)  :: npairs
  integer(is), intent(in)  :: iroots(npairs,2)
  
  ! Transition density matrix
  real(dp), intent(out)    :: rhoij(nmo,nmo,npairs)

  ! MRCI configuration derived types
  type(mrcfg)              :: cfgB,cfgK

  ! Scratch file numbers
  integer(is)              :: confscrB,vecscrB,confscrK,vecscrK

  ! Eigenpairs
  real(dp), allocatable    :: vecB(:,:),enerB(:),vecK(:,:),enerK(:)
  
  ! Everything else
  integer(is)              :: i,j,k
  integer(is)              :: nvecB,nvecK
  integer(is), allocatable :: iBra(:),iKet(:)
  integer(is), allocatable :: Bmap(:),Kmap(:)
  integer(is), allocatable :: ireadB(:),ireadK(:)
  
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  if (verbose) then
     write(6,'(/,52a)') ('-',i=1,52)
     write(6,'(3(x,a))') &
          'Transition density matrix calculation for the',&
          trim(irreplbl(irrepB,ipg)),' &',trim(irreplbl(irrepK,ipg)),&
          'subspaces'
     write(6,'(52a)') ('-',i=1,52)
  endif
     
!----------------------------------------------------------------------
! If C bindings are on, then convert the MRCI configuration and
! eigenvector file names from the C char type to the Fortran character
! type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(conffileB_in)
  call c2fstr(conffileB_in,conffileB,length)
  length=cstrlen(vecfileB_in)
  call c2fstr(vecfileB_in,vecfileB,length)
  length=cstrlen(conffileK_in)
  call c2fstr(conffileK_in,conffileK,length)
  length=cstrlen(vecfileK_in)
  call c2fstr(vecfileK_in,vecfileK,length)
#else
  conffileB=adjustl(trim(conffileB_in))
  vecfileB=adjustl(trim(conffileB_in))
  conffileK=adjustl(trim(conffileK_in))
  vecfileK=adjustl(trim(vecfileK_in))
#endif
  
!----------------------------------------------------------------------
! Register the configuration and eigenvector scratch files
!----------------------------------------------------------------------
  call register_scratch_file(confscrB,conffileB)
  call register_scratch_file(vecscrB,vecfileB)

  call register_scratch_file(confscrK,conffileK)
  call register_scratch_file(vecscrK,vecfileK)

!----------------------------------------------------------------------
! Set up the bra and ket configuration derived types
!----------------------------------------------------------------------
  call cfgB%initialise(irrepB,confscrB)
  call cfgK%initialise(irrepK,confscrK)

  if (verbose) then
     write(6,'(/,x,a,x,i0)') 'Bra CSF basis dimension:',cfgB%csfdim
     write(6,'(x,a,x,i0)') 'Ket CSF basis dimension:',cfgK%csfdim
  endif
  
!----------------------------------------------------------------------
! Merge the bra and ket reference spaces
!----------------------------------------------------------------------
  call merge_ref_space(cfgB,cfgK)

!----------------------------------------------------------------------
! Which eigenvectors are needed?
!----------------------------------------------------------------------
  !
  ! Bra and ket states appearing in the requested 1-TDMs
  !
  allocate(iBra(nrootsB), iKet(nrootsK))
  iBra=0; iKet=0
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
  allocate(Bmap(npairs), Kmap(npairs))
  Bmap=0; Kmap=0
  do i=1,npairs
     Bmap(i)=sum(iBra(1:iroots(i,1)))
     Kmap(i)=sum(iKet(1:iroots(i,2)))
  enddo

!----------------------------------------------------------------------
! Read in the bra eigenvectors
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(vecB(cfgB%csfdim,nvecB), enerB(nvecB))
  vecB=0.0d0; enerB=0.0d0

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
  call read_some_eigenpairs(vecscrB,vecB,enerB,cfgB%csfdim,nvecB,ireadB)
  
!----------------------------------------------------------------------
! Read in the ket eigenvectors
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(vecK(cfgK%csfdim,nvecK), enerK(nvecK))
  vecK=0.0d0; enerK=0.0d0

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
  call read_some_eigenpairs(vecscrK,vecK,enerK,cfgK%csfdim,nvecK,ireadK)
  
!----------------------------------------------------------------------
! Compute the 1-TDMs
!----------------------------------------------------------------------
  call tdm_mrci(cfgB,cfgK,cfgB%csfdim,cfgK%csfdim,nvecB,nvecK,vecB,&
       vecK,npairs,rhoij,Bmap,Kmap)
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  call cfgB%finalise
  call cfgK%finalise
  deallocate(vecB)
  deallocate(vecK)
  deallocate(enerB)
  deallocate(enerK)
  deallocate(iBra)
  deallocate(iKet)
  deallocate(Bmap)
  deallocate(Kmap)
  deallocate(ireadB)
  deallocate(ireadK)
  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
    
end subroutine transition_density_mrci
