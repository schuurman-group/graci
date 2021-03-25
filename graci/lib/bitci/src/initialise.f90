!######################################################################
! bitci_initialise: interface to the bitci library. Sets all globally
!                   accessible variables that are needed to perform
!                   a CI calculation:
!
!                   - Spin multiplicity
!
!                   - Number of electrons
!
!                   - MO basis information
!
!                   - Symmetry information, including the point group
!
!                   - Integrals interface (?)
!
!                   - Number of 64-bit integers needed to represent
!                     a single Slater determinant
!######################################################################
#ifdef CBINDING
subroutine bitci_initialise(imult1,nel1,nmo1,mosym1,moen1,pgroup1) &
     bind(c,name="bitci_initialise")
#else
  subroutine bitci_initialise(imult1,nel1,nmo1,mosym1,moen1,pgroup1)
#endif

  use constants
  use bitglobal
  use utils
  use iomod
  
  implicit none

  integer(is), intent(in) :: imult1,nel1,nmo1
  integer(ib), intent(in) :: mosym1(nmo1) !< Symmetries of the MOs
  real(dp), intent(in)    :: moen1(nmo1)
  character(len=*)        :: pgroup1
  integer(is)             :: k
  
!----------------------------------------------------------------------
! Quick sanity check on the number of electrons and spin multiplicity
!----------------------------------------------------------------------
  if (mod(imult1,2).eq.0.and.mod(nel1,2).eq.0 &
       .or.mod(imult1,2).ne.0.and.mod(nel1,2).ne.0) then
     write(errmsg,'(a,1x,i0,a,i0)') 'Nonsensical Nel/multiplicity:',&
          nel1,'/',imult1
     call error_control
  endif
     
!----------------------------------------------------------------------
! Set the spin multiplicity
!----------------------------------------------------------------------
  imult=imult1

!----------------------------------------------------------------------
! Set the number of electrons
!----------------------------------------------------------------------
  nel=nel1
  nel_beta=(nel-imult+1)/2
  nel_alpha=(nel-imult+1)/2+imult-1

!----------------------------------------------------------------------
! Set the number of spatial orbitals
!----------------------------------------------------------------------
  nmo=nmo1

!----------------------------------------------------------------------
! Set the MO information
!----------------------------------------------------------------------
  ! MO irreps
  allocate(mosym(nmo))
  mosym=mosym1

  ! MO energies
  allocate(moen(nmo))
  moen=moen1

!----------------------------------------------------------------------
! Point group label
!----------------------------------------------------------------------
  pgroup=pgroup1
  
  ! Convert the point group label to lowercase
  call lowercase(pgroup)

!----------------------------------------------------------------------
! Initialise the symmetry arrays
!----------------------------------------------------------------------
  call initialise_symmetry
  
!----------------------------------------------------------------------
! Set the bitstring integer array lengths
!----------------------------------------------------------------------
  n_int=(nmo-1)/64+1

!----------------------------------------------------------------------
! Generate the base determinant
!----------------------------------------------------------------------
  call generate_base_det
  
!----------------------------------------------------------------------
! Create the scratch directory
!----------------------------------------------------------------------
  scratchdir='bitscratch'
  call system('mkdir -p '//trim(scratchdir))
  
  return
    
end subroutine bitci_initialise

!######################################################################
! initialise_symmetry: sets point group and irrep labels, and the the
!                      point group multiplication array
!######################################################################
subroutine initialise_symmetry

  use constants
  use bitglobal
  use iomod
  
  implicit none

  integer(is)      :: i
  character(len=3) :: all_labels(8)

!----------------------------------------------------------------------
! Set the point group labels
!----------------------------------------------------------------------
  pglbls=(/'c1 ','ci ','c2 ','cs ','c2h','c2v','d2 ','d2h'/)

!----------------------------------------------------------------------
! Set the irrep labels
!----------------------------------------------------------------------
  ! C1
  nirrep(1)=1
  irreplbl(0,1)='a1'
  
  ! Ci
  nirrep(2)=2
  irreplbl(0,2)='ag'
  irreplbl(1,2)='au'
  
  ! C2
  nirrep(3)=2
  irreplbl(0,3)='a'
  irreplbl(1,3)='b'
  
  ! Cs
  nirrep(4)=2
  irreplbl(0,4)='a'''
  irreplbl(1,4)='a'''''
  
  ! C2h
  nirrep(5)=4
  irreplbl(0,5)='ag'
  irreplbl(1,5)='bg'
  irreplbl(2,5)='au'
  irreplbl(3,5)='bu'
  
  ! C2v
  nirrep(6)=4
  irreplbl(0,6)='a1'
  irreplbl(1,6)='a2'
  irreplbl(2,6)='b1'
  irreplbl(3,6)='b2'
  
  ! D2
  nirrep(7)=4
  irreplbl(0,7)='a1'
  irreplbl(1,7)='b1'
  irreplbl(2,7)='b2'
  irreplbl(3,7)='b3'

  ! D2h
  nirrep(8)=8
  irreplbl(0,8)='a1g'
  irreplbl(1,8)='b1g'
  irreplbl(2,8)='b2g'
  irreplbl(3,8)='b3g'
  irreplbl(4,8)='a1u'
  irreplbl(5,8)='b1u'
  irreplbl(6,8)='b2u'
  irreplbl(7,8)='b3u'
  
!----------------------------------------------------------------------
! Make sure that the point group label is 'legal', i.e., that it
! corresponds to a real Abelian point group
!----------------------------------------------------------------------
  ! Check the given point group label against the Abelian point groups
  ipg=-1
  do i=1,8
     if (pgroup.eq.pglbls(i)) then
        ipg=i
        exit
     endif
  enddo

  ! Exit if the given point group is not recognised
  if (ipg.eq.-1) then
     errmsg='Unrecognised Abelian point group: '//trim(pgroup)
     call error_control
  endif
  
  return
  
end subroutine initialise_symmetry
