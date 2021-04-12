!######################################################################
! bitci_initialise: Interface to the bitci library. Sets all globally
!                   accessible variables that are needed to perform
!                   a CI calculation.
!######################################################################
#ifdef CBINDING
subroutine bitci_initialise(imult1,nel1,nmo1,mosym1,moen1,ipg1,enuc1,&
     iham,label1) bind(c,name="bitci_initialise")
#else
subroutine bitci_initialise(imult1,nel1,nmo1,mosym1,moen1,ipg1,enuc1,&
     iham,label1)
#endif

  use constants
  use bitglobal
  use utils
  use iomod
  use baseconf
  use csf
  use spin_coupling
  use precompute
  use hparam
  use iso_c_binding, only: C_CHAR
  
  implicit none

  integer(is), intent(in)            :: imult1,nel1,nmo1,ipg1
  integer(ib), intent(in)            :: mosym1(nmo1)
  integer(is), intent(in)            :: iham
  real(dp), intent(in)               :: moen1(nmo1)
  real(dp), intent(in)               :: enuc1
  real(dp)                           :: s,smax

#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: label1(*)
  character(len=255)                 :: label
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: label1
  character(len=255)                 :: label
#endif

!----------------------------------------------------------------------
! If C bindings are on, then convert the calculation label from the
! C char type to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(label1)
  call c2fstr(label1,label,length)
#else
  label=adjustl(trim(label1))
#endif

!----------------------------------------------------------------------
! Quick sanity check on the number of electrons and spin multiplicity
!----------------------------------------------------------------------
  if (mod(imult1,2) == 0 .and. mod(nel1,2) == 0 &
       .or. mod(imult1,2) /= 0.and. mod(nel1,2) /= 0) then
     write(errmsg,'(a,1x,i0,a,i0)') 'Nonsensical Nel/multiplicity:',&
          nel1,'/',imult1
     call error_control
  endif

!----------------------------------------------------------------------
! Exit if the given point group is not recognised
!----------------------------------------------------------------------
  if (ipg1.lt.1.or.ipg1.gt.8) then
     write(errmsg,'(a,1x,i0)') 'Illegal point group index:',ipg1
     call error_control
  endif
  
!----------------------------------------------------------------------
! Set the spin multiplicity
!----------------------------------------------------------------------
  imult=imult1

!----------------------------------------------------------------------
! Exit if the requested spin multiplicity is not supported by the
! maximum number of openshells
!----------------------------------------------------------------------
  smax=dble(nomax)/2.0d0
  s=dble(imult-1)/2.0
  if (s > smax) then
     write(errmsg,'(a)') 'The requested multiplicity is incompatible'&
          //' with the current value of nomax.'
     call error_control
  endif
  
!----------------------------------------------------------------------
! Set the number of electrons
!----------------------------------------------------------------------
  nel=nel1
  nel_beta=(nel-imult+1)/2
  nel_alpha=(nel-imult+1)/2+imult-1
  nel_spin(1)=nel_alpha
  nel_spin(2)=nel_beta
  
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
! Initialise the symmetry arrays
!----------------------------------------------------------------------
  call initialise_symmetry
  
!----------------------------------------------------------------------
! Point group
!----------------------------------------------------------------------
  ! Point group index
  ipg=ipg1

  ! Point group label
  pgroup=pglbls(ipg)

  ! Number of irreps
  nirrep=pgdim(ipg)

!----------------------------------------------------------------------
! Nuclear repulsion energy
!----------------------------------------------------------------------
  enuc=enuc1
  
!----------------------------------------------------------------------
! Set the bitstring integer array lengths
!----------------------------------------------------------------------
  n_int=(nmo-1)/64+1

!----------------------------------------------------------------------
! Load the Hamiltonian parameters
!----------------------------------------------------------------------
  call load_hpar(iham)

!----------------------------------------------------------------------
! Generate the CSFs for the given spin multiplicity up to the maximum
! number of open shells
!----------------------------------------------------------------------
  call generate_csfs

!----------------------------------------------------------------------
! Generate the spin coupling coefficients for the given spin
! multiplicity
!----------------------------------------------------------------------
  call generate_coupling_coefficients
  
!----------------------------------------------------------------------
! Generate the base determinant and base configuration
!----------------------------------------------------------------------
  call generate_base_det
  call generate_base_conf

!----------------------------------------------------------------------
! Pre-compute integrals
!----------------------------------------------------------------------
  call precompute_integrals
  
!----------------------------------------------------------------------
! Create the scratch directory
!----------------------------------------------------------------------
  ! Top level scratch directory
  scratchdir='bitscratch'
  call system('mkdir -p '//trim(scratchdir))

  ! Calculation specific scratch sub-directory
  scratchdir=trim(scratchdir)//'/'//trim(label)
  call system('mkdir -p '//trim(scratchdir))
  
!----------------------------------------------------------------------
! Initialise the scratch unit counter and allocate the scratch unit
! arrays
!----------------------------------------------------------------------
  nscratch=0
  maxunits=200
  allocate(scrunit(maxunits))
  scrunit=0
  allocate(scrname(maxunits))
  scrname=''

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
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
  pgdim(1)=1
  irreplbl(0,1)='A'
  
  ! Ci
  pgdim(2)=2
  irreplbl(0,2)='Ag'
  irreplbl(1,2)='Au'
  
  ! C2
  pgdim(3)=2
  irreplbl(0,3)='A'
  irreplbl(1,3)='B'
  
  ! Cs
  pgdim(4)=2
  irreplbl(0,4)='A'''
  irreplbl(1,4)='A'''''
  
  ! C2h
  pgdim(5)=4
  irreplbl(0,5)='Ag'
  irreplbl(1,5)='Bg'
  irreplbl(2,5)='Au'
  irreplbl(3,5)='Bu'
  
  ! C2v
  pgdim(6)=4
  irreplbl(0,6)='A1'
  irreplbl(1,6)='A2'
  irreplbl(2,6)='B1'
  irreplbl(3,6)='B2'
  
  ! D2
  pgdim(7)=4
  irreplbl(0,7)='A1'
  irreplbl(1,7)='B1'
  irreplbl(2,7)='B2'
  irreplbl(3,7)='B3'

  ! D2h
  pgdim(8)=8
  irreplbl(0,8)='A1g'
  irreplbl(1,8)='B1g'
  irreplbl(2,8)='B2g'
  irreplbl(3,8)='B3g'
  irreplbl(4,8)='A1u'
  irreplbl(5,8)='B1u'
  irreplbl(6,8)='B2u'
  irreplbl(7,8)='B3u'
  
  return
  
end subroutine initialise_symmetry
