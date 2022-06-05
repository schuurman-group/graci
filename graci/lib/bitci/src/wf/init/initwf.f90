!######################################################################
! bitwf_initialise: Interface to the bitwf library. Sets all the
!                   globally accessible variables that are needed to
!                   evaluate matrix elements <psi_m|O|psi'_n> in the
!                   Slater determinant representation.
!######################################################################
#ifdef CBINDING
subroutine bitwf_initialise(imultB1,imultK1,nelB1,nelK1,nmoB1,nmoK1,&
     smo1,ipg1,calctype_in) bind(c,name='bitwf_initialise')
#else
subroutine bitwf_initialise(imultB1,imultK1,nelB1,nelK1,nmoB1,nmoK1,&
     smo1,ipg1,calctype_in)
#endif

  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use setsym
  use csf
  use iomod

  implicit none

  ! Calculation type
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: calctype_in(*)
  character(len=255)                 :: calctype
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: calctype_in
  character(len=255)                 :: calctype
#endif
  
  ! Bra and Ket spin multiplicities and numbers of electrons
  integer(is), intent(in) :: imultB1,imultK1,nelB1,nelK1

  ! Number of bra and ket MOs
  integer(is), intent(in) :: nmoB1,nmoK1

  ! Bra-ket MO overlap matrix
  real(dp), intent(in)    :: smo1(nmoB1,nmoK1)

  ! Point group indice
  integer(is), intent(in) :: ipg1
  
  ! Everything else
  integer(is)             :: i,j
  real(dp)                :: s,smax
  logical                 :: verbose
  
!----------------------------------------------------------------------
! If C bindings are on, then convert the calculation type character
! string from the C char type to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(calctype_in)
  call c2fstr(calctype_in,calctype,length)
#else
  calctype=adjustl(trim(calctype_in))
#endif

!----------------------------------------------------------------------
! Quick sanity check on the numbers of electrons and spin
! multiplicities
!----------------------------------------------------------------------
  ! Bra
  if (mod(imultB1,2) == 0 .and. mod(nelB1,2) == 0 &
       .or. mod(imultB1,2) /= 0.and. mod(nelB1,2) /= 0) then
     write(errmsg,'(a,1x,i0,a,i0)') &
          'Nonsensical bra Nel/multiplicity:',nelB1,'/',imultB1
     call error_control
  endif

  ! Ket
  if (mod(imultK1,2) == 0 .and. mod(nelK1,2) == 0 &
       .or. mod(imultK1,2) /= 0.and. mod(nelK1,2) /= 0) then
     write(errmsg,'(a,1x,i0,a,i0)') &
          'Nonsensical ket Nel/multiplicity:',nelK1,'/',imultK1
     call error_control
  endif

!----------------------------------------------------------------------
! Exit if the requested spin multiplicity is not supported by the
! maximum number of open shells
!----------------------------------------------------------------------
  ! S_max
  smax=dble(nomax)/2.0d0

  ! Bra
  s=dble(imultB1-1)/2.0
  if (s > smax) then
     write(errmsg,'(a)') &
          'The requested bra multiplicity is incompatible'&
          //' with the current value of nomax.'
     call error_control
  endif

  ! Ket
  s=dble(imultK1-1)/2.0
  if (s > smax) then
     write(errmsg,'(a)') &
          'The requested ket multiplicity is incompatible'&
          //' with the current value of nomax.'
     call error_control
  endif

!----------------------------------------------------------------------
! Exit if the given point group is not recognised
!----------------------------------------------------------------------
  if (ipg1 < 1 .or. ipg1 > 8) then
     write(errmsg,'(a,1x,i0)') 'Illegal bra point group index:',ipg1
     call error_control
  endif

!----------------------------------------------------------------------
! Set the bra and ket spin multiplicities
!----------------------------------------------------------------------
  imultB=imultB1
  imultK=imultK1

!----------------------------------------------------------------------
! Set the bra and ket numbers of electrons
!----------------------------------------------------------------------
  nelB=nelB1
  nelB_beta=(nelB-imultB+1)/2
  nelB_alpha=(nelB-imultB+1)/2+imultB-1
  nelB_spin(1)=nelB_alpha
  nelB_spin(2)=nelB_beta

  nelK=nelK1
  nelK_beta=(nelK-imultK+1)/2
  nelK_alpha=(nelK-imultK+1)/2+imultK-1
  nelK_spin(1)=nelK_alpha
  nelK_spin(2)=nelK_beta

!----------------------------------------------------------------------
! Set the number of spatial orbitals
!----------------------------------------------------------------------
  nmoB=nmoB1
  nmoK=nmoK1

!----------------------------------------------------------------------
! Set the bit string integer array lengths
!----------------------------------------------------------------------
  n_intB=(nmoB-1)/n_bits+1
  n_intK=(nmoK-1)/n_bits+1

!----------------------------------------------------------------------
! Save the MO overlap matrix
!----------------------------------------------------------------------
  allocate(smo(nmoB,nmoK))
  smo=smo1
  
!----------------------------------------------------------------------
! Initialise the symmetry arrays
!----------------------------------------------------------------------
  call initialise_symmetry
  
!----------------------------------------------------------------------
! Point group: note that we only currently support equal bra and ket
! point groups
!----------------------------------------------------------------------
  ! Point group index
  ipg=ipg1

  ! Point group label
  pgroup=pglbls(ipg)

  ! Number of irreps
  nirrep=pgdim(ipg)
  
!----------------------------------------------------------------------
! Generate the bra and ket CSFs
!----------------------------------------------------------------------
  ! Don't spam the user with info that's already been generated
  verbose=.false.

  ! Bra CSFs
  call generate_csfs(imultB,nocase2,ncsfsB,ndetsB,maxcsfB,&
       maxdetB,csfcoeB,detvecB,verbose)

  ! Ket CSFs
  call generate_csfs(imultK,nocase2,ncsfsK,ndetsK,maxcsfK,&
       maxdetK,csfcoeK,detvecK,verbose)

!----------------------------------------------------------------------
! Scratch directory: the top-level bitscratch directory is assumed to
! have already been created in a preceding bitci calculation
!----------------------------------------------------------------------
  scratchdir='bitscratch/wf'
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
  
end subroutine bitwf_initialise
