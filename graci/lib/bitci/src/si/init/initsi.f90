!######################################################################
! bitsi_initialise: Interface to the bitsi library. Sets all the
!                   globally accessible variables that are needed to
!                   perform a State Interaction calculation.
!######################################################################
#ifdef CBINDING
subroutine bitsi_intialise(imultB1,imultK1,nelB1,nelK1,nmo1,ipg1) &
     bind(c,name='bitsi_initialise')
#else
subroutine bitsi_intialise(imultB1,imultK1,nelB1,nelK1,nmo1,ipg1)
#endif

  use constants
  use bitglobal
  use setsym
  use csf
  use spin_coupling
  use iomod
  
  implicit none

  ! Bra and Ket spin multiplicities and numbers of electrons
  integer(is), intent(in) :: imultB1,imultK1,nelB1,nelK1

  ! Dimensions
  integer(is), intent(in) :: nmo1

  ! Point group index
  integer(is), intent(in) :: ipg1
  
  ! Calculation type
  logical                 :: samemult,soc,rdm,dyson
  
  ! Everything else
  real(dp)                :: s,smax
  logical                 :: verbose
  
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
! maximum number of openshells
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
  if (ipg1.lt.1.or.ipg1.gt.8) then
     write(errmsg,'(a,1x,i0)') 'Illegal point group index:',ipg1
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
  nmo=nmo1

!----------------------------------------------------------------------
! Set the bitstring integer array lengths
!----------------------------------------------------------------------
  n_int=(nmo-1)/64+1

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
! Is this a spin-orbit coupling calculation?
!----------------------------------------------------------------------
  soc=.false.
  rdm=.false.
  dyson=.false.
  
  if (abs(imultB-imultK) == 2) then

     ! SOC calculation
     soc=.true.
     samemult=.false.

  else if (abs(imultB-imultK) == 1) then

     ! Dyson orbital calculation
     dyson=.true.
     samemult=.false.
     
  else if (imultB == imultK) then

     ! 1-RDM or 1-TDM calculation
     rdm=.true.
     samemult=.true.

  else

     ! Non-sensical combination of bra and ket spin multiplicities
     write(errmsg,'(a,1x,i0,a,i0)') &
          'Non-sensical combination of bra/ket spin multiplicities:',&
          imultB,'/',imultK
     call error_control

  endif

!----------------------------------------------------------------------
! Generate the bra and ket CSFs
!----------------------------------------------------------------------
  ! 1-RDM or 1-TDM calculation: equal bra and ket spin multiplicities
  if (rdm) then
     verbose=.false.
     call generate_csfs(imultB,nocase2,ncsfs,ndets,maxcsf,maxdet,&
          csfcoe,detvec,verbose)
  endif
     
  ! SOC calculation: non-equal bra and ket spin multiplicities, equal
  ! numbers of electrons
  if (soc) then
     verbose=.true.
     errmsg='SOC calculations not yet supported in bitSI'
     call error_control
  endif

  ! Dyson orbital calculation: non-equal bra and ket spin
  ! multiplicities and non-equal numbers of electrons
  if (dyson) then
     verbose=.true.
     errmsg='Dyson orbital calculations not yet supported in bitSI'
     call error_control
  endif
  
!----------------------------------------------------------------------
! Compute the spin-coupling coefficients
!----------------------------------------------------------------------
  ! 1-RDM or 1-TDM calculation: equal bra and ket spin multiplicities
  if (rdm) then
     verbose=.false.
     call generate_coupling_coefficients(imultB,nocase1,nocase2,&
          maxcsf,maxdet,ncsfs,ndets,csfcoe,detvec,npattern1,&
          npattern2,nspincp,N1s,verbose,spincp,patternmap,offspincp)
  endif
     
  ! SOC calculation: non-equal bra and ket spin multiplicities, equal
  ! numbers of electrons
  if (soc) then
     verbose=.true.
     errmsg='SOC calculations not yet supported in bitSI'
     call error_control
  endif

  ! Dyson orbital calculation: non-equal bra and ket spin
  ! multiplicities and non-equal numbers of electrons
  if (dyson) then
     verbose=.true.
     errmsg='Dyson orbital calculations not yet supported in bitSI'
     call error_control
  endif

!----------------------------------------------------------------------
! Scratch directory: the top-level bitscratch directory is assumed to
! have already been created in a preceding bitci calculation
!----------------------------------------------------------------------
  scratchdir='bitscratch/interaction'
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
  
end subroutine bitsi_intialise
