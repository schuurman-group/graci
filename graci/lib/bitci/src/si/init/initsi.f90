!######################################################################
! bitsi_initialise: Interface to the bitsi library. Sets all the
!                   globally accessible variables that are needed to
!                   perform a State Interaction calculation.
!######################################################################
#ifdef CBINDING
subroutine bitsi_intialise(imultB1,imultK1,nelB1,nelK1,nmo1,ipg1,&
     calctype_in) bind(c,name='bitsi_initialise')
#else
subroutine bitsi_intialise(imultB1,imultK1,nelB1,nelK1,nmo1,ipg1,&
     calctype_in)
#endif

  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use setsym
  use csf
  use spin_coupling
  use spin_coupling_k1
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

  ! Dimensions
  integer(is), intent(in) :: nmo1

  ! Point group index
  integer(is), intent(in) :: ipg1
  
  ! Calculation type
  logical                 :: samemult
  
  ! Everything else
  real(dp)                :: s,smax,sb,sk
  real(dp)                :: tiny=1e-10_dp
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
! Make sure that the calculation type is recognised
!----------------------------------------------------------------------
  select case(calctype)

  case('tdm')
     continue

  case('soc')
     continue

  case default
     errmsg='Error in bitsi_intialise: unrecognised calculation type:'&
          //' '//trim(calctype)
     call error_control

  end select
        
!----------------------------------------------------------------------
! Sanity checks on the bra and ket spin multiplicities
!----------------------------------------------------------------------
  select case(calctype)

  case('tdm')
     ! Calculation of matrix elements over the singlet excitation
     ! operators: used in the calculation of 1-TDMs

     ! Bra and ket spin multiplicities have to be equal
     if (imultB /= imultK) then
        errmsg='Error in bitsi_intialise: same bra and ket '&
             //'multiplicities required for a 1-TDM calculation'
        call error_control
     endif
     
  case('soc')
     ! Calculation of matrix elements over the triplet excitation
     ! operators: used in the calculation of SOC matrix elemens

     ! SOC calculations don't make sense unless Delta S = 0, +1
     ! because we are only considering the matrix elements
     ! over T_pq^(1,k), k=0,+1 (the k=-1 terms can be generated
     ! from the k=+1 ones)
     sb=dble(imultB-1)/2.0
     sk=dble(imultk-1)/2.0
     if (sb-sk > 1.0d0+tiny .or. sb-sk < 0.0d0-tiny) then
        errmsg='Error in bitsi_intialise: |S_bra-S_ket| = 0 or 1'
        call error_control
     endif
     
  end select
     
!----------------------------------------------------------------------
! Generate the bra and ket CSFs
!----------------------------------------------------------------------
  select case(calctype)

  case('tdm')
     ! 1-TDM calculation: equal bra and ket spin multiplicities
     verbose=.false.
     call generate_csfs(imultB,nocase2,ncsfs,ndets,maxcsf,maxdet,&
          csfcoe,detvec,verbose)

  case('soc')
     ! SOC calculation: potentially non-equal bra and ket spin
     !                  multiplicities
     if (imultB /= imultK) then
        verbose=.true.
     else
        verbose=.false.
     endif
     ! Bra CSFs
     call generate_csfs(imultB,nocase2,ncsfsB,ndetsB,maxcsfB,maxdetB,&
          csfcoeB,detvecB,verbose)
     ! Ket CSFs
     call generate_csfs(imultK,nocase2,ncsfsK,ndetsK,maxcsfK,maxdetK,&
          csfcoeK,detvecK,verbose)

  end select
  
!----------------------------------------------------------------------
! Compute the spin-coupling coefficients
!----------------------------------------------------------------------
  select case(calctype)

  case('tdm')
     ! 1-TDM calculation: calculation of spin-coupling coefficients
     !                    over the singlet excitation operators E_pq
     verbose=.false.
     call generate_coupling_coefficients(imultB,nocase1,nocase2,&
          maxcsf,maxdet,ncsfs,ndets,csfcoe,detvec,nspincp,N1s,&
          verbose,spincp,patternmap,offspincp)

  case('soc')
     ! SOC calculation: calculation of spin-coupling coefficients for
     !                  the triplet excitation operators
     !                  T_pq^(1,k), k = 0 or +1
     verbose=.true.
     if (imultB == imultK) then
        ! Same spin multiplicities: compute spin-coupling coefficients
        ! for T_pq^(1,k=0)
        errmsg='The k=0 SCC code needs writing'
        call error_control
     else
        ! Different spin multiplicities: compute spin-coupling
        ! coefficients for T_pq^(1,k=+1)
        call scc_k1(imultB,imultK,&
             nocase1,nocase2,&
             maxcsfB,maxdetB,ncsfsB,ndetsB,csfcoeB,detvecB,&
             maxcsfK,maxdetK,ncsfsK,ndetsK,csfcoeK,detvecK,&
             nspincp,N1s,verbose,spincp,patternmap1,patternmap2a,&
             patternmap2b,offspincp)
     endif
     
  end select
  
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
