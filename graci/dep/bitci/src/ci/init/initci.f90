!######################################################################
! bitci_initialise: Interface to the bitci library. Sets all globally
!                   accessible variables that are needed to perform
!                   a CI calculation.
!######################################################################
#ifdef CBINDING
subroutine bitci_initialise(imult1,nel1,nmo1,mosym1,moen1,ipg1,&
     escf1,ham1,label1,verbose1) bind(c,name="bitci_initialise")
#else
subroutine bitci_initialise(imult1,nel1,nmo1,mosym1,moen1,ipg1,&
     escf1,ham1,label1,verbose1)
#endif

  use constants
  use bitglobal
  use utils
  use setsym
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
  real(dp), intent(in)               :: moen1(nmo1)
  real(dp), intent(in)               :: escf1
  logical, intent(in)                :: verbose1
  real(dp)                           :: s,smax

#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: label1(*),ham1(*)
  character(len=255)                 :: label,ham
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: label1,ham1
  character(len=255)                 :: label
#endif

  integer(is)                        :: i,iham(1)
  
!----------------------------------------------------------------------
! If C bindings are on, then convert the Hamiltonian and calculation
! labels from the C char type to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(label1)
  call c2fstr(label1,label,length)
  length=cstrlen(ham1)
  call c2fstr(ham1,ham,length)
#else
  label=adjustl(trim(label1))
  ham=adjustl(trim(ham1))
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
! Exit if the requested spin multiplicity is not supported by the
! maximum number of open shells
!----------------------------------------------------------------------
  smax=dble(nomax)/2.0d0
  s=dble(imult1-1)/2.0
  if (s > smax) then
     write(errmsg,'(a)') 'The requested multiplicity is incompatible'&
          //' with the current value of nomax.'
     call error_control
  endif
  
!----------------------------------------------------------------------
! Exit if the given point group is not recognised
!----------------------------------------------------------------------
  if (ipg1 < 1 .or. ipg1 > 8) then
     write(errmsg,'(a,1x,i0)') 'Illegal point group index:',ipg1
     call error_control
  endif

!----------------------------------------------------------------------
! Make sure that the requested Hamiltonian is supported
!----------------------------------------------------------------------
  ! Hamiltonian index
  iham=findloc(hlbl,value=trim(ham))

  ! Exit if the Hamiltonian label was not found
  if (iham(1) == 0) then
     write(errmsg,'(a,1x,a)') 'Unrecognised Hamiltonian:',trim(ham)
     call error_control
  endif

!----------------------------------------------------------------------
! Set the verbosity logical flag
!----------------------------------------------------------------------
  verbose=verbose1
  
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
! E_SCF: note that this is the SCF energy if this is an MRCI
!        calculation but is the DFT energy if we are performing a
!        DFT/MRCI calculation
!----------------------------------------------------------------------
  escf=escf1

!----------------------------------------------------------------------
! Determine the index of the last non-frozen virtual MO
!----------------------------------------------------------------------
  lastvirt=nmo
  do i=1,nmo
     if (moen(i) >= efreeze) then
        lastvirt=i-1
        exit
     endif
  enddo
  
!----------------------------------------------------------------------
! Set the bitstring integer array lengths
!----------------------------------------------------------------------
  n_int=(nmo-1)/n_bits+1

!----------------------------------------------------------------------
! Load the Hamiltonian parameters
!----------------------------------------------------------------------  
  call load_hpar(iham(1))

!----------------------------------------------------------------------
! Generate the CSFs for the given spin multiplicity up to the maximum
! number of open shells
!----------------------------------------------------------------------
  call generate_csfs(imult,nocase2,ncsfs,ndets,maxcsf,maxdet,&
       csfcoe,detvec,verbose)

!----------------------------------------------------------------------
! Generate the spin coupling coefficients for the given spin
! multiplicity
!----------------------------------------------------------------------
  call generate_coupling_coefficients(imult1,nocase1,nocase2,maxcsf,&
       maxdet,ncsfs,ndets,csfcoe,detvec,nspincp,N1s,verbose,spincp,&
       patmap,offspincp)

!----------------------------------------------------------------------
! Generate the base determinant and base configuration
!----------------------------------------------------------------------
  call generate_base_det

  call generate_base_conf

!----------------------------------------------------------------------
! Pre-compute integrals
!----------------------------------------------------------------------
! Note that this has to be called *after* load_hpar in order to know
! whether the full Fock matrix has to be computed
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
! bitci_int_initialize: Initialises the bitci integrals class object
!######################################################################
#ifdef CBINDING
subroutine bitci_int_initialize(integral_src, integral_method, &
     integral_precision, hcore_file, eri_file) bind(c,name="bitci_int_initialize")
#else
subroutine bitci_int_initialize(integral_src, integral_method, &
     integral_precision, hcore_file, eri_file)
#endif

  use constants
  use bitglobal
  use iomod
  use exactmod
  use dfmod

#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: integral_src(*)
  character(kind=C_CHAR), intent(in) :: integral_method(*)
  character(kind=C_CHAR), intent(in) :: integral_precision(*)
  character(kind=C_CHAR), intent(in) :: hcore_file(*)
  character(kind=C_CHAR), intent(in) :: eri_file(*)

  character(len=255)                 :: int_src
  character(len=255)                 :: int_method
  character(len=255)                 :: int_precision
  character(len=255)                 :: f_core
  character(len=255)                 :: f_eri

  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: integral_src
  character(len=*), intent(in)       :: integral_method
  character(len=*), intent(in)       :: integral_precision
  character(len=*), intent(in)       :: hcore_file
  character(len=*), intent(in)       :: eri_file

  character(len=255)                 :: int_src
  character(len=255)                 :: int_method
  character(len=255)                 :: int_precision
  character(len=255)                 :: hcore
  character(len=255)                 :: eri
#endif

!----------------------------------------------------------------------
! If C bindings are on, then convert the calculation label from the
! C char type to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(integral_src)
  call c2fstr(integral_src, int_src,length)
  length=cstrlen(integral_method)
  call c2fstr(integral_method, int_method,length)
  length=cstrlen(integral_precision)
  call c2fstr(integral_precision, int_precision, length)
  length=cstrlen(hcore_file)
  call c2fstr(hcore_file, f_core,length)
  length=cstrlen(eri_file)
  call c2fstr(eri_file, f_eri,length)
#else
  int_src       = adjustl(trim(integral_src))
  int_method    = adjustl(trim(integral_method))
  int_precision = adjustl(trim(integral_precision))
  f_core        = adjustl(trim(hcore_file))
  f_eri         = adjustl(trim(eri_file))
#endif

!----------------------------------------------------------------------
! first identify what kind of integral method we're using:
! exact, or density df
!----------------------------------------------------------------------
  if (trim(adjustl(int_method)) .eq. 'exact') then
     select case(trim(adjustl(int_precision)))
         case ('single')
           allocate(exact_sp::bitci_ints)
         case ('double')
           allocate(exact_dp::bitci_ints)
         case default
           allocate(exact_dp::bitci_ints)
     end select

  else if (trim(adjustl(int_method)) .eq. 'df') then
     select case(trim(adjustl(int_precision)))
         case ('single')
           print *,'allocating sp df...'
           allocate(df_sp::bitci_ints)
           print *,'done.'
         case ('double')
           allocate(df_dp::bitci_ints)
         case default
           allocate(df_dp::bitci_ints)
     end select

  else
     stop 'integral method not recognized in bitci_init_integrals: '&
          //trim(adjustl(int_method))

  endif

!----------------------------------------------------------------------
! second: initialize integral object using the appropriate
!         interface
!----------------------------------------------------------------------
  if (trim(adjustl(int_src)) .eq. 'pyscf') then
     call bitci_ints%init_pyscf(f_core, f_eri)
  else
     stop 'integral interface not recognized in bitci_init_integrals: '&
          //trim(adjustl(int_src))
  endif

  return
  
end subroutine bitci_int_initialize
