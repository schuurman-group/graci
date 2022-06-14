!######################################################################
! override_hparam: overrides the default DFT/MRCI Hamiltonian
!                  parameter values with a user-specified set
!######################################################################
#ifdef CBINDING
subroutine override_hparam(npar,par) bind(c,name="override_hparam")
#else
subroutine override_hparam(npar,par)
#endif

  use constants
  use hparam
  use iomod
  
  implicit none

  ! User-specified Hamiltonian parameter values
  integer(is), intent(in) :: npar
  real(dp), intent(in)    :: par(npar)

!----------------------------------------------------------------------
! First, some sanity checks
!----------------------------------------------------------------------
  ! Have the DFT/MRCI parameter arrays been initialised
  if (.not. allocated(hpar)) then
     errmsg='Error in override_hparam: ' &
          //'load_hpar needs to be called first'
     call error_control
  endif

  ! Has the correct no. parameters been passed
  if (npar /= nhpar) then
     errmsg='Error in override_hparam: incorrect number of parameters'
     call error_control
  endif

!----------------------------------------------------------------------
! Set the new parameter values
!----------------------------------------------------------------------
  hpar=par
  
  return
  
end subroutine override_hparam

!######################################################################
! get_hparam: returns the DFT/MRCI Hamiltonian parameters in use
!######################################################################
#ifdef CBINDING
subroutine get_hparam(npar,par) bind(c,name="get_hparam")
#else
subroutine get_hparam(npar,par)
#endif

  use constants
  use hparam
  use iomod
  
  implicit none

  ! Hamiltonian parameter values
  integer(is), intent(in) :: npar
  real(dp), intent(out)   :: par(npar)

!----------------------------------------------------------------------
! First, some sanity checks
!----------------------------------------------------------------------
  ! Have the DFT/MRCI parameter arrays been initialised
  if (.not. allocated(hpar)) then
     errmsg='Error in get_hparam: load_hpar needs to be called first'
     call error_control
  endif

  ! Has the correct no. parameters been passed
  if (npar /= nhpar) then
     errmsg='Error in get_hparam: incorrect number of parameters'
     call error_control
  endif

!----------------------------------------------------------------------
! Copy over the Hamiltonian parameter values
!----------------------------------------------------------------------
  par=hpar

  return
  
end subroutine get_hparam