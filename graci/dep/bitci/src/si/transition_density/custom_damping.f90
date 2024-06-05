!######################################################################
! override_damping_param: overrides the default 1-TDM damping function
!                         parameter values with a user-specified set
!######################################################################
#ifdef CBINDING
subroutine override_damping_param(npar,par) &
     bind(c,name="override_damping_param")
#else
subroutine override_damping_param(npar,par)
#endif

  use constants
  use tdm_param
  use iomod
  
  implicit none

  ! User-specified 1-TDM damping function parameter values
  integer(is), intent(in) :: npar
  real(dp), intent(in)    :: par(npar)

!----------------------------------------------------------------------
! First, some sanity checks
!----------------------------------------------------------------------
  ! Has the 1-TDM damping function parameter array been initialised?
  if (.not. allocated(damppar)) then
     errmsg='Error in override_damping_param: ' &
          //'load_damping needs to be called first'
     call error_control
  endif

  ! Has the correct no. parameters been passed (including desel)
  if (npar /= ndamppar) then
     errmsg='Error in override_damping_param: '&
          //'incorrect number of parameters'
     call error_control
  endif

!----------------------------------------------------------------------
! Set the new parameter values
!----------------------------------------------------------------------
  damppar=par
  
  return

end subroutine override_damping_param
