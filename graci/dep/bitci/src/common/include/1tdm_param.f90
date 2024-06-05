!**********************************************************************
! Parameters of the various DFT/MRCI 1-TDM damping functions
!**********************************************************************
module tdm_param

  use constants

  implicit none

  save

  ! Number of 1-TDM damping functions implemented
  integer(is), parameter :: ndamp = 2

  ! 1-TDM damping function labels
  character(len=20), parameter, dimension(ndamp) :: damplbl= &
       ['abinitio            ', &
       'qe8                  ']

  ! 1-TDM damping function integer label
  integer(is)            :: idamping

  ! 1-TDM damping function parameters
  real(dp), allocatable  :: damppar(:)

  ! Number of 1-TDM damping function parameters
  integer(is)            :: ndamppar

!----------------------------------------------------------------------
! QE8 damping function
!----------------------------------------------------------------------
  real(dp), parameter, dimension(3) :: qe8_damping= &
       [0.692173d0, & ! p1
       4.611269d0, &  ! p2
       8.0d0]         ! n

contains

!######################################################################
! load_hpar: loads the 1-TDM damping function parameters
!######################################################################
  subroutine load_damping(idamp)

    use constants
    use bitglobal
    use iomod

    implicit none

    ! 1-TDM damping function label
    integer(is), intent(in) :: idamp

!----------------------------------------------------------------------
! Set the 1-TDM damping function integer label
!----------------------------------------------------------------------
    idamping=idamp

!----------------------------------------------------------------------
! Load the 1-TDM damping function parameters
!----------------------------------------------------------------------
    select case(idamp)
    
    case(1)
       ! Ab initio 1-TDMs: do nothing
       ltdmdamp=.false.
       return

    case(2)
       ! QE8 damping function
       ltdmdamp=.true.
       ndamppar=3
       allocate(damppar(ndamppar))
       damppar=qe8_damping
       
    case default
       ! Unrecognised 1-TDM damping function
       write(errmsg,'(a,x,i0)') &
            'Error in load_damping: unrecognised damping '&
            //'function number',idamp
       call error_control
       
    end select
    
    return
    
  end subroutine load_damping

!######################################################################
  
end module tdm_param
