!**********************************************************************
! Parameters of the various DFT/MRCI 1-TDM damping functions
!**********************************************************************
module tdm_param

  use constants

  implicit none

  save

  ! Number of 1-TDM damping functions implemented
  integer(is), parameter :: ntdm = 2

  ! 1-TDM damping function labels
  character(len=20), parameter, dimension(ntdm) :: tdmlbl= &
       ['abinitio            ', &
       'qe8                  ']

  ! 1-TDM damping function integer label
  integer(is)            :: itdm

  ! 1-TDM damping function parameters
  real(dp), allocatable  :: tdmpar(:)

  ! Number of 1-TDM damping function parameters
  integer(is)            :: ntdmpar

!----------------------------------------------------------------------
! QE8 damping function
!----------------------------------------------------------------------
  real(dp), parameter, dimension(5) :: qe8= &
       [0.692173d0, & ! p1
       4.611269d0, &  ! p2
       8.0d0]         ! n

contains

!######################################################################
! load_hpar: loads the 1-TDM damping function parameters
!######################################################################
  subroutine load_tdm(idamp)

    use constants
    use bitglobal
    use iomod

    implicit none

    ! 1-TDM damping function label
    integer(is), intent(in) :: idamp

!----------------------------------------------------------------------
! Set the 1-TDM damping function integer label
!----------------------------------------------------------------------
    itdm=idamp

!----------------------------------------------------------------------
! Load the 1-TDM damping function parameters
!----------------------------------------------------------------------

    print*,''
    print*,'itdm:',itdm
    print*,''
    stop
    
    return
    
  end subroutine load_tdm

!######################################################################
  
end module tdm_param
