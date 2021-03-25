module iomod

  save

  !
  ! Error message to be written when terminating the program
  !
  character(len=200) :: errmsg

contains

!#######################################################################
! freeunit: determines the first free unit number
!#######################################################################
  subroutine freeunit(unit)

    use constants
      
    implicit none

    integer(is) :: unit,i
    logical     :: lopen
      
    !
    ! Find the first free unit
    ! N.B. Save the first 20 io units for standard files
    !
    do i=20,1000
       inquire(unit=i,opened=lopen)
       if (.not.lopen) then
          unit=i
          exit
       endif
    enddo

    return

  end subroutine freeunit

!#######################################################################
! error_control: writes the passed string to the screen and, if open,
!                the log file, then terminates the program
!#######################################################################
  subroutine error_control
      
    implicit none
      
    logical :: lopen

    !
    ! Write the error message to the screen and terminate
    ! the program
    !
    write(6,'(/,2x,a,/)') trim(errmsg)
    stop
    
  end subroutine error_control
  
!#######################################################################
! scratch_name: for a given string filename, returns the string
!               scratchdir/filename
!#######################################################################
  subroutine scratch_name(filename,path)

    use constants
    use bitglobal
    
    implicit none

    character(len=*) :: filename,path

    path=trim(scratchdir)//'/'//trim(filename)
    
    return
    
  end subroutine scratch_name
    
!#######################################################################
  
end module iomod
