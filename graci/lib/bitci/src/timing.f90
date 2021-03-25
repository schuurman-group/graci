module timing

  implicit none

contains

!######################################################################
! get_times: returns both the wall and CPU time
!######################################################################
  subroutine get_times(twall,tcpu)

    use constants

    implicit none
    
    real(dp), intent(out) :: twall,tcpu
    integer(is)           :: c,cr,cm

    ! cpu time
    call cpu_time(tcpu)
    
    ! wall time
    call system_clock(c,cr,cm)
    twall=dble(c)/dble(cr)
        
    return
    
  end subroutine get_times
  
!######################################################################

  subroutine report_times(twall,tcpu,routine)

    use constants
    
    implicit none

    real(dp), intent(in) :: twall,tcpu
    character(len=*)     :: routine
    
    write(6,'(/,a,1x,i0,a,i0,1x,a)') 'Times (Wall/CPU) for '&
         //trim(routine)//':',nint(twall),'/',nint(tcpu),'seconds'
    
    return
    
  end subroutine report_times
  
!######################################################################
  
end module timing
