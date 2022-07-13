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
    character(len=13)    :: awall,acpu

    write(awall,'(F13.3)') twall
    write(acpu,'(F13.3)') tcpu
    
    write(6,'(/,x,a,1x,3a,1x,a)') 'Times (Wall/CPU) for '&
         //trim(routine)//':',trim(adjustl(awall)),'/',&
         trim(adjustl(acpu)),'seconds'
    
    return
    
  end subroutine report_times
  
!######################################################################
  
end module timing
