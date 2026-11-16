!**********************************************************************
! grid_io: C-string conversion and timing
!**********************************************************************
module grid_io

  use grid_constants
  use iso_c_binding

  implicit none

contains

!######################################################################
! c2fstr: NUL-terminated C string -> left-adjusted Fortran string
!######################################################################
  subroutine c2fstr(cstr, fstr)

    character(kind=c_char), intent(in) :: cstr(*)
    character(len=*), intent(out)      :: fstr
    integer(is)                        :: i

    fstr = ''
    do i = 1, len(fstr)
       if (cstr(i) == c_null_char) exit
       fstr(i:i) = cstr(i)
    enddo
    fstr = adjustl(fstr)

  end subroutine c2fstr

!######################################################################
! wall_time: wall-clock time in seconds
!######################################################################
  function wall_time() result(t)

    real(dp)    :: t
    integer(ib) :: count, rate

    call system_clock(count, rate)
    t = real(count, dp) / real(rate, dp)

  end function wall_time

end module grid_io
