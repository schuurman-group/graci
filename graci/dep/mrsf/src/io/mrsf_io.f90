!**********************************************************************
! mrsf_io: C-string conversion, error handling and timing helpers
!**********************************************************************
module mrsf_io

  use mrsf_constants
  use iso_c_binding, only: c_char, c_null_char

  implicit none

contains

!######################################################################
! cstrlen: length of a NUL-terminated C string (excluding the NUL)
!######################################################################
  function cstrlen(cstr) result(n)

    character(kind=c_char), intent(in) :: cstr(*)
    integer(is)                        :: n

    n = 0
    do
       if (cstr(n+1) == c_null_char) exit
       n = n + 1
       if (n >= 255) exit
    enddo

  end function cstrlen

!######################################################################
! c2fstr: C string -> trimmed, left-adjusted Fortran string
!######################################################################
  subroutine c2fstr(cstr, fstr)

    character(kind=c_char), intent(in) :: cstr(*)
    character(len=*), intent(out)      :: fstr
    integer(is)                        :: i, n

    n = min(cstrlen(cstr), len(fstr))
    fstr = ''
    do i = 1, n
       fstr(i:i) = cstr(i)
    enddo
    fstr = adjustl(fstr)

  end subroutine c2fstr

!######################################################################
! mrsf_error: print an error message and stop
!######################################################################
  subroutine mrsf_error(msg)

    character(len=*), intent(in) :: msg

    write(6,'(/,2x,a,/)') 'MRSF error: '//trim(msg)
    flush(6)
    stop 1

  end subroutine mrsf_error

!######################################################################
! wall_time: wall-clock time in seconds
!######################################################################
  function wall_time() result(t)

    real(dp)    :: t
    integer(ib) :: count, rate

    call system_clock(count, rate)
    t = real(count, dp) / real(rate, dp)

  end function wall_time

!######################################################################
! freeunit: an unused Fortran unit number
!######################################################################
  subroutine freeunit(unit)

    integer(is), intent(out) :: unit
    logical                  :: lopen

    do unit = 20, 1000
       inquire(unit=unit, opened=lopen)
       if (.not. lopen) return
    enddo
    call mrsf_error('no free unit number')

  end subroutine freeunit

end module mrsf_io
