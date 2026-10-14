!**********************************************************************
! scf_io: C-string conversion, timing and the MKL thread interface
!**********************************************************************
module scf_io

  use scf_constants
  use iso_c_binding

  implicit none

#ifdef USE_MKL
  interface
     function mkl_threads_local(nt) bind(c, name='MKL_Set_Num_Threads_Local') result(prev)
       import :: c_int
       integer(c_int), value :: nt
       integer(c_int)        :: prev
     end function mkl_threads_local
  end interface
#endif

  logical :: verbose = .false.

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

!######################################################################
! blas_sequential / blas_restore: sequential MKL inside OpenMP regions
!######################################################################
  subroutine blas_sequential()
#ifdef USE_MKL
    integer(c_int) :: prev
    prev = mkl_threads_local(1_c_int)
#endif
  end subroutine blas_sequential

  subroutine blas_restore()
#ifdef USE_MKL
    integer(c_int) :: prev
    prev = mkl_threads_local(0_c_int)
#endif
  end subroutine blas_restore

end module scf_io
