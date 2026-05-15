!**********************************************************************
! Thin wrapper around the MKL service routines for runtime BLAS/LAPACK
! thread-count control. Compiles to no-ops when the library is not
! built against MKL (USE_MKL undefined), keeping the caller code
! identical.
!
! The MKL service routines are accessed through C bindings rather than
! `use mkl_service`. The MKL Fortran module location varies between
! oneAPI releases and is not always added to the include path by the
! compiler driver, whereas the C entry points are always linked in by
! -mkl / -qmkl.
!**********************************************************************
module mkl_compat

  use constants
#ifdef USE_MKL
  use, intrinsic :: iso_c_binding, only: c_int
#endif

  implicit none

#ifdef USE_MKL
  interface
     function mkl_get_max_threads() bind(C,name='MKL_Get_Max_Threads') &
          result(nthreads)
       import :: c_int
       integer(c_int) :: nthreads
     end function mkl_get_max_threads

     subroutine mkl_set_num_threads(nthreads) &
          bind(C,name='MKL_Set_Num_Threads')
       import :: c_int
       integer(c_int), value :: nthreads
     end subroutine mkl_set_num_threads
  end interface
#endif

contains

!######################################################################
! save_and_set_blas_threads: record the current MKL thread count, set
!                            it to n, and return the previous value
!                            for later restoration. No-op when not
!                            built against MKL.
!######################################################################
  function save_and_set_blas_threads(n) result(saved)

    implicit none

    integer(is), intent(in) :: n
    integer(is)             :: saved

#ifdef USE_MKL
    saved = int(mkl_get_max_threads(),kind=is)
    call mkl_set_num_threads(int(n,kind=c_int))
#else
    saved = 1
#endif

    return

  end function save_and_set_blas_threads

!######################################################################
! restore_blas_threads: restore the MKL thread count previously
!                       captured by save_and_set_blas_threads. No-op
!                       when not built against MKL.
!######################################################################
  subroutine restore_blas_threads(saved)

    implicit none

    integer(is), intent(in) :: saved

#ifdef USE_MKL
    call mkl_set_num_threads(int(saved,kind=c_int))
#endif

    return

  end subroutine restore_blas_threads

!######################################################################

end module mkl_compat
