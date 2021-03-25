!
!
! M. S. Schuurman, Jan 6., 2020
!
!
module accuracy 
  implicit none

  public
    integer, parameter :: sik         = selected_int_kind(4)       ! Small integers
    integer, parameter :: ik          = selected_int_kind(8)       ! 64-bit integers
    integer, parameter :: hik         = selected_int_kind(18)      ! "Pointer" integers - sufficient to mem address
    integer, parameter :: drk         = kind(1.d0)                 ! Double-precision real kind
    integer, parameter :: rk          = selected_real_kind(14,17)  ! Standard reals and complex
    integer, parameter :: ark         = selected_real_kind(14,17)  ! Possible increased precision reals and complex

    complex(drk),parameter        :: zero_c   = (0., 0.)
    complex(drk),parameter        :: one_c    = (1., 0.)
    complex(drk),parameter        :: I_drk    = (0., 1.)
    real(drk),parameter           :: zero_drk = 0.
    real(drk),parameter           :: one_drk  = 1.
    integer(ik),parameter         :: zero_ik  = 0
    integer(ik),parameter         :: one_ik   = 1
    real(drk),parameter           :: mp_drk   = 1.e-16

  !
 
end module accuracy 

