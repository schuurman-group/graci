module integrals

  use constants
  
  implicit none

  !
  ! Base integrals type
  !
  type eri
     integer(is)           :: nmo           ! number of MOs
     real(dp), allocatable :: h_core(:,:)

   contains
     procedure, public :: init_pyscf => init_pyscf
     procedure, public :: h_1e       => h_1e
     procedure, public :: mo_ints    => mo_ints
     procedure, public :: mo_int     => mo_int
     procedure, public :: indx_ut    => indx_ut
     procedure, public :: finalize   => finalize 
  end type eri

contains

  !
  ! Base initialisation routine: does nothing, but is required
  ! to set the interface that all other initialisation routines
  ! will use
  !
  subroutine init_pyscf(ints, core_file, eri_file)

    class(eri)              :: ints
    character(len=255)      :: core_file
    character(len=255)      :: eri_file

    return
  end subroutine init_pyscf
  
  ! 
  ! return core hamiltonian 1e integrals
  !
  function h_1e(ints, i, j) result(int_val)

    class(eri)              :: ints
    integer(is),intent(in)  :: i,j

    real(dp)                :: int_val

    int_val = ints%h_core(i,j)

    return
  end function h_1e

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  subroutine mo_ints(ints, indices, int_vec)

    class(eri)             :: ints
    integer(is),intent(in) :: indices(:,:)
    real(dp), intent(out)  :: int_vec(:)

    return
  end subroutine mo_ints
  
  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  function mo_int(ints, i, j, k, l) result(int_val)

    class(eri)             :: ints
    integer(is),intent(in) :: i, j, k, l

    real(dp)               :: int_val

    return
  end function mo_int

  !
  ! returns upper-triangle index
  !
  function indx_ut(ints, i, j) result(ut)
    class(eri)                           :: ints
    integer(is), intent(in)              :: i
    integer(is), intent(in)              :: j

    integer(is), dimension(2)            :: ij
    integer(is)                          :: ut

    ij = (/ max(i,j), min(i,j) /)
    ut = int(ij(1)*(ij(1) - 1)/2 + ij(2))

    return
  end function indx_ut

  !
  ! deallocate data structures
  !
  subroutine finalize(ints)

    class(eri)             :: ints

  end subroutine finalize

end module integrals
