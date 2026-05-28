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
     procedure, public :: init_pyscf  => init_pyscf_base
     procedure, public :: h_1e        => h_1e_base
     procedure, public :: mo_ints     => mo_ints_base
     procedure, public :: mo_int      => mo_int_base
     procedure, public :: mo_ints_lr  => mo_ints_lr_base
     procedure, public :: mo_int_lr   => mo_int_lr_base
     procedure, public :: indx_ut     => indx_ut_base
     procedure, public :: finalize    => finalize_base
  end type eri

contains

  !
  ! Base initialisation routine: does nothing, but is required
  ! to set the interface that all other initialisation routines
  ! will use
  !
  subroutine init_pyscf_base(ints, core_file, eri_file, eri_lr_file)

    class(eri)                               :: ints
    character(len=255)                       :: core_file
    character(len=255)                       :: eri_file
    character(len=255), optional, intent(in) :: eri_lr_file

    return
  end subroutine init_pyscf_base

  !
  ! return core hamiltonian 1e integrals
  !
  function h_1e_base(ints, i, j) result(int_val)

    class(eri)              :: ints
    integer(is),intent(in)  :: i,j

    real(dp)                :: int_val

    int_val = ints%h_core(i,j)

    return
  end function h_1e_base

  !
  ! Return a list of 2e integrals (full ERI): stub, overridden by subtype
  !
  subroutine mo_ints_base(ints, indices, int_vec)

    class(eri)             :: ints
    integer(is),intent(in) :: indices(:,:)
    real(dp),intent(out)   :: int_vec(:)

    return
  end subroutine mo_ints_base

  !
  ! Return a single 2e integral (full ERI): stub, overridden by subtype
  !
  function mo_int_base(ints, i, j, k, l) result(int_val)

    class(eri)             :: ints
    integer(is),intent(in) :: i, j, k, l

    real(dp)               :: int_val

    return
  end function mo_int_base

  !
  ! Return a list of LR 2e integrals: stub returning zero (non-RSH / exact types)
  !
  subroutine mo_ints_lr_base(ints, indices, int_vec)

    class(eri)             :: ints
    integer(is),intent(in) :: indices(:,:)
    real(dp),intent(out)   :: int_vec(:)

    int_vec = 0.0d0

    return
  end subroutine mo_ints_lr_base

  !
  ! Return a single LR 2e integral: stub returning zero (non-RSH / exact types)
  !
  function mo_int_lr_base(ints, i, j, k, l) result(int_val)

    class(eri)             :: ints
    integer(is),intent(in) :: i, j, k, l

    real(dp)               :: int_val

    int_val = 0.0d0

    return
  end function mo_int_lr_base

  !
  ! deallocate data structures
  !
  subroutine finalize_base(ints)

    class(eri)             :: ints

  end subroutine finalize_base

  !
  ! returns upper-triangle index
  !
  function indx_ut_base(ints, i, j) result(ut)
    class(eri)                           :: ints
    integer(is), intent(in)              :: i
    integer(is), intent(in)              :: j

    integer(is), dimension(2)            :: ij
    integer(is)                          :: ut

    ij = (/ max(i,j), min(i,j) /)
    ut = int(ij(1)*(ij(1) - 1)/2 + ij(2))

    return
  end function indx_ut_base

end module integrals
