module integrals

  use constants
  
  implicit none

  !
  ! Base integrals type
  !
  type eri
     private
     integer(is)       :: nmo
     character(len=60) :: name
   contains
     procedure, public :: init => init_base
  end type eri

contains

  !
  ! Base initialisation routine: does nothing, but is required
  ! to set the interface that all other initialisation routines
  ! will use
  !
  ! * Michael: we can change this to be whatever we think is best
  !            maybe just passing the name of a h5 file will be the
  !            way to do this? *
  !
  subroutine init_base(ints,nipar,ipar,nrpar,rpar)

    class(eri)              :: ints
    integer(is), intent(in) :: nipar,nrpar
    integer(is), intent(in) :: ipar(nipar)
    real(dp), intent(in)    :: rpar(nrpar)
        
    return
    
  end subroutine init_base
    
end module integrals
