module exactmod

  use constants
  use integrals
  
  implicit none

  !
  ! Exact integrals type
  !
  type, extends(eri) :: exact
   contains
     procedure   :: init => init_exact
  end type exact

contains

  !
  ! Eaxct integrals intitialisation routine
  !
  subroutine init_exact(ints,nipar,ipar,nrpar,rpar)

    class(exact)            :: ints
    integer(is), intent(in) :: nipar,nrpar
    integer(is), intent(in) :: ipar(nipar)
    real(dp), intent(in)    :: rpar(nrpar)

    ! Method name
    ints%name='Exact'

    ! No. MOs
    ints%nmo=ipar(1)
    
    return

  end subroutine init_exact

end module exactmod
