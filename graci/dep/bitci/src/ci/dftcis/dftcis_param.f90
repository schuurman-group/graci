!**********************************************************************
! DFT/CIS Hamiltonian parameters
!**********************************************************************
module dftcis_param

  use constants
  
  implicit none

  save

  ! Number of Hamiltonians implemented
  integer(is), parameter :: nham = 2

  ! Hamiltonian labels
  character(len=20), parameter, dimension(nham) :: hlbl= &
       ['grimme1_b3lyp',&
       'ottawa1_bhlyp']

  ! Hamiltonian integer label
  integer(is)            :: ihamiltonian

  ! Hamiltonian parameter array
  real(dp), allocatable  :: hpar(:)

  ! Number of Hamiltonian parameters
  integer(is)            :: nhpar

!----------------------------------------------------------------------
! Grimme's original DFT/CIS Hamiltonian for singlet states and
! the B3LYP functional
! Chem. Phys. Lett., 259, 128 (1996)
!----------------------------------------------------------------------
  real(dp), parameter, dimension(3) :: grimme1_b3lyp= &
       [0.317d0,  & ! c1
       0.0033d0,  & ! c2
       1.27e+7_dp]  ! c3

!----------------------------------------------------------------------
! Our parameterisation for singlet states and the BHLYP functional
!----------------------------------------------------------------------
  !***************************
  !* Temporary parameter set *
  !* These need optimising   *
  !***************************
  real(dp), parameter, dimension(3) :: ottawa1_bhlyp= &
       [0.59596d0, & ! c1
       0.0033d0,   & ! c2
       1.27e+7_dp]   ! c3
  
contains
  
!######################################################################
! load_hpar_dftcis: loads the DFT/CIS Hamiltonian parameters
!######################################################################
  subroutine load_hpar_dftcis(iham)

    use constants
    use bitglobal
    use iomod

    implicit none

    ! Hamiltonian integer label
    integer(is), intent(in) :: iham

!----------------------------------------------------------------------
! Set the Hamiltonian integer label
!----------------------------------------------------------------------
    ihamiltonian=iham

!----------------------------------------------------------------------
! Load the Hamiltonian parameters
!----------------------------------------------------------------------
    select case(iham)

    case(1)
       ! Grimme's original parameterisation
       ihamiltonian=1
       nhpar=3
       allocate(hpar(nhpar))
       hpar=grimme1_b3lyp

    case(2)
       ! Our singlet/BHLYP parameterisation
       ihamiltonian=2
       nhpar=3
       allocate(hpar(nhpar))
       hpar=ottawa1_bhlyp
       
    case default
       ! Unrecognised Hamiltonian
       write(errmsg,'(a,x,i0)') &
            'Error in load_hpar_dftcis: '&
            //'unrecognised Hamiltonian number',iham
       call error_control
       
    end select
    
    return
    
  end subroutine load_hpar_dftcis
  
!######################################################################
! unload_hpar_dftcis: Deallocation of saved allocated arrays
!######################################################################
  subroutine unload_hpar_dftcis

    use constants
    
    implicit none

    nhpar=0
    
    deallocate(hpar)
    
    return
    
  end subroutine unload_hpar_dftcis
    
!######################################################################
  
end module dftcis_param
