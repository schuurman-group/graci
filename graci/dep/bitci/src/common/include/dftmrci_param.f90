!**********************************************************************
! Parameters of the various DFT/MRCI Hamiltonians
!**********************************************************************
module hparam

  use constants
  
  implicit none

  save

  ! Number of Hamiltonians implemented
  integer(is), parameter :: nham=13
  
  ! Hamiltonian labels
  character(len=20), parameter, dimension(nham) :: hlbl= &
       ['abinitio            ', &
        'grimme              ', &
        'grimme_short        ', &
        'r2016               ', &
        'r2016_short         ', &
        'r2017               ', &
        'r2017_short         ', &
        'r2018               ', &
        'r2018_short         ', &
        'r2022               ', &
        'qe8                 ', &
        'qe8_short           ', &
        'cvs-qe8             ']

  ! Hamiltonian integer label
  integer(is)           :: ihamiltonian

  ! Hamiltonian parameter array
  real(dp), allocatable :: hpar(:)

  ! Number of Hamiltonian parameters
  integer(is)           :: nhpar

  ! Configuration selection energy cutoff
  real(dp)              :: desel
  
!----------------------------------------------------------------------
! Grimme's original DFT/MRCI Hamiltonian
! J. Chem. Phys., 111, 5645 (1999)
!----------------------------------------------------------------------
  ! Singlet parameters, delta E_sel = 1.0
  real(dp), parameter, dimension(5) :: grimme1= &
       [0.6195d0, & ! p1
       3.2719d0, &  ! p2
       0.5102d0, &  ! pJ
       0.5945d0, &  ! p[0]
       0.1058d0]    ! alpha
  
  ! Singlet parameters, delta E_sel = 0.8
  real(dp), parameter, dimension(5) :: grimme1_short= &
       [0.6290d0, & ! p1
        8.0000d0, & ! p2
        0.5030d0, & ! pJ
        0.6110d0, & ! p[0]
        0.1190d0]   ! alpha
  
  ! Triplet parameters, delta E_sel = 1.0
  real(dp), parameter, dimension(5) :: grimme3= &
       [0.6195d0, & ! p1
       3.2719d0, &  ! p2
       0.4930d0, &  ! pJ
       0.0000d0, &  ! p[0]
       0.0563d0]    ! alpha

  ! Triplet parameters, delta E_sel = 0.8
  real(dp), parameter, dimension(5) :: grimme3_short= &
       [0.6290d0, & ! p1
       8.0000d0, &  ! p2
       0.4860d0, &  ! pJ
       0.0000d0, &  ! p[0]
       0.0630d0]    ! alpha

!----------------------------------------------------------------------
! Lyskov's 2016 redesigned DFT/MRCI Hamiltonian
! J. Chem. Phys., 144, 034104 (2016)
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(4) :: r2016= &
       [0.507894d0, & ! pJ
       0.355895d0, &  ! pF
       0.568168d0, &  ! p1
       1.788d0]       ! p2
       
  ! delta E_sel = 0.8
  real(dp), parameter, dimension(4) :: r2016_short= &
       [0.503506d0, & ! pJ
       0.368122d0, &  ! pF
       0.579809d0, &  ! p1
       2.187d0]       ! p2

!----------------------------------------------------------------------
! Heil's 2017 DFT/MRCI Hamiltonian for odd and even electron numbers
! J. Chem. Phys., 147, 194104 (2017)
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(4) :: r2017= &
       [0.503001d0, & ! pJ
       0.358727d0, &  ! pF
       0.563893d0, &  ! p1
       1.8571d0]      ! p2

  ! delta E_sel = 0.8
  real(dp), parameter, dimension(4) :: r2017_short= &
       [0.500779d0, & ! pJ
       0.356986d0, &  ! pF
       0.573523d0, &  ! p1
       1.9266d0]      ! p2

!----------------------------------------------------------------------
! Heil's DFT/MRCI Hamiltonian for transition metal complexes
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(4) :: r2018= &
       [0.508918d0, & ! pJ
       0.362362d0, &  ! pF
       0.558411d0, &  ! p1
       4.47165d0]     ! p2

  ! delta E_sel = 0.8
  real(dp), parameter, dimension(4) :: r2018_short= &
       [0.505808d0, & ! pJ
       0.359626d0, &  ! pF
       0.577732d0, &  ! p1
       11.499113d0]   ! p2

!----------------------------------------------------------------------
! R2022 DFT/MRCI Hamiltonian
! J. Phys. Chem A, 127, 2011 (2023)
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(5) :: r2022= &
       [3.4673d0, & ! p2
       0.5085d0,  & ! pJ^he
       0.4649d0,  & ! pJ^hhee
       0.3426d0,  & ! px^he
       0.5416d0]    ! px^hhee
  
!----------------------------------------------------------------------
! QE8 Hamiltonians
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(5) :: qe8= &
       [0.425623d0, & ! pJ
       0.252259d0, &  ! pF
       0.692173d0, &  ! p1
       4.611269d0, &  ! p2
       8.0d0]         ! n

  ! delta E_sel = 0.8
  real(dp), parameter, dimension(5) :: qe8_short= &
       [0.419d0, & ! pJ
       0.258d0,  & ! pF
       0.712d0,  & ! p1
       4.69d0,   & ! p2
       8.0d0]      ! n
  
  real(dp), parameter, dimension(7) :: cvs_qe8= &
       [0.425623d0, & ! pJ_vv
       0.252259d0, &  ! pF_vv
       0.499646d0, &  ! p1
       0.214962d0, &  ! p2
       8.0d0, &       ! n
       0.560644d0, &  ! pJ_cv
       0.252259d0]    ! pF_cv
  
contains

!######################################################################
! load_hpar: loads the DFT/MRCI Hamiltonian parameters
!######################################################################
  subroutine load_hpar(iham)

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
       ! Canonical: do nothing
       ldftmrci=.false.
       desel=999.9d0
       return
       
    case(2)
       ! Grimme, standard
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       if (imult == 1) then
          hpar=grimme1
          desel=1.0d0
       else if (imult == 3) then
          hpar=grimme3
          desel=1.0d0
       else
          errmsg='Only singlet and triplet states are supported using' &
               //' Grimme''s Hamiltonian'
          call error_control
       endif
       
    case(3)
       ! Grimme, short
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       if (imult == 1) then
          hpar=grimme1_short
          desel=0.8d0
       else if (imult == 3) then
          hpar=grimme3_short
          desel=0.8d0
       else
          errmsg='Only singlet and triplet states are supported using' &
               //' Grimme''s Hamiltonian'
          call error_control
       endif
       
    case(4)
       ! Lyskov, standard
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2016
       desel=1.0d0
       
    case(5)
       ! Lyskov, short
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2016_short
       desel=0.8d0
       
    case(6)
       ! Heil17, standard
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2017
       desel=1.0d0
       
    case(7)
       ! Heil17, short
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2017_short
       desel=0.8d0
       
    case(8)
       ! Heil18, standard
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2018
       desel=1.0d0
              
    case(9)
       ! Heil18, short
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=r2018_short
       desel=0.8d0

    case(10)
       ! R2022
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       hpar=r2022
       desel=1.0d0
       
    case(11)
       ! QE8
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       hpar=qe8
       desel=1.0d0

    case(12)
       ! QE8, short
       ldftmrci=.true.
       nhpar=5
       allocate(hpar(nhpar))
       hpar=qe8_short
       desel=0.8d0
       
    case(13)
       ! CVS-QE8
       ldftmrci=.true.
       nhpar=7
       allocate(hpar(nhpar))
       hpar=cvs_qe8
       desel=1.0d0

    case default
       ! Unrecognised Hamiltonian
       write(errmsg,'(a,x,i0)') &
            'Error in load_hpar: unrecognised Hamiltonian number',&
            iham
       call error_control
       
    end select
       
    return
    
  end subroutine load_hpar

!######################################################################
  
end module hparam
