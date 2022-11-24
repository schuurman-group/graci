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
       ['canonical           ', &
        'grimme_standard     ', &
        'grimme_short        ', &
        'lyskov_standard     ', &
        'lyskov_short        ', &
        'heil17_standard     ', &
        'heil17_short        ', &
        'heil18_standard     ', &
        'heil18_short        ', &
        'cvs_standard        ', &
        'cvs_short           ', &
        'test_heil17        ', &
        'cvs_test_heil17    ']

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
  real(dp), parameter, dimension(5) :: grimme1_standard= &
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
  real(dp), parameter, dimension(5) :: grimme3_standard= &
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
  real(dp), parameter, dimension(4) :: lyskov_standard= &
       [0.507894d0, & ! pJ
       0.355895d0, &  ! pF
       0.568168d0, &  ! p1
       1.788d0]       ! p2
       
  ! delta E_sel = 0.8
  real(dp), parameter, dimension(4) :: lyskov_short= &
       [0.503506d0, & ! pJ
       0.368122d0, &  ! pF
       0.579809d0, &  ! p1
       2.187d0]       ! p2

!----------------------------------------------------------------------
! Heil's 2017 DFT/MRCI Hamiltonian for odd and even electron numbers
! J. Chem. Phys., 147, 194104 (2017)
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(4) :: heil17_standard= &
       [0.503001d0, & ! pJ
       0.358727d0, &  ! pF
       0.563893d0, &  ! p1
       1.8571d0]      ! p2

  ! delta E_sel = 0.8
  real(dp), parameter, dimension(4) :: heil17_short= &
       [0.500779d0, & ! pJ
       0.356986d0, &  ! pF
       0.573523d0, &  ! p1
       1.9266d0]      ! p2

!----------------------------------------------------------------------
! Heil's DFT/MRCI Hamiltonian for transition metal complexes
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(4) :: heil18_standard= &
       [0.508918d0, & ! pJ
       0.362362d0, &  ! pF
       0.558411d0, &  ! p1
       4.47165d0]     ! p2

  ! delta E_sel = 0.8
  real(dp), parameter, dimension(4) :: heil18_short= &
       [0.505808d0, & ! pJ
       0.359626d0, &  ! pF
       0.577732d0, &  ! p1
       11.499113d0]   ! p2

!----------------------------------------------------------------------
! Experimental CVS-DFT/MRCI Hamiltonian for K-edge core-excited
! states
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(6) :: cvs_standard= &
       [0.503001d0, & ! pJ^(vv)
       0.358727d0, &  ! pF^(vv)
       0.563893d0, &  ! p1
       1.8571d0, &    ! p2
       0.503001d0, &  ! pJ^(cv)
       0.358727d0]    ! pF^(cv)

  ! delta E_sel = 0.8
  real(dp), parameter, dimension(6) :: cvs_short= &
       [0.500779d0, & ! pJ
       0.356986d0, &  ! pF
       0.573523d0, &  ! p1
       1.9266d0, &    ! p2
       0.500779d0, &  ! pJ^(cv)
       0.356986d0]    ! pF^(cv)
 
!----------------------------------------------------------------------
! Experimental Hamiltonians testing new functionals. These are the 
! standard Heil17 parameterization and the experimental CVS 
! parameterization for K-edge core-excited states
!----------------------------------------------------------------------
  ! delta E_sel = 1.0
  real(dp), parameter, dimension(4) :: test_heil17= &
       [0.503001d0, & ! pJ
       0.358727d0, &  ! pF
       0.563893d0, &  ! p1
       1.8571d0]      ! p2

  ! delta E_sel = 1.0
  real(dp), parameter, dimension(6) :: cvs_test_heil17= &
       [0.503001d0, & ! pJ^(vv)
       0.358727d0, &  ! pF^(vv)
       0.563893d0, &  ! p1
       1.8571d0, &    ! p2
       0.503001d0, &  ! pJ^(cv)
       0.358727d0]    ! pF^(cv)

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
          hpar=grimme1_standard
          desel=1.0d0
       else if (imult == 3) then
          hpar=grimme1_standard
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
          hpar=grimme1_short
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
       hpar=lyskov_standard
       desel=1.0d0
       
    case(5)
       ! Lyskov, short
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=lyskov_short
       desel=0.8d0
       
    case(6)
       ! Heil17, standard
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=heil17_standard
       desel=1.5d0
       
    case(7)
       ! Heil17, short
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=heil17_short
       desel=0.8d0
       
    case(8)
       ! Heil18, standard
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=heil18_standard
       desel=1.0d0
              
    case(9)
       ! Heil18, short
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=heil18_short
       desel=0.8d0

    case(10)
       ! CVS, standard
       ldftmrci=.true.
       nhpar=6
       allocate(hpar(nhpar))
       hpar=cvs_standard
       desel=1.0d0

    case(11)
       ! CVS, short
       ldftmrci=.true.
       nhpar=6
       allocate(hpar(nhpar))
       hpar=cvs_short
       desel=0.8d0

    case(12)
       ! test standard
       ldftmrci=.true.
       nhpar=4
       allocate(hpar(nhpar))
       hpar=test_heil17
       desel=1.0d0

    case(13)
       ! CVS test, standard
       ldftmrci=.true.
       nhpar=6
       allocate(hpar(nhpar))
       hpar=cvs_test_heil17
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
