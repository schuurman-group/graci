module bitglobal

  use constants
  
  implicit none

  save

!----------------------------------------------------------------------
! Common variables
!----------------------------------------------------------------------
  !
  ! Name of the scratch directory
  !
  character(len=255) :: scratchdir

  !
  ! Scratch files
  !
  integer(is)                     :: maxunits   ! Maximum number of scratch
                                                ! units
  integer(is)                     :: nscratch   ! Number of scratch units in use
  integer(is), allocatable        :: scrunit(:) ! Values of the scratch units
  character(len=255), allocatable :: scrname(:) ! Names of the scratch units

  !
  ! Maximum number of open shells
  !
  integer(is), parameter :: nexmax=8        ! Max. excitation degree relative
                                            ! to the base configuration allowed
                                            ! in the CSF basis
  integer(is), parameter :: nomax=10        ! Max. no. open shells allowed in
                                            ! the CSF basis
  integer(is), parameter :: nocase1=nomax   ! Max. number of open shells
                                            ! for the Case 1 spin coupling
                                            ! coefficients
  integer(is), parameter :: nocase2=nomax+2 ! Max. number of open shells
                                            ! for the Case 2 spin coupling
                                            ! coefficients

  !
  ! Spin coupling coefficients
  !
  integer(is)              :: nspincp(2)      ! Number of unique spin coupling
                                              ! coefficients 
  integer(is)              :: maxpattern(2)   ! Maximum possible simplified
                                              ! spatial occupation pattern values
  integer(is)              :: npattern1,&     ! Number of Case 1 and Case 2
                              npattern2       ! patterns
  integer(is), allocatable :: patternmap1(:)  ! Case 1 pattern -> array index
                                              ! mapping
  integer(is), allocatable :: patternmap2(:)  ! Case 2 pattern -> array index
                                              ! mapping
  integer(ib), allocatable :: N1s(:)          ! Bit strings comprised of N 1's

  real(dp), allocatable    :: spincp1(:,:,:)  ! Case 1 spin coupling coefficients
  real(dp), allocatable    :: spincp2(:,:,:)  ! Case 2 spin coupling coefficients
  
  !
  ! Buffered I/O
  !
  integer(is), parameter :: bufsize=2097152 ! Buffer size
                                            ! (16 MB of 8-byte reals)
  
!----------------------------------------------------------------------
! Configuration interaction variables
!----------------------------------------------------------------------
  !
  ! Dimensions
  !
  integer(is) :: nmo,nel,nel_alpha,nel_beta,nel_spin(2)
  
  !
  ! Base determinant and configuration
  !
  integer(ib), allocatable :: det0(:,:),conf0(:,:)
  
  !
  ! MO information
  !
  integer(ib), allocatable :: mosym(:)
  real(dp), allocatable    :: moen(:)

  !
  ! Nuclear repulsion energy
  !
  real(dp) :: enuc
  
  !
  ! Symmetry
  !
  integer(is)      :: ipg             ! Index of the point group
  integer(is)      :: nirrep          ! Number of irreps in the point group
  integer(is)      :: pgdim(8)        ! Dimensions of every Abelian point group
  character(len=3) :: pgroup          ! Point group label
  character(len=3) :: irreplbl(0:7,8) ! Labels of the irreps in each Abelian
                                      ! point group
  character(len=3) :: pglbls(8)       ! Labels of every Abelian point group

  !
  ! Spin multiplicity
  !
  integer(is) :: imult
  
  !
  ! CSF information
  !
  integer(is)              :: ncsfs(0:nocase2) ! Number of CSFs as a function of
                                               ! the number of open shells
  integer(is)              :: maxcsf           ! Maximum number of CSFs across
                                               ! all numbers of open shells
  real(dp), allocatable    :: csfcoe(:,:,:)    ! CSF expansion coefficients
  integer(is)              :: ndets(0:nocase2) ! Number of determinants as a
                                               ! function of the number of open
                                               ! shells
  integer(is)              :: maxdet           ! Maximum number of determinants
                                               ! across all numbers of open
                                               ! shells
  integer(ib), allocatable :: detvec(:,:)      ! Encoding of the determinants
                                               ! contributing to the CSFs

  !
  ! Pre-computed integrals
  !
  real(dp), allocatable :: fii(:)    ! On-diagonal Fock matrix elements
  real(dp), allocatable :: fock(:,:) ! Fock matrix
  real(dp), allocatable :: Vc(:,:)   ! Coulomb integrals V_iijj
  real(dp), allocatable :: Vx(:,:)   ! Exchange integrals V_ijji
  real(dp)              :: Escf      ! SCF energy

  !
  ! DFT/MRCI logical flag
  !
  logical :: ldftmrci

  !
  ! Hamiltonian build & diagonalisation
  !
  real(dp), parameter :: epshij=1e-6_dp ! Threshold for neglect of Hamiltonian
                                        ! matrix elements in disk-based
                                        ! sigma-vector builds
  real(dp), parameter :: hshift=1e-6_dp ! Hamiltonian shift: useful if the
                                        ! lowest eigenvalue is vanishingly small

!----------------------------------------------------------------------
! State interaction variables
!----------------------------------------------------------------------
  !
  ! Bra and ket spin multiplicities
  !
  integer(is) :: imultB,imultK

  !
  ! Bra and ket numbers of electrons
  !
  integer(is) :: nelB,nelB_alpha,nelB_beta,nelB_spin(2)
  integer(is) :: nelK,nelK_alpha,nelK_beta,nelK_spin(2)
  
  
end module bitglobal

