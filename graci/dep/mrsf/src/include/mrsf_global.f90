!**********************************************************************
! mrsf_global: module-level state of the MRSF-TDDFT library
!
! Orbital classes (from the ROKS/ROHF triplet occupation vector):
!   C  = doubly occupied (nC), O1,O2 = singly occupied, V = virtual (nV)
!   hole set     H = C + [O1,O2]  (alpha occupied), nocca = nC+2
!   particle set P = [O1,O2] + V  (beta virtual),   nvirb = nV+2
! Amplitudes X(a,i), a in P (local index), i in H (local index),
! packed index ia = a + nvirb*(i-1).  Local positions of the SOMOs:
!   particle: pO1=1, pO2=2;   hole: hO1=nC+1, hO2=nC+2
!**********************************************************************
module mrsf_global

  use mrsf_constants

  implicit none
  save

  ! dimensions
  integer(is) :: nmo = 0, nel = 0, imult_ref = 3
  integer(is) :: nocca = 0, nC = 0, nV = 0, nvirb = 0, xdim = 0
  integer(is) :: naux = 0, nplane = 0, nQ = 0, nblk = 0, npblk = 0
  integer(is) :: nirrep = 1, ipg = 1
  integer(is) :: iO1 = 0, iO2 = 0

  ! flags
  logical     :: verbose = .true.
  logical     :: init_done = .false.
  logical     :: ints_loaded = .false.
  logical     :: store_sp = .false.
  logical     :: vv_full = .false.

  character(len=255) :: label = ''

  ! parameters
  real(dp)    :: chf = 1.0_dp
  real(dp)    :: spc(3) = 1.0_dp     ! kappa_coco, kappa_ovov, kappa_coov
  real(dp)    :: escf = 0.0_dp
  real(dp)    :: mem_budget = 1.0e9_dp

  ! orbital data
  real(dp), allocatable    :: occ(:), moen(:)
  integer(is), allocatable :: mosym(:)
  integer(is), allocatable :: Hmap(:), Pmap(:)   ! local -> MO index
  integer(is), allocatable :: Hinv(:), Pinv(:)   ! MO -> local index (0: none)
  integer(is), allocatable :: slot_irrep(:)      ! irrep of each packed slot
  real(dp), allocatable    :: FaHH(:,:), FbPP(:,:)

  ! three-index integrals
  real(dp), allocatable    :: Bvv(:,:,:)      ! (nvirb,nvirb,nplane|naux)
  real(sp), allocatable    :: Bvv_sp(:,:,:)
  real(dp), allocatable    :: Dall(:,:)       ! (naux,nmo) diagonals B^Q_pp
  real(dp), allocatable    :: Boo(:,:,:,:)    ! (nQ,nocca,nocca,nblk)
  real(sp), allocatable    :: Boo_sp(:,:,:,:)
  real(dp), allocatable    :: Bco(:,:,:)      ! (naux,nC,2)
  real(dp), allocatable    :: Bvo(:,:,:)      ! (naux,nV,2)

  ! pairing-strength integral blocks
  real(dp), allocatable    :: Gp(:,:,:,:)     ! (nC,nC,2,2)  (i O_x|j O_y)
  real(dp), allocatable    :: Mp(:,:,:,:)     ! (nV,nV,2,2)  (a O_x|b O_y)
  real(dp), allocatable    :: Np(:,:,:,:)     ! (nC,nV,2,2)  (i O_x|O_y a)
  real(dp), allocatable    :: Hp(:,:)         ! (nC,nV) N12 - N21
  real(dp)                 :: K12 = 0.0_dp    ! (O1 O2|O1 O2)

  ! multiplicity-independent part of the diagonal
  real(dp), allocatable    :: diag0(:,:)      ! (nvirb,nocca)

  ! work arrays
  real(dp), allocatable    :: Twork(:,:,:,:)  ! (nvirb,nQ,nocca,nvmax)
  real(dp), allocatable    :: plane_scr(:,:)  ! (nvirb,nvirb)
  real(dp), allocatable    :: Boo_scr(:,:,:)  ! (nQ,nocca,nocca)
  integer(is)              :: nvmax = 0

  ! key of the currently loaded integrals (for reuse across sections)
  character(len=255)       :: loaded_file = '', loaded_prec = '', loaded_vv = ''
  integer(is)              :: loaded_nmo = 0
  real(dp), allocatable    :: loaded_occ(:)

  ! timings and counters
  real(dp)    :: time_load = 0.0_dp, time_sigma = 0.0_dp, time_exch = 0.0_dp
  integer(is) :: nsigma_calls = 0, nsigma_vecs = 0

end module mrsf_global
