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
!
! Extended MRSF-TDDFT (EMRSF): the vector gains nC further columns, the
! core-to-virtual amplitudes y(a,i), a in P (SOMO rows masked), i in C,
! of the closed-shell configuration G = |C Cbar O1 O1bar|; a vector is
! the (nvirb x ncol) matrix [x | y], ncol = nocca + nC, xdim_tot =
! nvirb*ncol (xdim = nvirb*nocca remains the MRSF part).
!**********************************************************************
module mrsf_global

  use mrsf_constants

  implicit none
  save

  ! dimensions
  integer(is) :: nmo = 0, nel = 0, imult_ref = 3
  integer(is) :: nocca = 0, nC = 0, nV = 0, nvirb = 0, xdim = 0
  integer(is) :: ncol = 0, xdim_tot = 0, ncv = 0
  integer(is) :: naux = 0, nplane = 0, nQ = 0, nblk = 0, npblk = 0
  integer(is) :: nirrep = 1, ipg = 1
  integer(is) :: iO1 = 0, iO2 = 0

  ! flags
  logical     :: verbose = .true.
  logical     :: init_done = .false.
  logical     :: ints_loaded = .false.
  logical     :: store_sp = .false.
  logical     :: vv_full = .false.
  logical     :: extended = .false.     ! EMRSF: CV columns present
  logical     :: ext_ready = .false.    ! ext_initialise done
  ! frozen core: nfc doubly occupied MOs (the lowest by orbital energy)
  ! excluded from the hole set of the response space; their slots are
  ! masked. nfc_req is set before mrsf_initialise (mrsf_set_frozen_core)
  integer(is) :: nfc_req = 0, nfc = 0, nocca_act = 0, nC_act = 0
  ! exchange term: explicit (ij|ab) tensor (mrsf_etensor) or DF plane sweep;
  ! exchange_mode 0 = auto (tensor when it fits mem_budget), 1 = tensor, 2 = df
  integer(is) :: exchange_mode = 0
  logical     :: use_etensor = .false.

  character(len=255) :: label = ''

  ! parameters
  real(dp)    :: chf = 1.0_dp
  real(dp)    :: spc(3) = 1.0_dp     ! kappa_coco, kappa_ovov, kappa_coov
  real(dp)    :: escf = 0.0_dp
  real(dp)    :: mem_budget = 1.0e9_dp
  real(dp)    :: ccp = 1.0_dp        ! EMRSF coupling scale c_cp
  logical     :: use_kernel = .false.! EMRSF: XC kernel in the singlet CV block

  ! orbital data
  real(dp), allocatable    :: occ(:), moen(:)
  integer(is), allocatable :: mosym(:)
  integer(is), allocatable :: Hmap(:), Pmap(:)   ! local -> MO index
  integer(is), allocatable :: Hinv(:), Pinv(:)   ! MO -> local index (0: none)
  integer(is), allocatable :: slot_irrep(:)      ! irrep of each packed slot (xdim_tot)
  logical, allocatable     :: frozen_hole(:)     ! (nocca) hole excluded from the response
  integer(is), allocatable :: hact(:), cact(:)   ! active hole / core local indices
  real(dp), allocatable    :: FaHH(:,:), FbPP(:,:)

  ! three-index integrals
  real(dp), allocatable    :: Bvv(:,:,:)      ! (nvirb,nvirb,nplane|naux)
  real(sp), allocatable    :: Bvv_sp(:,:,:)
  real(dp), allocatable    :: Dall(:,:)       ! (naux,nmo) diagonals B^Q_pp
  real(dp), allocatable    :: Boo(:,:,:,:)    ! (nQ,nocca,nocca,nblk)
  real(sp), allocatable    :: Boo_sp(:,:,:,:)
  real(dp), allocatable    :: Bco(:,:,:)      ! (naux,nC,2)
  real(dp), allocatable    :: Bvo(:,:,:)      ! (naux,nV,2)
  real(dp), allocatable    :: Bcv(:,:,:)      ! (nvirb,nC,naux) B^Q_{a,i}, a in P, i in C (EMRSF)

  ! pairing-strength integral blocks
  real(dp), allocatable    :: Gp(:,:,:,:)     ! (nC,nC,2,2)  (i O_x|j O_y)
  real(dp), allocatable    :: Mp(:,:,:,:)     ! (nV,nV,2,2)  (a O_x|b O_y)
  real(dp), allocatable    :: Np(:,:,:,:)     ! (nC,nV,2,2)  (i O_x|O_y a)
  real(dp), allocatable    :: Hp(:,:)         ! (nC,nV) N12 - N21
  real(dp)                 :: K12 = 0.0_dp    ! (O1 O2|O1 O2)

  ! multiplicity-independent part of the diagonal
  real(dp), allocatable    :: diag0(:,:)      ! (nvirb,nocca)

  ! EMRSF data (ext_initialise)
  real(dp), allocatable    :: Fp_emb(:,:)     ! (nvirb,nvirb) F'_VV embedded (SOMO rows/cols 0)
  real(dp), allocatable    :: Fp_cc(:,:)      ! (nC,nC) F'_CC
  real(dp), allocatable    :: Fcv(:,:)        ! (nV,nC) F^DFT_{b j}
  real(dp), allocatable    :: fO2V(:)         ! (nV) F^DFT_{O2 b}
  real(dp), allocatable    :: fCO1(:)         ! (nC) F^DFT_{j O1}
  real(dp), allocatable    :: w1(:)           ! (nV) (O2O1|O2 b)
  real(dp), allocatable    :: w2(:)           ! (nC) (jO1|O2O1)
  real(dp), allocatable    :: Bo12(:)         ! (naux) B^Q_{O1 O2}
  real(dp), allocatable    :: diag_cv(:,:)    ! (nvirb,nC) CV diagonal without the shift
  real(dp)                 :: A_G = 0.0_dp    ! shift of the CV block

  ! work arrays
  real(dp), allocatable    :: Twork(:,:,:,:)  ! (nvirb,nQ,ncol,nvmax)
  real(dp), allocatable    :: plane_scr(:,:)  ! (nvirb,nvirb)
  real(dp), allocatable    :: Boo_scr(:,:,:)  ! (nQ,nocca,nocca)
  integer(is)              :: nvmax = 0

  ! key of the currently loaded integrals (for reuse across sections)
  character(len=255)       :: loaded_file = '', loaded_prec = '', loaded_vv = ''
  integer(is)              :: loaded_nmo = 0
  real(dp), allocatable    :: loaded_occ(:)

  ! timings and counters
  real(dp)    :: time_load = 0.0_dp, time_sigma = 0.0_dp, time_exch = 0.0_dp
  real(dp)    :: time_ext = 0.0_dp, time_xcs = 0.0_dp, time_inner = 0.0_dp
  integer(is) :: nsigma_calls = 0, nsigma_vecs = 0, nxcs_vecs = 0, ninner_vecs = 0

end module mrsf_global
