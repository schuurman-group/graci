!######################################################################
! Top-level routine to compute GVVPT2 QDPT diabatic states and energies
! using pre-computed A-vectors
!######################################################################
#ifdef CBINDING
subroutine gvvpt2_diab(irrep,nroots,nextra,ireg,regfac,n_intR0,&
     ndetR0,nrootsR0,detR0,vecR0,nmoR0,smoR0,ncore,icore,lfrzcore,&
     confscr,vec0scr,Ascr,adtR0) bind(c,name="gvvpt2_diab")
#else
subroutine gvvpt2_diab(irrep,nroots,nextra,ireg,regfac,n_intR0,&
     ndetR0,nrootsR0,detR0,vecR0,nmoR0,smoR0,ncore,icore,lfrzcore,&
     confscr,vec0scr,Ascr,adtR0)
#endif

  use constants
  use bitglobal
  use conftype
  use csf2det
  use utils
  use iomod
  use timing
  
  implicit none

  ! Irrep number
  integer(is), intent(in)  :: irrep

  ! Number of roots requested
  integer(is), intent(in)  :: nroots

  ! Number of buffer states
  integer(is), intent(in)  :: nextra

  ! Regularizer index: 1 <-> ISA
  !                    2 <-> sigma^p
  integer(is), intent(in)  :: ireg
  
  ! Regularization factor
  real(dp), intent(in)     :: regfac

  ! Eigenvectors to follow, expressed in a Slater determinant basis
  integer(is), intent(in)  :: n_intR0,ndetR0,nrootsR0
  integer(ib), intent(in)  :: detR0(n_intR0,2,ndetR0)
  real(dp), intent(in)     :: vecR0(ndetR0,nrootsR0)

   ! MO overlaps
  integer(is), intent(in)  :: nmoR0
  real(dp), intent(in)     :: smoR0(nmoR0,nmo)

  ! Frozen core orbital info
  integer(is), intent(in)  :: ncore
  integer(is), intent(in)  :: icore(ncore)
  logical(is), intent(in)  :: lfrzcore
  
  ! Array of MRCI configuration scratch file numbers
  integer(is), intent(in)  :: confscr(0:nirrep-1)

  ! Array of reference space eigenvector scratch file numbers
  integer(is), intent(in)  :: vec0scr(0:nirrep-1)

  ! A-vector scratch file number
  integer(is), intent(in)  :: Ascr

  ! ADT matrix of the previous geometry
  real(dp), intent(in)     :: adtR0(nrootsR0,nrootsR0)
  
  ! MRCI configuration derived types
  type(mrcfg)              :: cfg,cfg_ref

  ! Zeroth-order eigenpairs
  integer(is)              :: refdim
  real(dp), allocatable    :: E0(:)
  real(dp), allocatable    :: vec0(:,:)
  real(dp), allocatable    :: vec0_det(:,:)
  integer(ib), allocatable :: det_ref(:,:,:)
  
  ! A-vectors
  real(dp), allocatable    :: Avec(:,:)

  ! Wave function overlaps
  integer(is)              :: npairs
  integer(is), allocatable :: ipairs(:,:)
  real(dp), allocatable    :: Sij(:),Smat(:,:)
  real(dp), allocatable    :: smoT(:,:)
  real(dp)                 :: normthrsh
  logical                  :: lprint

  ! Prototype diabatic state expansion coefficients
  real(dp), allocatable    :: pre_coe(:,:)
  real(dp), allocatable    :: proto_coe(:,:)
  
  ! I/O variables
  integer(is)              :: iscratch
  
  ! Timing variables
  real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                              twall_end

  ! Everything else
  integer(is)              :: i,j,n
  integer(is)              :: nvec,ndet_ref
  integer(is)              :: dim1,dim2
  
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  ! Section header
  if (verbose) then
     write(6,'(/,52a)') ('-',i=1,52)
     write(6,'(x,a,/,x,a)') &
          'GVVPT2 diabatisation calculation for the ' &
          //trim(irreplbl(irrep,ipg)),'subspace'
     write(6,'(52a)') ('-',i=1,52)
  endif
  
!----------------------------------------------------------------------
! For now, we will only support nroots = nrootsR0
! i.e., same numbers of reference and displaced geometry states
!----------------------------------------------------------------------
  if (nroots /= nrootsR0) then
     errmsg='Error in gvvpt2_diab: nroots != nrootsR0'
     call error_control
  endif

!----------------------------------------------------------------------
! Output the regularizer being used
!----------------------------------------------------------------------
  select case(ireg)
  case(1)
     if (verbose) write(6,'(/,x,a)') 'Regularizer: ISA'
  case(2)
     if (verbose) write(6,'(/,x,a)') 'Regularizer: sigma^p'
  case default
     errmsg='Error in gvvpt2_diab: unrecognized regularizer index'
     call error_control
  end select
  
!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr(irrep))

  if (verbose) &
       write(6,'(/,x,a,x,i0)') 'CSF basis dimension:',cfg%csfdim

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  ! No. reference space states used to construct the prototype
  ! diabatic states (PSDs)
  ! nvec = no. roots + no. buffer states
  nvec=nroots+nextra

  ! Number of reference space CSFs
  refdim=cfg%csfs0h(cfg%n0h+1)-1

  allocate(Avec(cfg%csfdim,nvec))
  Avec=0.0d0

  allocate(E0(nvec))
  E0=0.0d0

  allocate(vec0(refdim,nvec))
  vec0=0.0d0

  npairs=nvec*nrootsR0
  allocate(Sij(npairs))
  allocate(Smat(nvec,nrootsR0))
  allocate(ipairs(npairs,2))
  Sij=0.0d0
  Smat=0.0d0
  ipairs=0

  allocate(smoT(nmo,nmoR0))
  smoT=0.0d0

  allocate(pre_coe(nvec,nrootsR0))
  pre_coe=0.0d0

  allocate(proto_coe(nvec,nrootsR0))
  proto_coe=0.0d0
  
!----------------------------------------------------------------------
! Read in the zeroth-order eigenpairs
!----------------------------------------------------------------------
  call read_all_eigenpairs(vec0scr(irrep),vec0,E0,refdim,nvec)
  
!----------------------------------------------------------------------
! Read in the A-vectors
!----------------------------------------------------------------------
  ! Open the A-vector file
  iscratch=scrunit(Ascr)
  open(iscratch,file=scrname(Ascr),form='unformatted',status='old')

  ! Check on the dimensions
  read(iscratch) dim2
  read(iscratch) dim1
  if (dim1 /= cfg%csfdim .or. dim2 /= nvec) then
     errmsg='Error in gvvpt2_diab: incorrect dimensions found when ' &
          //'parsing the A-vector file'
     call error_control
  endif

  ! Read in the A-vectors
  read(iscratch) Avec
  
  ! Close the A-vector file
  close(iscratch)

!----------------------------------------------------------------------
! Get the determinant representation of the zeroth-order eigenfunctions
!----------------------------------------------------------------------
  ! Set up a dummy MRCI configuration derived data type to hold the
  ! ref conf info
  allocate(cfg_ref%conf0h(cfg%n_int_I,2,cfg%n0h))
  allocate(cfg_ref%sop0h(cfg%n_int_I,2,cfg%n0h))
  cfg_ref%n0h=cfg%n0h
  cfg_ref%n_int_I=cfg%n_int_I
  cfg_ref%n1I=0
  cfg_ref%n2I=0
  cfg_ref%n1E=0
  cfg_ref%n2E=0
  cfg_ref%n1I1E=0
  cfg_ref%conf0h=cfg%conf0h
  cfg_ref%sop0h=cfg%sop0h
  
  ! Determine the dimension of the reference space determinant basis
  call get_detdim(cfg_ref,ndet_ref)

  ! Allocate arrays
  allocate(det_ref(n_int,2,ndet_ref))
  allocate(vec0_det(ndet_ref,nvec))
  det_ref=0_ib
  vec0_det=0.0d0
  
  ! Compute the determinant representation of the reference space
  ! wave functions
  call det_trans(cfg_ref,cfg%m2c,nvec,refdim,ndet_ref,vec0,vec0_det,&
       det_ref)
  
!----------------------------------------------------------------------
! Compute the overlaps of the reference space wave functions with
! the wave functions of the previous geometry
!----------------------------------------------------------------------
! In the following,
!
! Smat(i,j) = <psi_i^(0)(R_n)|psi_j(R_n-1)>
!
! i.e., Smat(:,j) is the projection of the previous geometry
!       wave functions in terms onto the space spanned by the current
!       geometry reference space wave functions
!----------------------------------------------------------------------
  ! Truncation threshold
  normthrsh=0.999d0

  ! Fill in the array of bra-ket overlaps required
  n=0
  do i=1,nvec
     do j=1,nrootsR0
        n=n+1
        ipairs(n,1)=i
        ipairs(n,2)=j
     enddo
  enddo

  ! MO overlaps
  smoT=transpose(smoR0)
  
  ! Compute the overlaps
  lprint=.true.
  call overlap(nmo,nmoR0,&
       n_int,n_intR0,&
       ndet_ref,ndetR0,&
       nvec,nrootsR0,&
       det_ref,detR0,&
       vec0_det,vecR0,&
       smoT,&
       normthrsh,&
       ncore,icore,lfrzcore,&
       npairs,Sij,ipairs,lprint)
  
  ! Put the overlaps into a more useful form
  do n=1,npairs
     i=ipairs(n,1)
     j=ipairs(n,2)
     Smat(i,j)=Sij(n)
  enddo
  
!----------------------------------------------------------------------
! Transform the overlaps using the previous geometry ADT matrix
! to yield the precursor states
!----------------------------------------------------------------------
  pre_coe=matmul(Smat,adtR0)
  
!----------------------------------------------------------------------
! Output the squared norms of projections of the previous geometry
! diabatic wave functions onto the space spanned by the current
! geometry reference space wave functions
!----------------------------------------------------------------------
  write(6,'(/,x,21a)') ('-', i=1,21)
  write(6,'(4x,a)') 'Precursor state'
  write(6,'(5x,a)') 'squared norms'
  write(6,'(x,21a)') ('-', i=1,21)
  write(6,'(2x,a)') 'I  | <phi#_I|phi#_I>'
  write(6,'(x,21a)') ('-', i=1,21)

  do i=1,nrootsR0
     write(6,'(x,i3,x,a,x,F10.7)') &
          i,'|',dot_product(pre_coe(:,i),pre_coe(:,i))
  enddo

  write(6,'(x,21a)') ('-', i=1,21)

!----------------------------------------------------------------------
! Lowdin's symmetric orthonormalisation of the precursor states to
! yield the prototype diabatic states
!----------------------------------------------------------------------
  proto_coe=pre_coe
  call symm_ortho(nvec,nrootsR0,proto_coe)

!----------------------------------------------------------------------
! Calculation of the complement space states via:
!
! (1) projection of the reference space states onto the space
!     orthogonal to the prototype diabatic states, followed by;
!
! (2) the dropping of the null space
!----------------------------------------------------------------------

  
!----------------------------------------------------------------------
! Block diagonalisation of the Hamiltonian matrix between the prototype
! diabatic state and complement space blocks
!
! The block diagonalisation transformation T will be such that
! ||T-1|| = min.
!
! The thus transformed prototype diabatic states will form our model
! space
!----------------------------------------------------------------------
  
  
  STOP
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  if (verbose) &
       call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'gvvpt2_diab')

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  STOP
  
  return
  
end subroutine gvvpt2_diab
