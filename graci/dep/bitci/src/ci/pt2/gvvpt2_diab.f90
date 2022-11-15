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
  use protodiab
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

  ! Previous geometry eigensates expressed in a Slater determinant
  ! basis
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
  type(mrcfg)              :: cfg

  ! Zeroth-order eigenpairs
  integer(is)              :: refdim
  real(dp), allocatable    :: E0(:)
  real(dp), allocatable    :: vec0(:,:)
  
  ! A-vectors
  real(dp), allocatable    :: Avec(:,:)

  ! Prototype diabatic states
  real(dp), allocatable    :: vec_pds(:,:)

  ! Complement states
  integer(is)              :: ncomp
  real(dp), allocatable    :: vec_comp(:,:)

  ! Model space states
  real(dp), allocatable    :: vec_mod(:,:)
  
  ! I/O variables
  integer(is)              :: iscratch
  
  ! Timing variables
  real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                              twall_end
  
  ! Everything else
  integer(is)              :: i,j,n
  integer(is)              :: nvec
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
  
!----------------------------------------------------------------------
! Read in the zeroth-order eigenpairs
!----------------------------------------------------------------------
  call read_all_eigenpairs(vec0scr(irrep),vec0,E0,refdim,nvec)

!----------------------------------------------------------------------
! Compute the prototype diabatic states
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(vec_pds(refdim,nrootsR0))
  vec_pds=0.0d0

  ! Compute the prototype diabatic states in the ref CSF basis
  call get_pds_basis(cfg,refdim,nvec,vec0,nmoR0,n_intR0,ndetR0,&
       nrootsR0,detR0,vecR0,smoR0,ncore,icore,lfrzcore,adtR0,vec_pds)

!----------------------------------------------------------------------
! Calculation of the complement space states
!----------------------------------------------------------------------
  ! Allocate arrays
  ncomp=nvec-nrootsR0
  allocate(vec_comp(refdim,ncomp))
  vec_comp=0.0d0

  ! Compute the complement states in the ref CSF basis
  call get_complement_basis(refdim,nrootsR0,nvec,ncomp,vec0,vec_pds,&
       vec_comp)

!----------------------------------------------------------------------
! Calculation of the model space states
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(vec_mod(refdim,nrootsR0))
  vec_mod=0.0d0

  ! Compute the model space states in the ref CSF basis
  call get_model_basis(nvec,nrootsR0,ncomp,refdim,E0,vec_pds,&
       vec_comp,vec_mod)
    
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
