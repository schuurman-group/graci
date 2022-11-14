!######################################################################
! Top-level routine to compute GVVPT2 QDPT diabatic states and energies
! using pre-computed A-vectors
!######################################################################
#ifdef CBINDING
subroutine gvvpt2_diab(irrep,nroots,nextra,ireg,regfac,n_intR0,&
     ndetR0,nrootsR0,detR0,vecR0,nmoR0,smoR0,ncore,icore,lfrzcore,&
     confscr,vec0scr,Ascr) bind(c,name="gvvpt2_diab")
#else
subroutine gvvpt2_diab(irrep,nroots,nextra,ireg,regfac,n_intR0,&
     ndetR0,nrootsR0,detR0,vecR0,nmoR0,smoR0,ncore,icore,lfrzcore,&
     confscr,vec0scr,Ascr)
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
  real(dp)                 :: normthrsh
  logical                  :: lprint
  
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

  npairs=nrootsR0*nvec
  allocate(Sij(npairs))
  allocate(Smat(nvec,nrootsR0))
  allocate(ipairs(npairs,2))
  Sij=0.0d0
  Smat=0.0d0
  ipairs=0
  
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
  allocate(vec0_det(refdim,nvec))
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

  ! Compute the overlaps
  lprint=.true.
  call overlap(nmoR0,nmo,n_intR0,n_int,ndetR0,ndet_ref,nrootsR0,nvec,&
       detR0,det_ref,vecR0,vec0_det,smoR0,normthrsh,ncore,icore,&
       lfrzcore,npairs,Sij,ipairs,lprint)

  ! Put the overlaps into a more useful form
  do n=1,npairs
     i=ipairs(n,1)
     j=ipairs(n,2)
     Smat(i,j)=Sij(n)
  enddo
  
!----------------------------------------------------------------------
! Transform the overlaps using the previous geometry ADT matrix
!----------------------------------------------------------------------

  print*,''
  do i=1,nrootsR0
     print*,i,dot_product(Smat(:,i),Smat(:,i))
  enddo

  print*,''
  print*,'The R0 ADT matrix needs passing...'
  
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
