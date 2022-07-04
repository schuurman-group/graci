!######################################################################
! Top-level routine to compute quasi-diabatic states using 2nd-order
! generalised van Vleck perturbation theory
!######################################################################
#ifdef CBINDING
subroutine gvvpt2_diab(irrep,nroots,nextra,shift,n_intR0,ndetR0,&
     nrootsR0,detR0,vecR0,nmoR0,smoR0,adtR0,ncore,icore,lfrzcore,&
     confscr,vecscr,vec0scr,Qscr,dspscr,adt) bind(c,name="gvvpt2_diab")
#else
subroutine gvvpt2_diab(irrep,nroots,nextra,shift,n_intR0,ndetR0,&
     nrootsR0,detR0,vecR0,nmoR0,smoR0,adtR0,ncore,icore,lfrzcore,&
     confscr,vecscr,vec0scr,Qscr,dspscr,adt)
#endif

  use constants
  use bitglobal
  use conftype
  use hii
  use gvvpt2_hamiltonian
  use csf2det
  use mrciutils
  use utils
  use iomod
  use timing

  implicit none
  
  ! Irrep number
  integer(is), intent(in)  :: irrep

  ! Number of roots requested
  integer(is), intent(in)  :: nroots

  ! Number of extra roots to include in the GVVPT2 calculation
  integer(is), intent(in)  :: nextra

  ! ISA shift
  real(dp), intent(in)     :: shift

  ! Eigenvectors to follow, expressed in a Slater determinant basis
  integer(is), intent(in)  :: n_intR0,ndetR0,nrootsR0
  integer(ib), intent(in)  :: detR0(n_intR0,2,ndetR0)
  real(dp), intent(in)     :: vecR0(ndetR0,nrootsR0)

   ! MO overlaps
  integer(is), intent(in)  :: nmoR0
  real(dp), intent(in)     :: smoR0(nmoR0,nmo)

  ! ADT matrix at the previous geometry
  real(dp), intent(in)     :: adtR0(nrootsR0,nrootsR0)
  
  ! Frozen core orbital info
  integer(is), intent(in)  :: ncore
  integer(is), intent(in)  :: icore(ncore)
  logical(is), intent(in)  :: lfrzcore
  
  ! Array of MRCI configuration scratch file numbers
  integer(is), intent(in)  :: confscr(0:nirrep-1)

  ! Array of reference space eigenvector scratch file numbers
  integer(is), intent(in)  :: vec0scr(0:nirrep-1)

  ! Eigenpair scratch file number
  integer(is), intent(out) :: vecscr

  ! Q-space info scratch file number
  integer(is), intent(out) :: Qscr

  ! Damped strong perturber scratch file number
  integer(is), intent(out) :: dspscr

  ! ADT matrix
  real(dp), intent(out)    :: adt(nroots,nroots)
  
  ! MRCI configuration derived type
  type(mrcfg)              :: cfg,cfg_ref

  ! On-diagonal Hamiltonian matrix elements
  real(dp), allocatable    :: hdiag(:)

  ! On-diagonal Hamiltonian matrix elements
  ! averaged over spin couplings
  real(dp), allocatable    :: averageii(:)
  
  ! 1st-order corrected wavefunctions
  real(dp), allocatable    :: Avec(:,:)

  ! ENPT2 energy corrections (do we still need this?)
  real(dp), allocatable    :: E2(:)

  ! Zeroth-order eigenpairs
  integer(is)              :: refdim
  real(dp), allocatable    :: E_ref(:)
  real(dp), allocatable    :: vec_ref(:,:)
  integer(ib), allocatable :: det_ref(:,:,:)
  real(dp), allocatable    :: vec_ref_det(:,:)
  
  ! QDPT2 energies and mixing coefficients
  real(dp), allocatable    :: EQD(:),mix(:,:),work(:,:)

  ! Determinant representation of the 1st-order corrected
  ! wave functions
  integer(is)              :: ndet
  integer(ib), allocatable :: det(:,:,:)
  real(dp), allocatable    :: Avec_det(:,:)

  ! Wave function overlaps
  integer(is)              :: npairs
  integer(is), allocatable :: ipairs(:,:)
  real(dp), allocatable    :: Sij(:),Smat(:,:)
  real(dp)                 :: normthrsh

  ! Representation of the prototype diabatic states in the basis
  ! of the ref space eigenfunctions
  real(dp), allocatable    :: Cmat(:,:)

  ! Representation of the prototype diabatic states in the basis
  ! of the ref space eigenfunctions
  real(dp), allocatable    :: vec_proto(:,:)
  
  ! I/O variables
  integer(is)              :: iscratch
  character(len=60)        :: vecfile,Qfile
  character(len=2)         :: amult,airrep
  
  ! Everything else
  integer(is)              :: i,j,n,imax(1)
  integer(is)              :: n_int_I,nconf_ref,ndet_ref
  integer(is)              :: nvec
  integer(is), allocatable :: indx(:),iphase(:)
  real(dp), allocatable    :: Esort(:)
  real(dp), allocatable    :: Qnorm(:),Qener(:)
  real(dp)                 :: norm
  logical                  :: lprint
  
  ! Timing variables
  real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                              twall_end

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
! i.e., same numbers of roots and roots to follow
!----------------------------------------------------------------------
  if (nroots /= nrootsR0) then
     errmsg='Error in gvvpt2_follow: nroots != nrootsR0'
     call error_control
  endif
  
!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr(irrep))

  if (verbose) &
       write(6,'(/,x,a,x,i0)') 'CSF basis dimension:',cfg%csfdim

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  ! Total number of roots for which the ENPT2 corrections will be
  ! calculated
  nvec=nroots+nextra

  ! Number of reference space CSFs
  refdim=cfg%csfs0h(cfg%n0h+1)-1
  
  allocate(hdiag(cfg%csfdim))
  hdiag=0.0d0

  allocate(averageii(cfg%confdim))
  averageii=0.0d0

  allocate(Avec(cfg%csfdim,nroots))
  Avec=0.0d0

  allocate(E_ref(nvec))
  E_ref=0.0d0

  allocate(vec_ref(refdim,nvec))
  vec_ref=0.0d0
  
  allocate(E2(nroots))
  E2=0.0d0

  allocate(Qnorm(nroots))
  Qnorm=0.0d0

  allocate(Qener(nroots))
  Qener=0.0d0

  allocate(EQD(nroots))
  EQD=0.0d0
  
  allocate(mix(nroots,nroots))
  mix=0.0d0
  
  allocate(work(cfg%csfdim,nroots))
  work=0.0d0

  allocate(indx(nroots))
  indx=0.0d0

  allocate(iphase(nroots))
  iphase=0
  
  allocate(Smat(nroots,nroots))
  Smat=0.0d0

  allocate(Esort(nroots))
  Esort=0.0d0
  
!----------------------------------------------------------------------
! Read in the zeroth-order eigenpairs
!----------------------------------------------------------------------
  call read_all_eigenpairs(vec0scr(irrep),vec_ref,E_ref,refdim,nvec)

!----------------------------------------------------------------------
! Get the determinant representation of the zeroth-order eigenstates
!----------------------------------------------------------------------
  ! Set up a dummy MRCI configuration derived data type to hold
  ! the ref conf info
  n_int_I=cfg%n_int_I
  nconf_ref=cfg%n0h
  allocate(cfg_ref%conf0h(n_int_I,2,nconf_ref))
  allocate(cfg_ref%sop0h(n_int_I,2,nconf_ref))
  cfg_ref%n0h=nconf_ref
  cfg_ref%n_int_I=n_int_I
  cfg_ref%n1I=0
  cfg_ref%n2I=0
  cfg_ref%n1E=0
  cfg_ref%n2E=0
  cfg_ref%n1I1E=0  
  cfg_ref%conf0h=cfg%conf0h
  cfg_ref%sop0h=cfg%sop0h

  ! Determine the dimension of the determinant basis
  call get_detdim(cfg_ref,ndet_ref)

  ! Allocate arrays
  allocate(det_ref(n_int,2,ndet_ref))
  allocate(vec_ref_det(ndet_ref,nvec))
  det_ref=0_ib
  vec_ref_det=0.0d0

  ! Compute the determinant representation of the wave functions
  call det_trans(cfg_ref,cfg%m2c,nvec,refdim,ndet_ref,vec_ref,&
       vec_ref_det,det_ref)
  
!----------------------------------------------------------------------
! Compute the overlaps of the quasi-diabatic states from the
! previous geometry with the zeroth-order eigenstates at the current
! geometry
!----------------------------------------------------------------------
  ! Allocate arrays
  npairs=nrootsR0*nvec
  allocate(Sij(npairs))
  allocate(ipairs(npairs,2))
  allocate(Cmat(nvec,nrootsR0))
  Sij=0.0d0
  ipairs=0
  Cmat=0.0d0

  ! Tight truncation threshold
  normthrsh=0.999d0

  ! Fill in the array of required bra-ket overlaps
  n=0
  do i=1,nrootsR0
     do j=1,nvec
        n=n+1
        ipairs(n,1)=i
        ipairs(n,2)=j
     enddo
  enddo

  ! Compute the overlaps
  lprint=.false.
  call overlap(nmoR0,nmo,n_intR0,n_int,ndetR0,ndet_ref,nrootsR0,nvec,&
       detR0,det_ref,vecR0,vec_ref_det,smoR0,normthrsh,ncore,icore,&
       lfrzcore,npairs,Sij,ipairs,lprint)

  ! Put the overlaps into a more computationally useful form
  n=0
  do i=1,nrootsR0
     do j=1,nvec
        n=n+1
        Cmat(j,i)=Sij(n)
     enddo
  enddo

  ! Deallocate arrays
  deallocate(Sij)
  deallocate(ipairs)
  
!----------------------------------------------------------------------
! Construct the prototype diabatic states
!----------------------------------------------------------------------
  ! Transform using the previous geometry ADT matrix to yield
  ! the prototype diabatic states in the basis of the current geometry
  ! ref space eigenstates
  Cmat=matmul(Cmat,adtR0)

  ! Lowdin's symmetric orthonormalisation of the prototype diabatic
  ! states
  call symm_ortho(nvec,nrootsR0,Cmat)

  ! Prototype diabatic states in the basis of the ref space CSFs
  allocate(vec_proto(refdim,nrootsR0))
  vec_proto=matmul(vec_ref,Cmat)

!----------------------------------------------------------------------
! Compute the on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
  call hmat_diagonal(hdiag,cfg%csfdim,averageii,cfg%confdim,cfg)
  
!----------------------------------------------------------------------
! Construct the GVVPT2 effective Hamiltonian and determine its
! eigenpairs
! Also returns the 1st-order perturbed model functions in the Avec
! array
!----------------------------------------------------------------------
  call gvvpt2_heff(irrep,cfg,hdiag,averageii,cfg%csfdim,cfg%confdim,&
       vec0scr(irrep),Cmat,Avec,E2,nroots,nvec,shift,dspscr,EQD,mix)

!----------------------------------------------------------------------
! ADT matrix
!----------------------------------------------------------------------
  adt=transpose(mix)
  
!----------------------------------------------------------------------
! Add in the prototype diabatic wave functions to obtain the 1st-order
! diabatic wave functions
!----------------------------------------------------------------------
  do i=1,nroots
     Avec(1:refdim,i)=vec_proto(:,i)
  enddo

!----------------------------------------------------------------------
! Get the determinant representation of the 1st-order diabatic
! wave functions
!----------------------------------------------------------------------
  ! Determine the dimension of the determinant basis
  call get_detdim(cfg,ndet)

  ! Allocate arrays
  allocate(det(n_int,2,ndet))
  allocate(Avec_det(ndet,nvec))
  det=0_ib
  Avec_det=0.0d0

  ! Compute the determinant representation of the wave functions
  call det_trans(cfg,cfg%m2c,nroots,cfg%csfdim,ndet,Avec,Avec_det,det)

  ! Normalisation (we will let the overlap routine do the
  ! orthogonalisation)
  do i=1,nroots
     norm=sqrt(dot_product(Avec_det(:,i),Avec_det(:,i)))
     Avec_det(:,i)=Avec_det(:,i)/norm
  enddo
  
!----------------------------------------------------------------------
! Get the overlaps of the diabatic states with those from the previous
! geometry
!----------------------------------------------------------------------
  ! Allocate arrays
  npairs=nrootsR0*nroots
  allocate(Sij(npairs))
  allocate(ipairs(npairs,2))
  Sij=0.0d0
  ipairs=0

  ! Loose truncation threshold
  normthrsh=0.99d0

  ! Fill in the array of required bra-ket overlaps
  n=0
  do i=1,nrootsR0
     do j=1,nroots
        n=n+1
        ipairs(n,1)=i
        ipairs(n,2)=j
     enddo
  enddo

  ! Compute the overlaps between the previous geometry adiabatic states
  ! and the 1st-order diabatic states
  lprint=.false.
  call overlap(nmoR0,nmo,n_intR0,n_int,ndetR0,ndet,nrootsR0,nroots,&
       detR0,det,vecR0,Avec_det,smoR0,normthrsh,ncore,icore,&
       lfrzcore,npairs,Sij,ipairs,lprint)

  ! Put the overlaps into a more computationally useful form
  do n=1,npairs
     i=ipairs(n,1)
     j=ipairs(n,2)
     Smat(i,j)=Sij(n)
  enddo

  ! Transform the matrix of overlaps using the previous geometry ADT
  ! matrix to get the overlaps between diabatic states
  Smat=matmul(transpose(adtR0),Smat)

!----------------------------------------------------------------------
! Re-phasing
!----------------------------------------------------------------------
  ! Get the phase factors
  do i=1,nroots
     imax=maxloc(abs(Smat(i,:)))
     j=imax(1)
     if (Smat(i,j) < 0.0d0) then
        iphase(j)=-1
     else
        iphase(j)=1
     endif
     Smat(:,j)=0.0d0
  enddo

  ! Re-phasing of the ADT matrix
  do i=1,nroots
     adt(i,:)=iphase(i)*adt(i,:)
  enddo
  
!----------------------------------------------------------------------
! Compute the 1st-order adiabtic wave functions
!----------------------------------------------------------------------
  work=Avec
  Avec=matmul(work,mix)

!----------------------------------------------------------------------
! Sort the 2nd-order energies
!----------------------------------------------------------------------
  ! Sorting
  call dsortindxa1('A',nroots,EQD,indx)

  ! Lowest-lying energies
  do i=1,nroots
     Esort(i)=EQD(indx(i))
  enddo

!----------------------------------------------------------------------
! Write the Q-space information to disk
!----------------------------------------------------------------------
  ! Norms of the 1st-order wave functions projected onto the Q-space
  do i=1,nroots
     Qnorm(i)=sqrt(dot_product(Avec(refdim+1:cfg%csfdim,indx(i)),&
          Avec(refdim+1:cfg%csfdim,indx(i))))
  enddo

  ! ENPT2 energy corrections
  ! (do we still need these?)
  do i=1,nroots
     Qener(i)=E2(indx(i))
  enddo

  ! Register the scratch file
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name('Qinfo'//'.mult'//trim(amult)//&
       '.sym'//trim(airrep),Qfile)
  call register_scratch_file(Qscr,Qfile)

  ! Open the scratch file
  iscratch=scrunit(Qscr)
  open(iscratch,file=scrname(Qscr),form='unformatted',&
       status='unknown')

  ! Number of roots
  write(iscratch) nroots
  
  ! Norms of the wave function corrections
  write(iscratch) Qnorm

  ! Energy corrections
  write(iscratch) Qener
  
  ! Close the scratch file
  close(iscratch)

!----------------------------------------------------------------------
! Lowdin's symmetric orthonormalisation of the 1st-order corrected
! adiabatic wave functions
!----------------------------------------------------------------------
  call symm_ortho(cfg%csfdim,nroots,Avec)
  
!----------------------------------------------------------------------
! Output the Q-space norms
!----------------------------------------------------------------------
  if (verbose) then
     
     ! Table header
     write(6,'(/,x,19a)') ('-',i=1,19)
     write(6,'(4x,a)') 'Q-space info'
     write(6,'(x,19a)') ('-',i=1,19)
     write(6,'(2x,a)') 'State   ||psi_Q||'
     write(6,'(x,19a)') ('-',i=1,19)

     ! A-vector norms
     do i=1,nroots
        write(6,'(2x,i4,3x,F9.6)') &
             i,Qnorm(i)
     enddo

     ! Table footer
     write(6,'(x,19a)') ('-',i=1,19)
     
  endif

!----------------------------------------------------------------------
! Write the 2nd-order energies and 1st-order wave functions to disk
!----------------------------------------------------------------------
  ! Register the scratch file
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name('gvvpt2vec'//'.mult'//trim(amult)//&
       '.sym'//trim(airrep),vecfile)
  call register_scratch_file(vecscr,vecfile)
  
  ! Open the scratch file
  iscratch=scrunit(vecscr)
  open(iscratch,file=scrname(vecscr),form='unformatted',&
       status='unknown')
  
  ! No. CSFs
  write(iscratch) cfg%csfdim
  
  ! No. roots
  write(iscratch) nroots
  
  ! 2nd-order energies
  write(iscratch) Esort
    
  ! 1st-order wave functions
  do i=1,nroots
     write(iscratch) Avec(:,indx(i))
  enddo
  
  ! Close the scratch file
  close(iscratch)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(hdiag)
  deallocate(averageii)
  deallocate(Avec)
  deallocate(Avec_det)
  deallocate(E_ref)
  deallocate(vec_ref)
  deallocate(E2)
  deallocate(Qnorm)
  deallocate(Qener)
  deallocate(EQD)
  deallocate(mix)
  deallocate(work)
  deallocate(indx)
  deallocate(Smat)
  deallocate(Esort)
  deallocate(det_ref)
  deallocate(vec_ref_det)
  deallocate(Cmat)
  deallocate(vec_proto)
  deallocate(det)
  deallocate(Sij)
  deallocate(ipairs)
  deallocate(iphase)
  
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

  return
  
end subroutine gvvpt2_diab
  
