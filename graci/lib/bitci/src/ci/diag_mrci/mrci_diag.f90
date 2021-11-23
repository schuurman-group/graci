!**********************************************************************
! Routines for the diagonalisation if the MRCI Hamiltonian
!**********************************************************************

!######################################################################
! diag_mrci: Controls the diagonalisation of the MRCI Hamiltonian for
!            a given symmetry subspace
!######################################################################
#ifdef CBINDING
subroutine diag_mrci(irrep,nroots,confscr,vecscr,ialg,tol,niter,&
     blocksize,ldeflate,iguess,vec0scr) bind(c,name="diag_mrci")
#else
subroutine diag_mrci(irrep,nroots,confscr,vecscr,ialg,tol,niter,&
     blocksize,ldeflate,iguess,vec0scr)
#endif

  use constants
  use bitglobal
  use conftype
  use hii
  use hij_disk
  use hbuild_double
  use mrci_guess
  use full_diag
  use blockdav
  use gendav
  use mrciutils
  use iomod
  
  implicit none
  
  ! Irrep number
  integer(is), intent(in)  :: irrep

  ! Number of roots requested
  integer(is), intent(in)  :: nroots
  
  ! Array of MRCI configuration scratch file numbers
  integer(is), intent(in)  :: confscr(0:nirrep-1)

  ! Array of reference space eigenvector scratch file numbers
  integer(is), intent(in)  :: vec0scr(0:nirrep-1)

  ! Eigenpair scratch file number
  integer(is), intent(out) :: vecscr

  ! Davidson parameters
  integer(is), intent(in)  :: ialg,blocksize,niter,iguess
  real(dp), intent(in)     :: tol
  logical, intent(in)      :: ldeflate
  integer(is)              :: maxvec
  
  ! On-diagonal Hamiltonian matrix elements
  real(dp), allocatable    :: hdiag(:)

  ! On-diagonal Hamiltonian matrix elements
  ! averaged over spin couplings
  real(dp), allocatable    :: averageii(:)

  ! MRCI configuration derived type
  type(mrcfg)              :: cfg
  
  ! Direct-mode logical flag
  logical                  :: direct

  ! Dimension of the CSF basis past which we will switch to
  ! the direct mode of operations
  !integer(is), parameter   :: disk_lim=750000
  integer(is), parameter   :: disk_lim=7500000
  
  ! Dimension of the CSF basis past which full diagonalisation
  ! will not be used
  integer(is), parameter   :: full_lim=1000
  
  ! Subspace dimension to be used in the generation of the
  ! guess vectors
  integer(is)              :: guessdim
  
  ! Disk-based sigma-vector build variables
  integer(is)              :: hamscr,nrec

  ! Guess vector scratch file number
  integer(is)              :: guessscr

  ! Preconditioner integer flag
  integer(is)              :: ipre
  
  ! Everything else
  integer(is)              :: i,k
  integer(is)              :: isigma(3)
  real(dp)                 :: mem
  character(len=20)        :: vecstem

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  write(6,'(/,52a)') ('-',i=1,52)
  write(6,'(3(x,a))') 'MRCI space diagonalisation in the',&
       trim(irreplbl(irrep,ipg)),'subspace'
  write(6,'(52a)') ('-',i=1,52)

!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr(irrep))

  write(6,'(/,x,a,x,i0)') 'CSF basis dimension:',cfg%csfdim

!----------------------------------------------------------------------
! Memory required to store the configuration bit strings
!----------------------------------------------------------------------
  mem=32*n_int*cfg%confdim/1024.0d0**2
  if (mem < 1000.0d0) then
     write(6,'(/,x,a,x,F7.2,x,a)') &
          'Configuration bit strings require',mem,'MB'
  else
     write(6,'(/,x,a,x,F7.2,x,a)') &
          'Configuration bit strings require',mem/1024.0d0,'GB'
  endif

!----------------------------------------------------------------------
! Determine whether or not this will be a direct-mode calculation
! For now, we will just use a hard-wired limit on the CSF basis
! size
!----------------------------------------------------------------------
  if (cfg%csfdim > disk_lim) then
     direct=.true.
     write(6,'(/,x,a)') 'Using direct sigma-vector builds'
  else
     direct=.false.
     write(6,'(/,x,a)') 'Using disk-based sigma-vector builds'
  endif

!----------------------------------------------------------------------
! Compute the on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(hdiag(cfg%csfdim))
  allocate(averageii(cfg%confdim))
  hdiag=0.0d0
  averageii=0.0d0
  
  ! Compute the diagonal elements
  call hmat_diagonal(hdiag,cfg%csfdim,averageii,cfg%confdim,cfg)

!----------------------------------------------------------------------
! If we are using disk-based sigma-vector builds, then save the
! non-zero off-diagonal elements of the Hamiltonian matrix to disk
!----------------------------------------------------------------------
  if (.not. direct) &
       call save_hij(hamscr,nrec,irrep,averageii,cfg%confdim,cfg)
  
!----------------------------------------------------------------------
! Eigenpair scratch file stem
!----------------------------------------------------------------------
  vecstem='mrcivec'
  
!----------------------------------------------------------------------
! Full diagonalisation
!----------------------------------------------------------------------
  if (cfg%csfdim <= full_lim) &
       call diag_full(hamscr,nrec,cfg%csfdim,hdiag,irrep,nroots,&
       vecscr,vecstem)

!----------------------------------------------------------------------
! Iterative diagonalisation
!----------------------------------------------------------------------
  if (cfg%csfdim > full_lim) then

     ! Maximum subspace dimension
     maxvec=6*blocksize

     ! Dimension of the subspace used to generate the guess vectors
     guessdim=750
     if (guessdim > cfg%csfdim) guessdim=int(cfg%csfdim*0.9d0)

     ! Generate the guess vectors
     select case(iguess)
     case(1) ! Diagonalisation within a small subspace
        call mrci_guess_subspace(guessscr,blocksize,cfg,cfg%csfdim,&
             guessdim,hdiag,averageii,cfg%confdim)
     case(2) ! ENPT2 1st-order corrected wave functions
        call mrci_guess_enpt2(guessscr,blocksize,cfg,cfg%csfdim,&
             hdiag,averageii,cfg%confdim,vec0scr(irrep))
     case default
        errmsg='Error in diag_mrci: unrecognised value of iguess'
        call error_control
     end select

     ! Set the sigma-vector algorithm information
     if (direct) then
        isigma=0
        errmsg='Direct-mode diagonalisation not yet implemented'
        call error_control
     else
        isigma(1)=1
        isigma(2)=hamscr
        isigma(3)=nrec
     endif

     ! Generalised Davidson preconditioner: DPR is generally faster
     ! for small numbers of CSFs
     if (ialg == 1) then
        if (cfg%csfdim <= 50000) then
           ! Diagonal preconditiond residue
           ipre=1
        else
           ! Generalised Davidson preconditioner
           ipre=2
        endif
     endif
     
     ! ENPT2 guess vectors and the generalised Davidson preconditioner
     ! are not yet available due to laziness
     if (iguess == 2 .and. ipre == 2) ipre=1
     
     ! Perform the iterative diagonalisation
     select case(ialg)
     case(1)
        ! Generalised Davidson
        call generalised_davidson(irrep,isigma,cfg,cfg%csfdim,hdiag,&
             guessscr,vecscr,vecstem,nroots,blocksize,maxvec,tol,&
             niter,ipre)
     case(2)
        ! Block Davidson
        call block_davidson(irrep,isigma,cfg,cfg%csfdim,hdiag,&
             guessscr,vecscr,vecstem,nroots,blocksize,maxvec,ldeflate,&
             tol,niter)
     end select
        
  endif

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  call cfg%finalise
  deallocate(hdiag,averageii)
  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine diag_mrci
