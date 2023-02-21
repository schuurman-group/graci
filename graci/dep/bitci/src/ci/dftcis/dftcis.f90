!**********************************************************************
! Semi-empirical DFT/CIS calculations.
! See Chem. Phys. Lett., 259, 128 (1996) for details of the method.
!**********************************************************************
! For now, we will just use a disk-based algorithm in which the
! non-zero Hamiltonian matrix elements are pre-computed and written
! to disk. This way, the MRCI Davidson module can be used.
!**********************************************************************

!######################################################################
! diag_dftcis: Controls the diagonalisation of the DFT/CIS Hamiltonian
!              for a given symmetry subspace
!######################################################################
#ifdef CBINDING
subroutine diag_dftcis(irrep,nroots,icvs,vecscr,loose,iham) &
     bind(c,name='diag_dftcis')
#else
subroutine diag_dftcis(irrep,nroots,icvs,vecscr,loose,iham)
#endif

  use constants
  use bitglobal
  use dftcis_param
  use dftcis_conf
  use dftcis_hbuild
  use dftcis_guess
  use full_diag
  use gendav
  use iomod
  use conftype
  
  implicit none

  ! Irrep number
  integer(is), intent(in)  :: irrep

  ! No. roots
  integer(is), intent(in)  :: nroots

  ! CVS-MRCI: core MOs
  integer(is), intent(in)    :: icvs(nmo)
  
  ! Eigenpair scratch file number
  integer(is), intent(out) :: vecscr

  ! Loose integral screening: only to be used for the generation
  ! of RAS spaces
  logical, intent(in)      :: loose

  ! DFT/CIS Hamiltonian index
  integer(is), intent(in)  :: iham
  
  ! No. CIS CSFs
  integer(is)              :: ncsf
  
  ! CIS particle/hole indices
  integer(is), allocatable :: iph(:,:)

  ! Hamiltonian scratch file number
  integer(is)              :: hscr

  ! No. records on the Hamiltonian scratch file
  integer(is)              :: nrec

  ! On-diagonal Hamiltonian matrix elements
  real(dp), allocatable    :: hii(:)

  ! Dimension of the CSF basis past which full diagonalisation
  ! will not be used
  integer(is), parameter   :: full_lim=750

  ! Guess vector scratch file number
  integer(is)              :: guessscr

  ! Guess vector subspace dimension
  integer(is)              :: guessdim
  
  ! Davidson variables
  integer(is)              :: blocksize,maxvec,ipre,niter
  real(dp)                 :: tol
  
  ! Everything else
  integer(is)              :: i
  integer(is)              :: isigma(3)

  ! Dummy MRCI configuration derived type: temporarily required
  ! to be passed to the Davidson routines
  type(mrcfg)              :: cfg

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  if (verbose) then
     write(6,'(/,52a)') ('-',i=1,52)
     write(6,'(3(x,a))') 'DFT/CIS calculation in the',&
          trim(irreplbl(irrep,ipg)),'subspace'
     write(6,'(52a)') ('-',i=1,52)
  endif
     
!----------------------------------------------------------------------
! Checks on the no. open shells/multiplicity
!----------------------------------------------------------------------
  ! No. open shells in the base configuration
  if (mod(nel,2) /= 0) then
     errmsg='Error in diag_dftcis: ROCIS calculations are not yet'&
          //' supported'
     call error_control
  endif

  ! Triplets not yet coded up
  if (imult == 3) then
     errmsg='Error in diag_dftcis: triplet calculations are not yet'&
          //' supported'
     call error_control
  endif
  
!----------------------------------------------------------------------
! Load the DFT/MRCI Hamiltonian parameters
!----------------------------------------------------------------------
  ! Load the parameters
  call load_hpar_dftcis(iham)
    
!----------------------------------------------------------------------
! Generate the CIS configuration information
!----------------------------------------------------------------------
  call generate_cis_confs(irrep,ncsf,icvs,iph)
  
  if (verbose) write(6,'(/,x,a,x,i0)') 'N_CSF:',ncsf
  
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(hii(ncsf))
  hii=0.0d0

!----------------------------------------------------------------------
! Calculate the on-diagonal elements of the DFT/CIS Hamiltonian
! matrix
!----------------------------------------------------------------------
  call hii_dftcis(ncsf,iph,hii)
  
!----------------------------------------------------------------------
! Calculate and save the non-zero off-diagonal elements of the DFT/CIS
! Hamiltonian matrix
!----------------------------------------------------------------------
  call save_hij_dftcis(irrep,hscr,nrec,ncsf,iph,loose)
  
!----------------------------------------------------------------------
! Full diagonalisation
!----------------------------------------------------------------------
  if (ncsf <= full_lim) call diag_full(hscr,nrec,ncsf,hii,irrep,&
       nroots,vecscr,'dftcis')

!----------------------------------------------------------------------
! Iterative diagonalisation
!----------------------------------------------------------------------
  if (ncsf > full_lim) then

     ! Temporary hard wiring of parameters
     blocksize=nroots*2
     maxvec=4*blocksize
     niter=100
     tol=5e-4_dp

     ! Dimension of the subspace used to generate the guess vectors
     guessdim=750
     if (guessdim > ncsf) guessdim=int(ncsf*0.9d0)

     ! Guess vector generation
     call dftcis_guess_subspace(ncsf,guessdim,hii,iph,blocksize,&
          guessscr)
     
     ! Set the sigma-vector algorithm information: disk-based only
     ! for now
     isigma(1)=1
     isigma(2)=hscr
     isigma(3)=nrec

     ! Set the preconditioner: DPR makes sense due to the small
     ! no. CSFs in a DFT/CIS calculation
     ipre=1
     
     ! Generalised Davidson diagonalisation
     call generalised_davidson(irrep,isigma,cfg,ncsf,hii,guessscr,&
          vecscr,'dftcisvec',nroots,blocksize,maxvec,tol,&
          niter,ipre)
     
  endif
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  ! Unload the DFT/MRCI Hamiltonian parameter arrays
  call unload_hpar_dftcis
  
  deallocate(hii)
  
  return
  
end subroutine diag_dftcis

!######################################################################
! ras_guess_dftcis: Returns a suggestion for RAS1 and RAS3 MO indices
!                   based on the dominant CSFs in a set of DFT/CIS
!                   eigenfunctions
!######################################################################
#ifdef CBINDING
subroutine ras_guess_dftcis(irrep,nroots,icvs,vecscr,domph) &
     bind(c,name='ras_guess_dftcis')
#else
subroutine ras_guess_dftcis(irrep,nroots,icvs,vecscr,domph)
#endif

  use constants
  use bitglobal
  use dftcis_conf
  use utils
  use iomod

  implicit none

  ! Irrep number
  integer(is), intent(in)  :: irrep

  ! No. roots
  integer(is), intent(in)  :: nroots

  ! CVS-MRCI: core MOs
  integer(is), intent(in)  :: icvs(nmo)
  
  ! Eigenpair scratch file number
  integer(is), intent(out) :: vecscr

  ! Dominant particle/hole indices
  integer(is), intent(out) :: domph(nmo)
  
  ! No. CIS CSFs
  integer(is)              :: ncsf
  
  ! CIS particle/hole indices
  integer(is), allocatable :: iph(:,:)

  ! Eigenpairs
  real(dp), allocatable    :: vec(:,:),ener(:)

  ! MO selection threshold
  real(dp), parameter      :: targ=0.9025d0

  ! Limit on the size of the RAS spaces
  integer(is), parameter   :: maxras=60
  
  ! Everything else
  integer(is)              :: k,n,np,nh,ihomo
  integer(is), allocatable :: indx(:)
  real(dp), allocatable    :: abscoe(:)
  real(dp)                 :: sumsq
  
!----------------------------------------------------------------------
! Generate the CIS configuration information
!----------------------------------------------------------------------
  call generate_cis_confs(irrep,ncsf,icvs,iph)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(vec(ncsf,nroots))
  vec=0.0d0

  allocate(ener(nroots))
  ener=0.0d0

  allocate(indx(ncsf))
  indx=0

  allocate(abscoe(ncsf))
  abscoe=0.0d0
  
!----------------------------------------------------------------------
! Read in the eigenpairs
!----------------------------------------------------------------------
  call read_all_eigenpairs(vecscr,vec,ener,ncsf,nroots)
  
!----------------------------------------------------------------------
! Determine the particle/hole indices for the dominant CSFs
!----------------------------------------------------------------------
  ! Initialisation
  domph=0

  ! Loop over roots
  do k=1,nroots

     ! Sort the CSFs
     abscoe=abs(vec(:,k))
     call dsortindxa1('D',ncsf,abscoe,indx)

     sumsq=0.0d0
     
     ! Loop over CSFs
     do n=1,ncsf

        ! Squared norm
        sumsq=sumsq+vec(indx(n),k)**2
        
        ! Particle index
        domph(iph(1,indx(n)))=1
        
        ! Hole index
        domph(iph(2,indx(n)))=1

        ! Exit if we have reached the squared norm target
        if (sumsq >= targ) exit
        
     enddo
     
  enddo

!----------------------------------------------------------------------
! Ensure that the RAS spaces are reasonably sized
!----------------------------------------------------------------------
  ! Base conf HOMO index (note that we are taking M=S)
  ihomo=nel_alpha

  ! Trim the particle MOs  
  np=sum(domph(ihomo+1:nmo))
  if (np > maxras) then
     n=np
     do k=nmo,ihomo+1,-1
        if (domph(k) == 1) then
           domph(k)=0
           n=n-1
           if (n == maxras) exit
        endif
     enddo
  endif

  ! Trim the hole MOs, excluding any core ones flagged in the icvs
  ! array
  nh=sum(domph(1:ihomo))
  if (nh > maxras) then
     n=nh
     do k=1,ihomo
        if (domph(k) == 1 .and. icvs(k) /= 1) then
           domph(k)=0
           n=n-1
           if (n == maxras) exit
        endif
     enddo
  endif
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(vec)
  deallocate(ener)
  deallocate(indx)
  deallocate(abscoe)
  
  return
  
end subroutine ras_guess_dftcis

!######################################################################
