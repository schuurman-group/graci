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
subroutine diag_dftcis(irrep,nroots,vecscr) &
     bind(c,name='diag_dftcis')
#else
subroutine diag_dftcis(irrep,nroots,vecscr)
#endif

  use constants
  use bitglobal
  use dftcis_param
  use dftcis_conf
  use dftcis_hbuild
  use full_diag
  use iomod
  
  implicit none

  ! Irrep number
  integer(is), intent(in)  :: irrep

  ! No. roots
  integer(is), intent(in)  :: nroots

  ! Eigenpair scratch file number
  integer(is), intent(out) :: vecscr
  
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
  integer(is), parameter   :: full_lim=2500

  ! Subspace dimension to be used in the generation of the
  ! guess vectors
  integer(is)              :: subdim

  ! Guess vector scratch file number
  integer(is)              :: guessscr
  
  ! Everything else
  integer(is)              :: i
  integer(is)              :: iham
  integer(is)              :: blocksize,maxvec

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  write(6,'(/,52a)') ('-',i=1,52)
  write(6,'(3(x,a))') 'DFT/CIS calculation in the',&
       trim(irreplbl(irrep,ipg)),'subspace'
  write(6,'(52a)') ('-',i=1,52)
  
!----------------------------------------------------------------------
! Checks on the no. open shells/multiplicity
!----------------------------------------------------------------------
  ! No. open shells in the base configuration
  if (mod(nel,2) /= 0) then
     errmsg='Error in diag_dftcis: ROCIS calculations are not yet'&
          //' supported'
     call error_control
  endif

!----------------------------------------------------------------------
! Load the DFT/MRCI Hamiltonian parameters
!----------------------------------------------------------------------
  ! Temporary hardwiring of the Hamiltonian index
  iham=2
  
  ! Load the parameters
  call load_hpar_dftcis(iham)
    
!----------------------------------------------------------------------
! Generate the CIS configuration information
!----------------------------------------------------------------------
  call generate_cis_confs(irrep,ncsf,iph)
  
  write(6,'(/,x,a,x,i0)') 'N_CSF:',ncsf
  
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
  call save_hij_dftcis(irrep,hscr,nrec,ncsf,iph)
  
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
     subdim=1500
     if (subdim > ncsf) subdim=int(ncsf*0.9d0)

     print*,'The iterative DFT/CIS diagonalisation interface needs writing...'
     stop
     
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
! ras_guess_dftcis: Returns a suggestion for RAS1 and RAS2 MO indices
!                   based on the dominant CSFs in a set of DFT/CIS
!                   eigenfunctions
!######################################################################
#ifdef CBINDING
subroutine ras_guess_dftcis(irrep,nroots,vecscr,domph) &
     bind(c,name='ras_guess_dftcis')
#else
subroutine ras_guess_dftcis(irrep,nroots,vecscr,domph)
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
  real(dp), parameter      :: phthrsh=0.3d0
  
  ! Everything else
  integer(is)              :: k,n
  integer(is), allocatable :: indx(:)
  real(dp), allocatable    :: abscoe(:)
  
!----------------------------------------------------------------------
! Generate the CIS configuration information
!----------------------------------------------------------------------
  call generate_cis_confs(irrep,ncsf,iph)

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

     ! Loop over CSFs
     do n=1,ncsf

        ! Save the above threshold CSF information
        if (abs(vec(indx(n),k)) >= phthrsh) then

           ! Particle index
           domph(iph(1,indx(n)))=1

           ! Hole index
           domph(iph(2,indx(n)))=1
           
        endif

     enddo
     
  enddo
  
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
