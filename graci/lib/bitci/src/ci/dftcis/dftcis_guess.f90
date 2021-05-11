!**********************************************************************
! Routines for the generation of DFT/CIS guess vectors
!**********************************************************************

module dftcis_guess

contains

!######################################################################
! dftcis_guess_unit: Generation of single CSF guess vectors
!######################################################################
  subroutine dftcis_guess_unit(ncsf,hii,nguess,guessscr)

    use constants
    use bitglobal
    use utils
    use iomod
    
    implicit none

    ! No. CSFs
    integer(is), intent(in)  :: ncsf

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: hii(ncsf)

    ! No. guess vectors
    integer(is), intent(in)  :: nguess

    ! Guess vector scratch file number
    integer(is), intent(out) :: guessscr
    
    ! Guess vectors
    real(dp), allocatable    :: guessvec(:,:)
    
    ! Working arrays
    integer(is), allocatable :: indx(:)

    ! I/O variables
    integer(is)              :: iscratch
    character(len=60)        :: guessfile
    
    ! Everything else
    integer(is)              :: i
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(indx(ncsf))
    indx=0

    allocate(guessvec(ncsf,nguess))
    guessvec=0.0d0
    
!----------------------------------------------------------------------
! Sort the on-diagonal Hamiltonian matrix elements in order of
! increasing value
!----------------------------------------------------------------------
    call dsortindxa1('A',ncsf,hii,indx)

!----------------------------------------------------------------------
! Construct the guess vectors
!----------------------------------------------------------------------
    do i=1,nguess
       guessvec(indx(i),i)=1.0d0
    enddo

!----------------------------------------------------------------------
! Register the guess vector scratch file
!----------------------------------------------------------------------
    call scratch_name('guessvecs',guessfile)
    call register_scratch_file(guessscr,guessfile)
    
!----------------------------------------------------------------------
! Write the guess vectors to disk
!----------------------------------------------------------------------
    ! Open the scratch file
    iscratch=scrunit(guessscr)
    open(iscratch,file=scrname(guessscr),form='unformatted',&
         status='unknown')
    
    ! Dimensions
    write(iscratch) ncsf
    write(iscratch) nguess
    
    ! Guess vectors
    do i=1,nguess
       write(iscratch) guessvec(:,i)
    enddo
    
    ! Close the scratch file
    close(iscratch)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(indx)
    deallocate(guessvec)
    
    return
    
  end subroutine dftcis_guess_unit
  
!######################################################################
! dftcis_guess_subspace: Generation of guess vectors via subspace
!                        diagonalisation
!######################################################################
  subroutine dftcis_guess_subspace(dim,subdim,hii,iph,nguess,guessscr)

    use constants
    use bitglobal
    use dftcis_hbuild, only: hij_dftcis
    use utils
    use iomod
    
    implicit none

    ! No. CSFs in the full and subspaces
    integer(is), intent(in)  :: dim,subdim

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: hii(dim)

    ! CIS particle/hole indices
    integer(is), intent(in)  :: iph(2,dim)
    
    ! No. guess vectors
    integer(is), intent(in)  :: nguess

    ! Guess vector scratch file number
    integer(is), intent(out) :: guessscr

    ! Indices of the diagonals of the full Hamiltonian in
    ! ascending order
    integer(is), allocatable :: indx(:)

    ! Working arrays
    integer(is), allocatable :: isub(:)

    ! Subspace Hamiltonian matrix
    real(dp), allocatable    :: subhmat(:,:)

    ! Subspace eigenpairs
    real(dp), allocatable    :: subvec(:,:),subeig(:)

    ! Guess vectors
    real(dp), allocatable    :: guessvec(:,:)

    ! I/O variables
    integer(is)              :: iscratch
    character(len=60)        :: guessfile
    
    ! Everything else
    integer(is)              :: bcsf,kcsf,icsf,n
    integer(is)              :: i,j,a,b

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(indx(dim))
    indx=0

    allocate(isub(dim))
    isub=0

    allocate(subhmat(subdim,subdim))
    subhmat=0.0d0

    allocate(subvec(subdim,subdim))
    subvec=0.0d0

    allocate(subeig(subdim))
    subeig=0.0d0

    allocate(guessvec(dim,nguess))
    guessvec=0.0d0
    
!----------------------------------------------------------------------
! Sort the on-diagonal matrix elements in order of increasing value
!----------------------------------------------------------------------
    call dsortindxa1('A',dim,hii,indx)

!----------------------------------------------------------------------
! Flag the subspace CSFs
!----------------------------------------------------------------------
    do i=1,subdim
       isub(indx(i))=1
    enddo

!----------------------------------------------------------------------
! On-diagonal elements of the subspace Hamiltonian
!----------------------------------------------------------------------
    do i=1,subdim
       subhmat(i,i)=hii(indx(i))
    enddo

!----------------------------------------------------------------------
! Off-diagonal elements of the subspace Hamiltonian
!----------------------------------------------------------------------
    ! Loop over bra CSFs
    do bcsf=1,subdim-1

       ! Bra particle and hole indices
       a=iph(1,indx(bcsf))
       i=iph(2,indx(bcsf))

       ! Loop over ket CSFs
       do kcsf=bcsf+1,subdim

          ! Ket particle and hole indices
          b=iph(1,indx(kcsf))
          j=iph(2,indx(kcsf))

          ! Matrix element
          subhmat(bcsf,kcsf)=hij_dftcis(i,a,j,b)
          subhmat(kcsf,bcsf)=subhmat(bcsf,kcsf)
          
       enddo

    enddo

!----------------------------------------------------------------------
! Diagonalise the subspace Hamiltonian
!----------------------------------------------------------------------
    call diag_matrix_real(subhmat,subeig,subvec,subdim)

!----------------------------------------------------------------------
! Construct the guess vectors
!----------------------------------------------------------------------
    ! Loop over guess vectors
    do n=1,nguess

       ! Loop over subspace CSFs
       do icsf=1,subdim

          guessvec(indx(icsf),n)=subvec(icsf,n)
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Register the guess vector scratch file
!----------------------------------------------------------------------
    call scratch_name('guessvecs',guessfile)
    call register_scratch_file(guessscr,guessfile)

!----------------------------------------------------------------------
! Write the guess vectors to disk
!----------------------------------------------------------------------
    ! Open the scratch file
    iscratch=scrunit(guessscr)
    open(iscratch,file=scrname(guessscr),form='unformatted',&
         status='unknown')

    ! Dimensions
    write(iscratch) dim
    write(iscratch) nguess

    ! Guess vectors
    do n=1,nguess
       write(iscratch) guessvec(:,n)
    enddo

    ! Subspace eigenpairs (needed if the generalised Davidson
    ! eigensolver is being used)
    write(iscratch) subdim
    write(iscratch) indx(1:subdim)
    write(iscratch) subeig
    write(iscratch) subvec
    
    ! Close the scratch file
    close(iscratch)
    
    return

  end subroutine dftcis_guess_subspace

!######################################################################
  
end module dftcis_guess
