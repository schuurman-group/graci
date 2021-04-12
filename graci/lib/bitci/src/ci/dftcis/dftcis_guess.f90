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
  
end module dftcis_guess
