!**********************************************************************
! Module for the calculation of ref space guess vectors for use in
!  an interative diagonalisation calculation
!**********************************************************************
module ref_guess

  implicit none
  
contains

!######################################################################
! ref_guess_subspace: Generation of ref space guess vectors via
!                     a subspace diagonalisation
!######################################################################
  subroutine ref_guess_subspace(guessscr,nguess,dim,subdim,hii,nrec,&
       hscr)

    use constants
    use bitglobal
    use utils
    use iomod
    
    implicit none

    ! Guess vector scratch file number
    integer(is), intent(out) :: guessscr

    ! Number of guess vectors
    integer(is), intent(in)  :: nguess
    
    ! Full space and subspace dimensions
    integer(is), intent(in)  :: dim,subdim

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: hii(dim)

    ! Hamiltonian scratch file no. and no. records
    integer(is), intent(in)  :: hscr,nrec

    ! Indices of the diagonals of the full Hamiltonian in
    ! ascending order
    integer(is), allocatable :: indx(:)

    ! Working arrays
    integer(is), allocatable :: isub(:),f2s(:),s2f(:)

    ! Subspace Hamiltonian matrix
    real(dp), allocatable    :: subhmat(:,:)

    ! Subspace eigenpairs
    real(dp), allocatable    :: subvec(:,:),subeig(:)

    ! Guess vectors
    real(dp), allocatable    :: guessvec(:,:)
    
    ! I/O variables
    integer(is)              :: iscratch
    character(len=250)       :: guessfile
    character(len=2)         :: amult,airrep
    
    ! Everthing else
    real(dp), allocatable    :: hbuffer(:)
    integer(is), allocatable :: ibuffer(:,:)
    integer(is)              :: irec,nbuf,n,i,j,i1,j1
    integer(is)              :: counter

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(indx(dim))
    indx=0

    allocate(isub(dim))
    isub=0

    allocate(f2s(dim))
    f2s=0

    allocate(s2f(subdim))
    s2f=0

    allocate(subhmat(subdim,subdim))
    subhmat=0.0d0

    allocate(subvec(subdim,subdim))
    subvec=0.0d0

    allocate(subeig(subdim))
    subeig=0.0d0

    allocate(guessvec(dim,nguess))
    guessvec=0.0d0
    
    allocate(ibuffer(2,bufsize))
    ibuffer=0
    
    allocate(hbuffer(bufsize))
    hbuffer=0.0d0
    
!----------------------------------------------------------------------
! Sort the on-diagonal Hamiltonian matrix element in order of
! increasing value
!----------------------------------------------------------------------
    call dsortindxa1('A',dim,hii,indx)

!----------------------------------------------------------------------
! Flag the subspace CSFs
!----------------------------------------------------------------------
    do i=1,subdim
       isub(indx(i))=1
    enddo
    
!----------------------------------------------------------------------
! Full-to-subspace space mapping
!----------------------------------------------------------------------
    counter=0
    do i=1,dim
       if (isub(i) == 1) then
          counter=counter+1
          f2s(i)=counter
          s2f(counter)=i
       endif
    enddo
    
!----------------------------------------------------------------------
! On-diagonal elements
!----------------------------------------------------------------------
    do i=1,dim
       if (isub(i) == 1) then
          i1=f2s(i)
          subhmat(i1,i1)=hii(i)
       endif
    enddo
    
!----------------------------------------------------------------------
! Off-diagonal elements
!----------------------------------------------------------------------
    ! Open the scratch file
    iscratch=scrunit(hscr)
    open(iscratch,file=scrname(hscr),form='unformatted',status='old')

    ! Loop over records
    do irec=1,nrec

       ! Read in the next record
       read(iscratch) hbuffer,ibuffer,nbuf
       
       ! Fill in the subspace Hamiltonian matrix
       do n=1,nbuf
          
          i=ibuffer(1,n)
          j=ibuffer(2,n)

          if (isub(i) == 1 .and. isub(j) == 1) then
             i1=f2s(i)
             j1=f2s(j)
             subhmat(i1,j1)=hbuffer(n)
             subhmat(j1,i1)=subhmat(i1,j1)
          endif
             
       enddo
       
    enddo
    
    ! Close the scratch
    close(iscratch)

!----------------------------------------------------------------------
! Diagonalise the reference space Hamiltonian matrix
!----------------------------------------------------------------------
    ! Diagonalisation
    call diag_matrix_real(subhmat,subeig,subvec,subdim)

!----------------------------------------------------------------------
! Construct the guess vectors
!----------------------------------------------------------------------
    ! Loop over guess vectors
    do i=1,nguess

       ! Loop over subspace CSFs
       do j=1,subdim

          ! Fill in the guess vector
          j1=s2f(j)
          guessvec(j1,i)=subvec(j,i)
          
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
    do i=1,nguess
       write(iscratch) guessvec(:,i)
    enddo

    ! Close the scratch file
    close(iscratch)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(indx)
    deallocate(isub)
    deallocate(f2s)
    deallocate(s2f)
    deallocate(subhmat)
    deallocate(subvec)
    deallocate(subeig)
    deallocate(guessvec)
    deallocate(ibuffer)
    deallocate(hbuffer)
    
    return
    
  end subroutine ref_guess_subspace

!######################################################################
  
end module ref_guess
