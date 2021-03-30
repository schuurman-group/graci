!**********************************************************************
! Sigma-vector routines
!**********************************************************************
module sigma_mrci

contains

!######################################################################
! sigma_disk: Calculation of a batch of sigma vectors using
!             off-diagonal Hamiltonian matrix elements saved to disk
!######################################################################
  subroutine sigma_disk(hamscr,nrec,matdim,blocksize,bvec,sigvec,hdiag)

    use constants
    use bitglobal
    use iomod
    
    implicit none

    ! Hamiltonian scratch file no. and no. records
    integer(is), intent(in)  :: hamscr,nrec

    ! Dimensions
    integer(is), intent(in)  :: matdim,blocksize

    ! Subspace vectors
    real(dp), intent(in)     :: bvec(matdim,blocksize)

    ! Sigma vectors
    real(dp), intent(out)    :: sigvec(matdim,blocksize)

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: hdiag(matdim)

    ! Working arrays
    real(dp), allocatable    :: hbuffer(:)
    integer(is), allocatable :: ibuffer(:,:)

    ! Everything else
    integer(is)              :: iscratch,irec,nbuf,m,n,i,j
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(ibuffer(2,bufsize))
    ibuffer=0
    
    allocate(hbuffer(bufsize))
    hbuffer=0.0d0
    
!----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!----------------------------------------------------------------------
    sigvec=0.0d0
    do m=1,matdim
       sigvec(m,:)=hdiag(m)*bvec(m,:)
    enddo
    
!----------------------------------------------------------------------
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!----------------------------------------------------------------------
    ! Open the scratch file
    iscratch=scrunit(hamscr)
    open(iscratch,file=scrname(hamscr),form='unformatted',status='old')

    ! Loop over records
    do irec=1,nrec
    
       ! Read in the next record
       read(iscratch) hbuffer,ibuffer,nbuf

       ! Loop over elements in the buffer
       do n=1,nbuf

          ! Hij
          i=ibuffer(1,n)
          j=ibuffer(2,n)

          ! Contribution to the sigma vectors
          sigvec(i,:)=sigvec(i,:)+hbuffer(n)*bvec(j,:)
          sigvec(j,:)=sigvec(j,:)+hbuffer(n)*bvec(i,:)
                       
       enddo
       
    enddo
    
    ! Close the scratch file
    close(iscratch)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ibuffer)
    deallocate(hbuffer)

    return
    
  end subroutine sigma_disk
  
!######################################################################
! sigma_disk_old: Calculation of a batch of sigma vectors using
!                 off-diagonal Hamiltonian matrix elements saved to
!                 disk
!######################################################################
  subroutine sigma_disk_old(hamscr,nrec,matdim,maxvec,currdim,wmat,&
       vmat,hdiag)

    use constants
    use bitglobal
    use iomod
    
    implicit none

    ! Hamiltonian scratch file no. and no. records
    integer(is), intent(in)  :: hamscr
    integer(is), intent(in)  :: nrec

    ! Hamiltonian matrix dimension
    integer(is), intent(in)  :: matdim

    ! Maximum subspace dimension
    integer(is), intent(in)  :: maxvec
    
    ! Current subspace dimension
    integer(is), intent(in)  :: currdim

    ! Matrix of subspace vectors
    real(dp), intent(in)     :: vmat(matdim,maxvec)
      
    ! Matrix-vector products
    real(dp), intent(out)    :: wmat(matdim,maxvec)

    ! On-diagonal matrix elements
    real(dp), intent(in)     :: hdiag(matdim)

    ! Working arrays
    real(dp), allocatable    :: hbuffer(:)
    integer(is), allocatable :: ibuffer(:,:)

    ! Everything else
    real(dp)                 :: hij
    integer(is)              :: iscratch,irec,nbuf,m,n,i,j

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(ibuffer(2,bufsize))
    allocate(hbuffer(bufsize))
    
!----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!----------------------------------------------------------------------
    wmat=0.0d0
    do n=1,currdim
       do m=1,matdim
          wmat(m,n)=hdiag(m)*vmat(m,n)
       enddo
    enddo
    
!----------------------------------------------------------------------
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!----------------------------------------------------------------------
    ! Open the scratch file
    iscratch=scrunit(hamscr)
    open(iscratch,file=scrname(hamscr),form='unformatted',status='old')

    ! Loop over records
    do irec=1,nrec
    
       ! Read in the next record
       read(iscratch) hbuffer,ibuffer,nbuf

       ! Loop over elements in the buffer
       do n=1,nbuf
          
          ! Hij
          i=ibuffer(1,n)
          j=ibuffer(2,n)
          hij=hbuffer(n)

          ! Loop over vectors in the subspace
          wmat(i,1:currdim)=wmat(i,1:currdim)+hij*vmat(j,1:currdim)
          wmat(j,1:currdim)=wmat(j,1:currdim)+hij*vmat(i,1:currdim)
          
       enddo
       
    enddo
    
    ! Close the scratch file
    close(iscratch)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ibuffer)
    deallocate(hbuffer)

    return
    
  end subroutine sigma_disk_old

!######################################################################
  
end module sigma_mrci
