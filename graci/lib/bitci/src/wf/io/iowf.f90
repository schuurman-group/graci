!**********************************************************************
! I/O routines specific to the bitwf library, i.e., for the parsing
! of bitwf determinant scratch files
!**********************************************************************
module iowf

  implicit none

contains

!######################################################################
! read_ndet: reads the number of determinants from the scratch file
!            numbered wfscr
!######################################################################
  subroutine read_ndet(wfscr,ndet)

    use constants
    use bitglobal
    
    implicit none

    ! Scratch file number
    integer(is), intent(in)  :: wfscr

    ! No. determinants
    integer(is), intent(out) :: ndet

    ! Everything else
    integer(is)              :: iscratch,idum

    !
    ! Open the scratch file
    !
    iscratch=scrunit(wfscr)
    open(iscratch,file=scrname(wfscr),form='unformatted',status='old')

    !
    ! Read past the spin multiplicity and n_int
    !
    read(iscratch) idum
    read(iscratch) idum
    
    !
    ! Read the no. determinants
    !
    read(iscratch) ndet
    
    !
    ! Close the scratch file
    !
    close(iscratch)
    
    return
    
  end subroutine read_ndet
  
!######################################################################
! read_detwf: reads the determinant wave function information
!             (eigenvectors and determinant bit strings) from the
!             scratch file numbered wfscr
!
!             only the subset of eigenvectors specified in the iroots
!             array are read
!######################################################################
  subroutine read_detwf(wfscr,ndet,nvec,n_intX,det,vec,iroots)

    use constants
    use bitglobal
    use iomod
    
    implicit none

    ! Scratch file number
    integer(is), intent(in)  :: wfscr

    ! Dimensions
    integer(is), intent(in)  :: ndet,nvec,n_intX

    ! Determinants
    integer(ib), intent(out) :: det(n_intX,2,ndet)

    ! Eigenvectors
    real(dp), intent(out)    :: vec(ndet,nvec)

    ! Requested eigenvectors
    integer(is), intent(in)  :: iroots(nvec)
    
    ! Everything else
    integer(is)              :: iscratch,idum,i,n,nentries
    integer(is)              :: ipos(1)
    real(dp), allocatable    :: vec1(:)
    
    !
    ! Allocate working arrays
    !
    allocate(vec1(ndet))
    
    !
    ! Open the scratch file
    !
    iscratch=scrunit(wfscr)
    open(iscratch,file=scrname(wfscr),form='unformatted',status='old')
    
    !
    ! Read past the spin multiplicity
    !
    read(iscratch) idum

    !
    ! Sanity check on dimensions
    !
    ! n_intB/K
    read(iscratch) idum
    if (idum /= n_intX) then
       errmsg='Error in read_detwf: incorrect value of n_intX'
       call error_control
    endif
    
    ! ndet
    read(iscratch) idum
    if (idum /= ndet) then
       errmsg='Error in read_detwf: incorrect value of ndet'
       call error_control
    endif
    
    !
    ! Total no. roots saved to disk
    !
    read(iscratch) nentries

    !
    ! Read past the no. MOs
    !
    read(iscratch) idum

    !
    ! Eigenvectors
    !
    n=0
    do i=1,nentries
       ! Read in the next eigenvector
       read(iscratch) vec1(1:ndet)
       ! Save the eigenvector if it is in the subset of requested roots
       ipos=findloc(iroots,value=i)
       if (ipos(1) /= 0) then
          n=n+1
          vec(1:ndet,n)=vec1(1:ndet)
       endif
    enddo
    
    !
    ! Determinant bit strings
    !
    read(iscratch) det

    !
    ! Close the scratch file
    !
    close(iscratch)

    return
    
  end subroutine read_detwf

!######################################################################
  
end module iowf
