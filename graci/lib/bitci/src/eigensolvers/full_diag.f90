!**********************************************************************
! Routines for the full diagonalisation of the MRCI Hamiltonian
!**********************************************************************
module full_diag

  implicit none

contains

!######################################################################
! diag_mrci_full: Full diagonalisation of the MRCI Hamiltonian
!######################################################################
  subroutine diag_full(hamscr,nrec,csfdim,hdiag,irrep,nroots,&
       vecscr,vecstem)

    use constants
    use bitglobal
    use conftype
    use utils
    use iomod
    
    implicit none

    ! Hamiltonian scratch file number
    integer(is), intent(in)      :: hamscr

    ! Number of records on in the Hamiltonian scratch file
    integer(is), intent(in)      :: nrec

    ! CSF basis dimension
    integer(is), intent(in)      :: csfdim

    ! On-diagonal matrix elements
    real(dp), intent(in)         :: hdiag(csfdim)

    ! Irrep number
    integer(is), intent(in)      :: irrep

    ! Number of roots requested
    integer(is), intent(in)      :: nroots
    
    ! Eigenpair scratch file number and file stem
    integer(is), intent(out)     :: vecscr
    character(len=*), intent(in) :: vecstem
    
    ! Hamiltonian matrix and eigenpairs
    real(dp), allocatable        :: hmat(:,:),heig(:),hvec(:,:)

    ! I/O variables
    integer(is)                  :: iscratch
    character(len=60)            :: vecfile
    character(len=2)             :: amult,airrep
    
    ! Everthing else
    real(dp), allocatable        :: hbuffer(:)
    integer(is), allocatable     :: ibuffer(:,:)
    integer(is)                  :: irec,nbuf,n,i,j
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(hmat(csfdim,csfdim))
    hmat=0.0d0

    allocate(heig(csfdim))
    heig=0.0d0
    
    allocate(hvec(csfdim,csfdim))
    hvec=0.0d0

    allocate(ibuffer(2,bufsize))
    ibuffer=0
    
    allocate(hbuffer(bufsize))
    hbuffer=0.0d0

!----------------------------------------------------------------------
! Fill in the on-diagonal matrix elements
!----------------------------------------------------------------------
    do i=1,csfdim
       hmat(i,i)=hdiag(i)
    enddo
    
!----------------------------------------------------------------------
! Read in the non-zero off-diagonal elements from disk
!----------------------------------------------------------------------
    ! Open the scratch file
    iscratch=scrunit(hamscr)
    open(iscratch,file=scrname(hamscr),form='unformatted',status='old')

    ! Loop over records
    do irec=1,nrec

       ! Read in the next record
       read(iscratch) hbuffer,ibuffer,nbuf
       
       ! Fill in the Hamiltonian matrix
       do n=1,nbuf
          i=ibuffer(1,n)
          j=ibuffer(2,n)
          hmat(i,j)=hbuffer(n)
          hmat(j,i)=hbuffer(n)
       enddo
       
    enddo
    
    ! Close the scratch file
    close(iscratch)

!----------------------------------------------------------------------
! Diagonalise the reference space Hamiltonian matrix
!----------------------------------------------------------------------
    ! Diagonalisation
    call diag_matrix_real(hmat,heig,hvec,csfdim)
  
    ! Add on E_SCF + E_nuc
    heig=heig+escf+enuc

!----------------------------------------------------------------------
! Save the eigenpairs to disk
!----------------------------------------------------------------------
    ! Register the scratch file
    write(amult,'(i0)') imult
    write(airrep,'(i0)') irrep
    call scratch_name(trim(vecstem)//'.mult'//trim(amult)//&
         '.sym'//trim(airrep),vecfile)
    call register_scratch_file(vecscr,vecfile)
    
    ! Open the scratch file
    iscratch=scrunit(vecscr)
    open(iscratch,file=scrname(vecscr),form='unformatted',&
         status='unknown')
    
    ! No. CSFs
    write(iscratch) csfdim
    
    ! No. roots
    write(iscratch) nroots
    
    ! Eigenvalues
    write(iscratch) heig(1:nroots)
    
    ! Eigenvectors
    do i=1,nroots
       write(iscratch) hvec(:,i)
    enddo
    
    ! Close the scratch file
    close(iscratch)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(hmat)
    deallocate(heig)
    deallocate(hvec)
    deallocate(ibuffer)
    deallocate(hbuffer)
    
    return
    
  end subroutine diag_full
    
!######################################################################

end module full_diag
