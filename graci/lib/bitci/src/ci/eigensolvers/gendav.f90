!**********************************************************************
! Generalised Davidson diagonalisation module
! Adheres to the notation used in J. Comput. Chem. 22, 1574 (2001)
!**********************************************************************
module gendav

  use constants
  
  implicit none

  ! Subspace vectors
  real(dp), allocatable    :: bvec(:,:)

  ! Sigma vectors
  real(dp), allocatable    :: sigvec(:,:)

  ! Subspace Hamiltonian matrix and eigenpairs
  real(dp), allocatable    :: Gmat(:,:),alpha(:,:),rho(:),rho1(:)

  ! Residual norms
  real(dp), allocatable    :: rnorm(:)

  ! Generalised davidson preconditioner arrays
  integer(is)              :: subdim
  integer(is), allocatable :: s2f(:)
  integer(is), allocatable :: idiag(:)
  real(dp), allocatable    :: subeig(:),subvec(:,:)
  real(dp), allocatable    :: work(:)
  
  ! Counters, etc.
  integer(is)              :: currdim,nconv,nnew,nsigma
  integer(is), allocatable :: iconv(:)
  
contains

!######################################################################
! generalised_davidson: Master routine for a generalised Davidson
!                       calculation
!######################################################################
  subroutine generalised_davidson(irrep,isigma,cfg,matdim,hdiag,&
       guessscr,vecscr,vecstem,nstates,blocksize,maxvec,tol,niter,&
       ipre)

    use constants
    use bitglobal
    use conftype
    use iomod
    use timing
    
    implicit none

    ! Irrep number
    integer(is), intent(in)  :: irrep
    
    !
    ! Sigma-vector algorithm information:
    !
    ! isigma(1)   = 1 <-> disk-based algorithm
    !               0 <-> direct algorithm
    !
    ! isigma(2:3) = Hamiltonian scratch file number and no. records
    !               in the case of a disk-based algorithm being used
    !
    integer(is), intent(in)      :: isigma(3)

    ! MRCI configuration derived type
    type(mrcfg), intent(in)      :: cfg
      
    ! Hamiltonian matrix dimension
    integer(is), intent(in)      :: matdim

    ! On-diagonal matrix elements
    real(dp), intent(in)         :: hdiag(matdim)

    ! Guess vector scratch file number
    integer(is), intent(out)     :: guessscr
    
    ! Eigenpair scratch file number and file stem
    integer(is), intent(out)     :: vecscr
    character(len=*), intent(in) :: vecstem
    
    ! Number of roots requested
    integer(is), intent(in)      :: nstates

    ! Block size
    integer(is), intent(in)      :: blocksize

    ! Maximum subspace dimension
    integer(is), intent(in)      :: maxvec

    ! Convergence tolerance
    real(dp), intent(in)         :: tol

    ! Maximum number of iterations
    integer(is), intent(in)      :: niter

    ! Preconditioner number
    integer(is), intent(in)      :: ipre
    
    ! Everything else
    integer(is)                  :: i
    
    ! Timing variables
    real(dp)                     :: tcpu_start,tcpu_end,twall_start,&
                                    twall_end
    
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    call gdavinit(matdim,blocksize,maxvec,ipre)

!----------------------------------------------------------------------
! Read the guess vectors and, optionally, the guess subspace
! eigenvectors from disk (for use in the construction of the
! generalised Davidson preconditioner)
!----------------------------------------------------------------------
    call load_guess_vectors(guessscr,matdim,blocksize,ipre)

!----------------------------------------------------------------------
! Fill in the indices of the CSFs that will *not* be coupled in the
! construction of the preconditioner
!----------------------------------------------------------------------
    if (ipre == 2) then
       idiag=1
       do i=1,subdim
          idiag(s2f(i))=0
       enddo
    endif
       
!----------------------------------------------------------------------
! Perform the generalised Davidson iterations
!----------------------------------------------------------------------
    call run_gendav(isigma,matdim,blocksize,maxvec,nstates,niter,tol,&
         hdiag,ipre)

!----------------------------------------------------------------------
! Save the eigenpairs to disk
!----------------------------------------------------------------------
    call save_eigenpairs(irrep,vecscr,vecstem,matdim,nstates)
    
!----------------------------------------------------------------------
! Finalisation
!----------------------------------------------------------------------
    call gdavfinalise
    
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'generalised_davidson')

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
    flush(6)

    return
    
  end subroutine generalised_davidson
  
!######################################################################
! gdavinit: Allocation and initialisation of the working arrays used
!           in the generalised Davidson iterations
!######################################################################
  subroutine gdavinit(matdim,blocksize,maxvec,ipre)

    use constants
    use iomod
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: matdim,blocksize,maxvec

    ! Preconditioner number
    integer(is), intent(in) :: ipre
    
!----------------------------------------------------------------------
! Make sure that the maximum subspace dimension is an integer multiple
! of the block size
!----------------------------------------------------------------------
    if (mod(maxvec,blocksize) /= 0) then
       errmsg='The maximum subspace dimension is an integer multiple'&
            //' of the block size'
       call error_control
    endif

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Subspace vectors
    !allocate(bvec(maxvec,matdim))
    allocate(bvec(matdim,maxvec))
    bvec=0.0d0
    
    ! Sigma vectors
    !allocate(sigvec(maxvec,matdim))
    allocate(sigvec(matdim,maxvec))
    sigvec=0.0d0

    ! Subspace Hamiltonian matrix
    allocate(Gmat(maxvec,maxvec))
    Gmat=0.0d0

    ! Subspace eigenvectors
    allocate(alpha(maxvec,maxvec))
    alpha=0.0d0

    ! Subspace eigenvalues
    allocate(rho(maxvec), rho1(maxvec))
    rho=0.0d0; rho1=0.0d0

    ! Residual norms
    allocate(rnorm(blocksize))
    rnorm=0.0d0

    ! Convergence flags
    allocate(iconv(blocksize))
    iconv=0

    ! Indices of the CSFs that will *not* be coupled in the
    ! construction of the preconditioner
    allocate(idiag(matdim))
    idiag=0

    ! Generalised Davidson correction vector work array
    if (ipre == 2) then
       allocate(work(matdim))
       work=0.0d0
    endif
    
    return
    
  end subroutine gdavinit

!######################################################################
! gdavfinalise: Finalisation and deallocation of working arrays
!######################################################################
  subroutine gdavfinalise

    implicit none

    deallocate(bvec)
    deallocate(sigvec)
    deallocate(Gmat)
    deallocate(alpha)
    deallocate(rho)
    deallocate(rho1)
    deallocate(rnorm)
    deallocate(iconv)
    if (allocated(idiag)) deallocate(idiag)
    if (allocated(s2f)) deallocate(s2f)
    if (allocated(subeig)) deallocate(subeig)
    if (allocated(subvec)) deallocate(subvec)
    if (allocated(work)) deallocate(work)
    
    return
    
  end subroutine gdavfinalise
  
!######################################################################
! load_guess_vectors: Reads the guess vectors from disk
!######################################################################
  subroutine load_guess_vectors(guessscr,matdim,blocksize,ipre)

    use constants
    use bitglobal
    use iomod
    
    implicit none

    ! Guess vector scratch file number
    integer(is), intent(in) :: guessscr

    ! Dimensions
    integer(is), intent(in) :: matdim,blocksize

    ! Preconditioner number
    integer(is), intent(in) :: ipre
    
    ! Everything else
    integer(is)             :: iscratch,dim,nvec,i
    
    !
    ! Open the scratch file
    !
    iscratch=scrunit(guessscr)
    open(iscratch,file=scrname(guessscr),form='unformatted',&
         status='unknown')

    !
    ! Dimensions
    !
    read(iscratch) dim
    read(iscratch) nvec

    !
    ! Sanity check
    !
    if (dim /= matdim) then
       errmsg=&
            'Inconsistent Hamiltonian matrix dimension encountered in '&
            //trim(scrname(guessscr))
       call error_control
    endif
    if (nvec /= blocksize) then
       errmsg=&
            'Inconsistent blocksize encountered in '&
            //trim(scrname(guessscr))
       call error_control
    endif

    !
    ! Guess vectors
    !
    do i=1,blocksize
       read(iscratch) bvec(:,i)
    enddo

    !
    ! Subspace eigenpairs and mappings: needed for the construction
    ! of the generalised Davidson preconditioner
    !
    if (ipre == 2) then
       read(iscratch) subdim
       allocate(s2f(subdim),subeig(subdim),subvec(subdim,subdim))
       read(iscratch) s2f
       read(iscratch) subeig
       read(iscratch) subvec
       subeig=subeig-enuc-escf
    endif
       
    !
    ! Close the scratch file
    !
    close(iscratch)

    return
    
  end subroutine load_guess_vectors
    
!######################################################################
! run_gendav: Performs the generalised Davidson iterations
!######################################################################
  subroutine run_gendav(isigma,matdim,blocksize,maxvec,nstates,niter,&
       tol,hdiag,ipre)

    use constants
    use iomod
    
    implicit none

    !
    ! Sigma-vector algorithm information:
    !
    ! isigma(1)   = 1 <-> disk-based algorithm
    !               0 <-> direct algorithm
    !
    ! isigma(2:3) = Hamiltonian scratch file number and no. records
    !               in the case of a disk-based algorithm being used
    !
    integer(is), intent(in) :: isigma(3)
    
    ! Dimensions
    integer(is), intent(in) :: matdim,blocksize,maxvec,nstates

    ! Davidson parameters
    integer(is), intent(in) :: niter
    real(dp), intent(in)    :: tol

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: hdiag(matdim)

    ! Preconditioner number
    integer(is), intent(in) :: ipre
    
    ! Everything else
    integer(is)             :: k,i

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
! currdim: the current dimension of the subspace
! nnew:    the no. new subspace vectors added in a given iteration
! nconv:   the no. of converged roots
! nsigma:  the total no. sigma vectors calculated
!----------------------------------------------------------------------
    currdim=blocksize
    nnew=blocksize
    nconv=0
    nsigma=0
    
!----------------------------------------------------------------------
! Perform the generalised Davidson iterations
!----------------------------------------------------------------------
    ! Loop over iterations
    do k=1,niter

       ! Calculate the sigma vectors
       call sigma_vectors(matdim,blocksize,hdiag,isigma)

       ! Compute the new elements in the subspace Hamiltonian
       ! matrix
       call subspace_hamiltonian(matdim,blocksize)

       ! Compute the eigenpairs of the subspace Hamiltonian
       call subspace_diag(matdim,blocksize)
       
       ! Compute the residual vectors. Note that these will be stored
       ! in the bvec array and subsequently transformed in place
       ! to obtain the correction vectors, and then the new
       ! subspace vectors
       call residual_vectors(matdim,blocksize,tol)

       ! Print the report for this iteration
       call print_report(k,nstates)
       
       ! Stop here if all the roots are converged
       if (converged(nstates,tol)) exit
       
       ! Compute the correction vectors
       call correction_vectors(matdim,blocksize,hdiag,ipre)

       ! Compute the new subspace vectors
       call subspace_vectors(matdim,blocksize,tol)
       
       ! Subspace collapse?
       if (currdim+nnew*2 > maxvec) then
          call subspace_collapse(matdim,blocksize)
       else
          currdim=currdim+nnew
       endif

       ! If we are here and this is the last iteration, then
       ! we failed to converge all roots
       if (k == niter) then
          errmsg='Not all roots converged'
          call error_control
       endif
       
    enddo
    
    ! If we are here then convergence has been achieved
    write(6,'(/,x,a)') 'All roots converged'
    write(6,'(/,x,a,x,i0)') 'N_sigma:',nsigma
    
    return
    
  end subroutine run_gendav

!######################################################################
! sigma_vectors: Calculation of a batch of sigma vectors
!######################################################################
  subroutine sigma_vectors(matdim,blocksize,hdiag,isigma)

    use constants
    use sigma_mrci
    use iomod
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: matdim,blocksize

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: hdiag(matdim)

    !
    ! Sigma-vector algorithm information:
    !
    ! isigma(1)   = 1 <-> disk-based algorithm
    !               0 <-> direct algorithm
    !
    ! isigma(2:3) = Hamiltonian scratch file number and no. records
    !               in the case of a disk-based algorithm being used
    !
    integer(is), intent(in) :: isigma(3)

    ! Everything else
    integer(is)             :: ki,kf

!----------------------------------------------------------------------
! Update the total no. sigma vector calculations
!----------------------------------------------------------------------
    nsigma=nsigma+nnew
       
!----------------------------------------------------------------------
! Indices of the first and last subspace vectors for which
! sigma-vectors are required
!----------------------------------------------------------------------
    ki=currdim-nnew+1
    kf=currdim

!----------------------------------------------------------------------
! Compute the sigma vectors
!----------------------------------------------------------------------
    select case(isigma(1))
    
    case(0)
       ! Direct algorithm
       errmsg='Direct sigma-vector builds not yet implemented'
       call error_control
    case(1)
       ! Disk-based algorithm
       call sigma_disk(isigma(2),isigma(3),matdim,nnew,&
            bvec(:,ki:kf),sigvec(:,ki:kf),hdiag)
    end select
    
    return
    
  end subroutine sigma_vectors

!######################################################################
! subspace_hamiltonian: Calculation of the elements of the subspace
!                       Hamiltonian
!######################################################################
  subroutine subspace_hamiltonian(matdim,blocksize)

    use constants
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: matdim,blocksize

    ! Everything else
    integer(is)             :: i,j,j1

    ! Loop over subspace vectors
    do i=1,currdim

       ! Loop over new sigma vectors
       do j=1,nnew
       
          ! Index of the subspace and sigma vectors
          j1=currdim-nnew+j
          
          ! G_ij = < b_i | H | b_j > = < b_i | sigma_j >
          Gmat(i,j1)=dot_product(bvec(:,i),sigvec(:,j1))
          Gmat(j1,i)=Gmat(i,j1)
          
       enddo

    enddo
    
    return
    
  end subroutine subspace_hamiltonian

!######################################################################
! subspace_diag: Calculation of the eigenpairs of the subspace
!                Hamiltonian matrix
!######################################################################
  subroutine subspace_diag(matdim,blocksize)

    use constants
    use utils
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: matdim,blocksize

    !
    ! Diagonalise the subspace Hamiltonian
    !
    call diag_matrix_real(Gmat(1:currdim,1:currdim),rho(1:currdim),&
         alpha(1:currdim,1:currdim),currdim)
    
    return
    
  end subroutine subspace_diag
    
!######################################################################
! residual_vectors: Calculation of the residual vectors
!######################################################################
  subroutine residual_vectors(matdim,blocksize,tol)

    use constants
    use iomod
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: matdim,blocksize

    ! Residual norm convergence threshold
    real(dp), intent(in)    :: tol

    ! Working arrays
    real(dp), allocatable   :: rvec(:)
        
    ! Everything else
    integer(is)             :: ki,kf,k,k1,i,inc
    
    !******************************************************************
    ! r_k = Sum_i alpha_ik * (sigma_k - rho_k b_i),
    !
    ! where alpha_k and rho_k are the subspace eigenpairs, and sigma
    ! and b are the sigma and subspace vectors.    
    !******************************************************************

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(rvec(matdim))
    rvec=0
    
!----------------------------------------------------------------------
! Compute the residual vectors
! r_k = Sum_i alpha_ik * (sigma_k - rho_k b_i),
!----------------------------------------------------------------------
    ! Initialisation
    ki=currdim+1
    kf=currdim+nnew
    bvec(:,ki:kf)=0.0d0

    ! Loop over roots
    nnew=0
    do k=1,blocksize

       ! Compute the residual vector
       rvec=0.0d0
       do i=1,currdim
          rvec=rvec+alpha(i,k)*(sigvec(:,i)-rho(k)*bvec(:,i))          
       enddo

       ! Residual norm
       rnorm(k)=sqrt(dot_product(rvec,rvec))

       ! Update the convergence information
       if (rnorm(k) < tol) iconv(k)=1

       ! Save the residual vector and corresponding eigenvalue if it
       ! corresponds to an unconverged root
       if (iconv(k) == 0) then
          nnew=nnew+1
          bvec(:,ki-1+nnew)=rvec
          rho1(nnew)=rho(k)
       endif
       
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(rvec)
    
    return
    
  end subroutine residual_vectors
  
!######################################################################
! correction_vectors: Calculation of the correction vectors via the
!                     operation of the preconditioner on the residual
!                     vectors
!######################################################################
  subroutine correction_vectors(matdim,blocksize,hdiag,ipre)

    use constants

    implicit none

    ! Dimensions
    integer(is), intent(in)  :: matdim,blocksize

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: hdiag(matdim)

    ! Preconditioner number
    integer(is), intent(in)  :: ipre
    
    ! Working arrays
    real(dp), allocatable    :: subpre(:,:),subvec1(:,:),subvecT(:,:)
    
    ! Everything else
    integer(is)              :: ki,kf,k,k1,i,i1,j,j1,l
    
!----------------------------------------------------------------------
! Indices of the positions in the bvec array in which the correction
! vectors will be stored
!----------------------------------------------------------------------
    ki=currdim+1
    kf=currdim+nnew
    
!----------------------------------------------------------------------
! Diagonal preconditioned residue correction vectors
!----------------------------------------------------------------------
    if (ipre == 1) then
    
       ! Loop over correction vectors
       k1=0
       do k=ki,kf
          k1=k1+1
          
          ! Loop over elements of the correction vector
          do i=1,matdim
             bvec(i,k)=-bvec(i,k)/(hdiag(i)-rho1(k1))
          enddo

       enddo

    endif
       
!----------------------------------------------------------------------
! Generalised Davidson correction vectors
!----------------------------------------------------------------------
    if (ipre == 2) then

       ! Allocate arrays
       allocate(subpre(subdim,subdim), subvec1(subdim,subdim), &
            subvecT(subdim,subdim))
       subpre=0.0d0; subvec1=0.0d0
       
       ! Transpose of the coupled CSF subspace eigenvector matrix
       subvecT=transpose(subvec)
       
       ! Loop over correction vectors
       k1=0
       do k=ki,kf
          k1=k1+1
          
          ! Residual vector
          work=bvec(:,k)
          
          ! Initialise the correction vector to zero
          bvec(:,k)=0.0d0
          
          ! Uncoupled CSF part of the correction vector
          do i=1,matdim
             if (idiag(i) == 1) &
                  bvec(i,k)=-work(i)/(hdiag(i)-rho1(k1))
          enddo
          
          ! Coupled CSF preconditioner
          do i=1,subdim
             subvec1(:,i)=-subvec(:,i)/(subeig(i)-rho1(k1))
          enddo
          call dgemm('N','N',subdim,subdim,subdim,1.0d0,subvec1,&
               subdim,subvecT,subdim,0.0d0,subpre,subdim)
          
          ! Coupled CSF part of the correction vector
          do i=1,subdim
             i1=s2f(i)
             do j=1,subdim
                j1=s2f(j)
                bvec(i1,k)=bvec(i1,k) + subpre(i,j) * work(j1)
             enddo
          enddo
          
       enddo

       ! Deallocate arrays
       deallocate(subpre, subvec1, subvecT)

    endif
       
    return
    
  end subroutine correction_vectors

!######################################################################
! subspace_vectors: Calculation of the new subspace vectors via the
!                   orthogonalisation of the correction vectors
!                   against the previous subspace vectors
!######################################################################
  subroutine subspace_vectors(matdim,blocksize,tol)

    use constants
    
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: matdim,blocksize

    ! Convergence threshold
    real(dp), intent(in)     :: tol

    ! Orthogonalised correction vector norms
    real(dp)                 :: bnorm(blocksize)
      
    ! Everything else
    integer(is)              :: ki,kf,k,k1,i,n

!----------------------------------------------------------------------
! Indices of the positions in the bvec array in which the new subspace
! vectors will be stored
!----------------------------------------------------------------------
    ki=currdim+1
    kf=currdim+nnew
    
!----------------------------------------------------------------------
! Compute the (non-normalised) new subspace vectors
!----------------------------------------------------------------------
    ! Loop over correction vectors
    k1=0
    do k=ki,kf
       k1=k1+1

       ! Orthogonalise the correction vector against all previous
       ! subspace vectors
       do i=1,k-1
          bvec(:,k)=bvec(:,k) &
               -dot_product(bvec(:,k),bvec(:,i))*bvec(:,i) &
               /dot_product(bvec(:,i),bvec(:,i))
       enddo

       ! Norm of the orthogonalised correction vector
       bnorm(k1)=sqrt(dot_product(bvec(:,k),bvec(:,k)))
       
    enddo
    
!----------------------------------------------------------------------
! Expand the subspace by adding orthonormalised correction vectors 
!----------------------------------------------------------------------
    ! Loop over orthogonalised correction vectors
    n=0
    do k=ki,kf
       
       ! Add the orthonormalised correction vector to the subspace
       n=n+1
       bvec(:,k)=bvec(:,k)/bnorm(n)

    enddo
    
    return
    
  end subroutine subspace_vectors

!######################################################################
! subspace_collapse: Collapses the subspace
!######################################################################
  subroutine subspace_collapse(matdim,blocksize)

    use constants
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: matdim,blocksize

    ! Everything else
    integer(is)             :: i,k

!----------------------------------------------------------------------
! Compute the Ritz vectors. We will use the sigvec array as a working
! array here to store the Ritz vectors
!----------------------------------------------------------------------
    ! Initialisation
    sigvec=0.0d0

    ! Loop over Ritz vectors
    do k=1,blocksize

       ! Compute the kth Ritz vector 
       do i=1,currdim
          sigvec(:,k)=sigvec(:,k)+alpha(i,k)*bvec(:,i)
       enddo
       
    enddo

!----------------------------------------------------------------------
! Collapse the subspace to be spanned by the lowest-lying Ritz vectors
!----------------------------------------------------------------------
    ! New subspace dimension
    currdim=blocksize
    nnew=blocksize

    ! Save the Ritz vectors as the new subspace vectors
    do k=1,blocksize
       bvec(:,k)=sigvec(:,k)
    enddo

    return
    
  end subroutine subspace_collapse

!######################################################################
! converged: Returns .true. if all the roots are converged, else
!            returns .false.
!######################################################################
  function converged(nstates,tol)

    use constants

    implicit none

    logical :: converged

    ! Number of roots
    integer(is), intent(in) :: nstates

    ! Convergence threshold
    real(dp), intent(in)    :: tol
    
    nconv=sum(iconv(1:nstates))

    if (nconv >= nstates) then
       converged=.true.
    else
       converged=.false.
    endif
    
    return
    
  end function converged
  
!######################################################################
! print_report: Outputs the energies, residual norms, etc. for the
!               current iteration
!######################################################################
  subroutine print_report(k,nstates)

    use constants
    use bitglobal
    
    implicit none

    ! Iteration number
    integer(is), intent(in) :: k

    ! Dimensions
    integer(is), intent(in) :: nstates
    
    ! Everything else
    integer(is)             :: i
    
!----------------------------------------------------------------------
! Table header
!----------------------------------------------------------------------
    if (k == 1) then
       write(6,'(/,43a)') ('*',i=1,43)
       write(6,'(x,a,2x,a,3x,a,7x,a)') &
            'Iteration','Nvec','Max rnorm','Nconv'
       write(6,'(43a)') ('*',i=1,43)
    endif
    
!----------------------------------------------------------------------
! Information for the current iteration
!----------------------------------------------------------------------
    write(6,'(x,i4,7x,i4,3x,E13.7,3x,i4)') &
         k,currdim,maxval(rnorm(1:nstates)),sum(iconv(1:nstates))

!----------------------------------------------------------------------    
! Table footer
!----------------------------------------------------------------------
    if (sum(iconv(1:nstates)) == nstates) &
         write(6,'(43a)') ('*',i=1,43)
    
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
    flush(6)
    
    return
    
  end subroutine print_report
    
!######################################################################
! save_eigenpairs: Saves the converged eigenpairs to disk
!######################################################################
  subroutine save_eigenpairs(irrep,vecscr,vecstem,matdim,nstates)

    use constants
    use bitglobal
    use iomod
    
    implicit none

    ! Irrep number
    integer(is), intent(in)      :: irrep
    
    ! Eigenpair scratch file number and file stem
    integer(is), intent(out)     :: vecscr
    character(len=*), intent(in) :: vecstem
    
    ! Dimensions
    integer(is), intent(in)      :: matdim,nstates

    ! Everything else
    integer(is)                  :: i,k
    integer(is)                  :: iscratch
    character(len=60)            :: vecfile
    character(len=2)             :: amult,airrep
    
!----------------------------------------------------------------------
! Compute the Ritz vectors. We will use the sigvec array as a working
! array here to store the Ritz vectors
!----------------------------------------------------------------------
    ! Initialisation
    sigvec=0.0d0

    ! Loop over Ritz vectors
    do k=1,nstates

       ! Compute the kth Ritz vector 
       do i=1,currdim
          sigvec(:,k)=sigvec(:,k)+alpha(i,k)*bvec(:,i)
       enddo
       
    enddo

!----------------------------------------------------------------------
! Add on E_SCF and E_nuc to the eigenvalues
!----------------------------------------------------------------------
    rho=rho+escf+enuc
    
!----------------------------------------------------------------------
! Register the scratch file
!----------------------------------------------------------------------
    write(amult,'(i0)') imult
    write(airrep,'(i0)') irrep
    call scratch_name(trim(vecstem)//'.mult'//trim(amult)//&
         '.sym'//trim(airrep),vecfile)
    call register_scratch_file(vecscr,vecfile)

!----------------------------------------------------------------------
! Save the eigenpairs to disk
!----------------------------------------------------------------------
    ! Open the scratch file
    iscratch=scrunit(vecscr)
    open(iscratch,file=scrname(vecscr),form='unformatted',&
         status='unknown')
    
    ! No. CSFs
    write(iscratch) matdim
    
    ! No. roots
    write(iscratch) nstates
    
    ! Eigenvalues
    write(iscratch) rho(1:nstates)
    
    ! Eigenvectors
    do i=1,nstates
       write(iscratch) sigvec(:,i)
    enddo
     
    ! Close the scratch file
    close(iscratch)
    
    return
    
  end subroutine save_eigenpairs
  
!######################################################################
  
end module gendav
