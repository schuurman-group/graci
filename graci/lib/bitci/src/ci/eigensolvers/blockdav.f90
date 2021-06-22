!**********************************************************************
! Block Davidson diagonalisation module
!**********************************************************************
module blockdav

contains

!######################################################################
! block_davidson: Master routine for a block Davidson calculation
!######################################################################
  subroutine block_davidson(irrep,isigma,cfg,matdim,hdiag,guessscr,&
       vecscr,vecstem,nstates,blocksize,maxvec,ldeflate,tol,niter)

    use constants
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

    ! Subspace deflation
    logical, intent(in)          :: ldeflate

    ! Convergence tolerance
    real(dp), intent(in)         :: tol

    ! Maximum number of iterations
    integer(is), intent(in)      :: niter
    
    ! Working arrays
    real(dp), allocatable        :: vmat(:,:),wmat(:,:),rmat(:,:)
    real(dp), allocatable        :: ritzvec(:,:),res(:,:),reigvec(:,:)
    real(dp), allocatable        :: reigval(:),norm(:)
    
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
    call davinitialise(matdim,blocksize,maxvec,vmat,wmat,rmat,&
         ritzvec,res,reigvec,reigval,norm)

!----------------------------------------------------------------------
! Read the guess vectors from disk
!----------------------------------------------------------------------
    call read_guessvec(guessscr,matdim,maxvec,blocksize,vmat)

!----------------------------------------------------------------------
! Perform the block Davidson iterations
!----------------------------------------------------------------------
    call run_block_davidson(irrep,isigma,matdim,hdiag,maxvec,&
         blocksize,nstates,ldeflate,tol,niter,vmat,wmat,rmat,ritzvec,&
         res,reigvec,reigval,norm,vecscr,vecstem)
    
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'block_davidson')

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
    flush(6)
    
    return
    
  end subroutine block_davidson
  
!######################################################################
! davinit: initialisation of variables and arrays used in the block
!          Davidson iterations
!######################################################################
  subroutine davinitialise(matdim,blocksize,maxvec,vmat,wmat,rmat,&
       ritzvec,res,reigvec,reigval,norm)

    use constants
        
    implicit none

    ! Hamiltonian matrix dimension
    integer(is), intent(in) :: matdim

    ! Block size
    integer(is), intent(in) :: blocksize

    ! Maximum subspace dimension
    integer(is), intent(in) :: maxvec
    
    ! Working arrays
    real(dp), allocatable   :: vmat(:,:),wmat(:,:),rmat(:,:)
    real(dp), allocatable   :: ritzvec(:,:),res(:,:),reigvec(:,:)
    real(dp), allocatable   :: reigval(:),norm(:)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Matrix of subspace vectors
    allocate(vmat(matdim,maxvec))
    vmat=0.0d0
      
    ! Matrix-vector product
    allocate(wmat(matdim,maxvec))
    wmat=0.0d0
    
    ! Rayleigh matrix
    allocate(rmat(maxvec,maxvec))
    rmat=0.0d0
    
    ! n=blocksize lowest eigenpairs of the Rayleigh matrix
    allocate(reigvec(maxvec,blocksize))
    allocate(reigval(blocksize))
    reigvec=0.0d0
    reigval=0.0d0
    
    ! Ritz vectors
    allocate(ritzvec(matdim,blocksize))
    ritzvec=0.0d0
    
    ! Residual vectors
    allocate(res(matdim,blocksize))
    res=0.0d0
    
    ! Norms of the residual vectors
    allocate(norm(blocksize))
    norm=0.0d0
    
    return
    
  end subroutine davinitialise

!######################################################################
! read_guessvec: reads the guess vectors from disk
!######################################################################
  subroutine read_guessvec(guessscr,matdim,maxvec,blocksize,vmat)

    use constants
    use bitglobal
    use iomod
    
    implicit none

    ! Guess vector scratch file number
    integer(is), intent(in) :: guessscr

    ! Hamiltonian matrix dimension
    integer(is), intent(in) :: matdim

    ! Maximum subspace dimension
    integer(is), intent(in) :: maxvec
    
    ! Block size
    integer(is), intent(in) :: blocksize

    ! Subspace vectors
    real(dp), intent(out)   :: vmat(matdim,maxvec)
    
    ! Everything else
    integer(is)             :: iscratch,dim,nvec,i,j

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
       read(iscratch) vmat(:,i)
    enddo

    !
    ! Close the scratch file
    !
    close(iscratch)
    
    return
    
  end subroutine read_guessvec
    
!######################################################################
! run_block_davidson: Runs the block Davidson iterations
!######################################################################
  subroutine run_block_davidson(irrep,isigma,matdim,hdiag,maxvec,&
       blocksize,nstates,ldeflate,tol,niter,vmat,wmat,rmat,ritzvec,&
       res,reigvec,reigval,norm,vecscr,vecstem)

    use constants
    use iomod
    
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
    
    ! Hamiltonian matrix dimension
    integer(is), intent(in)      :: matdim

    ! On-diagonal matrix elements
    real(dp), intent(in)         :: hdiag(matdim)
    
    ! Maximum subspace dimension
    integer(is), intent(in)      :: maxvec
    
    ! Block size
    integer(is), intent(in)      :: blocksize

    ! Number of roots requested
    integer(is), intent(in)      :: nstates

    ! Subspace deflation
    logical, intent(in)          :: ldeflate

    ! Convergence tolerance
    real(dp), intent(in)         :: tol

    ! Maximum number of iterations
    integer(is), intent(in)      :: niter

    ! Working arrays
    real(dp), intent(inout)      :: vmat(matdim,maxvec)
    real(dp), intent(inout)      :: wmat(matdim,maxvec)
    real(dp), intent(inout)      :: rmat(maxvec,maxvec)
    real(dp), intent(inout)      :: ritzvec(matdim,blocksize)
    real(dp), intent(inout)      :: res(matdim,blocksize)
    real(dp), intent(inout)      :: reigvec(maxvec,blocksize)
    real(dp), intent(inout)      :: reigval(blocksize)
    real(dp), intent(inout)      :: norm(blocksize)

    ! Eigenpair scratch file number and file stem
    integer(is), intent(out)     :: vecscr
    character(len=*), intent(in) :: vecstem
    
    ! Preconditioner (hard-coded for now)
    integer(is), parameter       :: ipre=1
    
    ! Everything else
    integer(is)                  :: currdim,blocksize_curr,maxvec_curr
    integer(is)                  :: nstates_curr,nconv,nconv_prev
    integer(is)                  :: k
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
! currdim:        the current dimension of the subspace
!
! blocksize_curr: the current blocksize
!
! maxvec_curr:    the current maximum subspace dimension
!
! nstates_curr:   the current no. of states that we are solving for,
!                 i.e., the current no. of unconverged roots
!
! nconv:          the no. of converged roots
!----------------------------------------------------------------------
! Note that the blocksize, maximum subspace dimension and nstates will
! only change throughout the Davidson iterations if the converged
! vectors are removed from from the subspace, i.e., if ldeflate=.true.
!----------------------------------------------------------------------
    currdim=blocksize
    blocksize_curr=blocksize
    maxvec_curr=maxvec
    nstates_curr=nstates
    nconv=0
    nconv_prev=0
    
!----------------------------------------------------------------------
! Perform the Davidson iterations
!----------------------------------------------------------------------
    ! Loop over iterations
    do k=1,niter

       ! Calculate the matrix-vector products
       call hxvec(isigma,matdim,maxvec,currdim,wmat,vmat,hdiag)

       ! Calculate the Rayleigh matrix
       call calcrmat(matdim,maxvec,currdim,wmat,vmat,rmat)

       ! Diagonalise the Rayleigh matrix
       call diagrmat(maxvec,blocksize,blocksize_curr,currdim,rmat,&
            reigvec,reigval)

       ! Calculate the Ritz vectors
       call calcritzvec(matdim,maxvec,blocksize,blocksize_curr,&
            currdim,vmat,reigvec,ritzvec)

       ! Calculate the residuals
       call calcres(matdim,maxvec,blocksize,blocksize_curr,currdim,&
            nstates_curr,wmat,reigvec,reigval,ritzvec,res,norm,&
            ldeflate,nconv,tol)

       ! Output progress
       call wrtable(k,nstates_curr,blocksize,norm,reigval,tol)

       ! Exit if we have converged all roots
       if (nconv.eq.nstates) then
          write(6,'(/,2x,a,/)') 'All roots converged'
          exit
       endif

       ! Expand the subspace
       call subspace_expansion(matdim,maxvec,blocksize,&
            currdim,maxvec_curr,blocksize_curr,nconv,nconv_prev,&
            nstates,nstates_curr,vmat,ritzvec,reigval,res,norm,hdiag,&
            ipre,ldeflate,tol)

       ! Keep track of the no. of roots converged so far
       nconv_prev=nconv
       
    enddo

!----------------------------------------------------------------------
! Die here if we haven't converged all eigenpairs
!----------------------------------------------------------------------
    if (nconv.ne.nstates) then
       errmsg='Not all vectors have converged...'
       call error_control
    endif

!----------------------------------------------------------------------
! Write the eigenpairs to disk
!----------------------------------------------------------------------
    call wreigenpairs(irrep,vecscr,vecstem,matdim,blocksize,reigval,&
         ritzvec,nstates,nstates_curr,nconv_prev,ldeflate)

    return
    
  end subroutine run_block_davidson
  
!######################################################################
! hxvec: calculation of a batch of matrix-vector products
!######################################################################
  subroutine hxvec(isigma,matdim,maxvec,currdim,wmat,vmat,hdiag)
    
    use constants
    use sigma_mrci
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
    
    ! Hamiltonian matrix dimension
    integer(is), intent(in) :: matdim

    ! Maximum subspace dimension
    integer(is), intent(in) :: maxvec
    
    ! Current subspace dimension
    integer(is), intent(in) :: currdim

    ! Matrix of subspace vectors
    real(dp), intent(in)    :: vmat(matdim,maxvec)
      
    ! Matrix-vector products
    real(dp), intent(out)   :: wmat(matdim,maxvec)

    ! On-diagonal matrix elements
    real(dp), intent(in)    :: hdiag(matdim)

    select case(isigma(1))

    case(0) ! Direct algorithm
       errmsg='Direct sigma-vector builds not yet implemented'
       call error_control
       
    case(1) ! Disk-based algorithm

       call sigma_disk_old(isigma(2),isigma(3),matdim,maxvec,currdim,&
            wmat,vmat,hdiag)
       
    end select
       
    return
    
  end subroutine hxvec

!######################################################################
! calcrmat: calculation of the Rayleigh matrix, i.e., the projection
!           of the Hamiltonian onto the space spanned by the subspace
!           vectors
!######################################################################
  subroutine calcrmat(matdim,maxvec,currdim,wmat,vmat,rmat)

    use constants
    
    implicit none

    ! Hamiltonian matrix dimension
    integer(is), intent(in) :: matdim

    ! Maximum subspace dimension
    integer(is), intent(in) :: maxvec
    
    ! Current subspace dimension
    integer(is), intent(in) :: currdim

    ! Matrix of subspace vectors
    real(dp), intent(in)    :: vmat(matdim,maxvec)
      
    ! Matrix-vector products
    real(dp), intent(out)   :: wmat(matdim,maxvec)

    ! Rayleigh matrix
    real(dp), intent(inout) :: rmat(maxvec,maxvec)
    
    !
    ! R = V^T H V
    !
    call dgemm('T','N',currdim,currdim,matdim,1.0d0,&
         vmat(:,1:currdim),matdim,wmat(:,1:currdim),matdim,0.0d0,&
         rmat(1:currdim,1:currdim),currdim)

    return
    
  end subroutine calcrmat
  
!######################################################################
! diagrmat: diagonalisation of the Rayleugh matrix
!######################################################################
  subroutine diagrmat(maxvec,blocksize,blocksize_curr,currdim,rmat,&
       reigvec,reigval)

    use constants
    use iomod

    use utils
    
    implicit none

    ! Maximum subspace dimension
    integer(is), intent(in) :: maxvec

    ! Block size
    integer(is), intent(in) :: blocksize,blocksize_curr
    
    ! Current subspace dimension
    integer(is), intent(in) :: currdim

    ! Rayleigh matrix and it's eigenpairs
    real(dp), intent(inout) :: rmat(maxvec,maxvec)
    real(dp), intent(out)   :: reigvec(maxvec,blocksize)
    real(dp), intent(out)   :: reigval(blocksize)

    ! Everything else
    integer(is)             :: error,e2,i
    real(dp), allocatable   :: work(:)
    real(dp), allocatable   :: val(:)
    real(dp), allocatable   :: vec(:,:)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(work(3*currdim))
    allocate(val(currdim))
    allocate(vec(currdim,currdim))

!----------------------------------------------------------------------
! Diagonalise the Rayleigh matrix
!----------------------------------------------------------------------
    error=0
    e2=3*currdim
    vec=rmat(1:currdim,1:currdim)
    
    call dsyev('V','U',currdim,vec,currdim,val,work,e2,error)
    
    if (error.ne.0) then
       errmsg='Diagonalisation of the Rayleigh matrix in &
            subroutine diagrmat failed'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Save the n=blocksize_curr lowest eigenpairs to be used in the
! calculation of the Ritz vectors and residuals
!----------------------------------------------------------------------
    reigvec(1:currdim,1:blocksize_curr)=vec(1:currdim,1:blocksize_curr)
    reigval(1:blocksize_curr)=val(1:blocksize_curr)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(work)
    deallocate(val)
    deallocate(vec)
    
    return
    
  end subroutine diagrmat

!######################################################################
! calcritzvec: calculation of the Ritz vectors
!######################################################################
  subroutine calcritzvec(matdim,maxvec,blocksize,blocksize_curr,&
       currdim,vmat,reigvec,ritzvec)

    use constants
    
    implicit none

    ! Hamiltonian matrix dimension
    integer(is), intent(in) :: matdim
    
    ! Maximum subspace dimension
    integer(is), intent(in) :: maxvec

    ! Block size
    integer(is), intent(in) :: blocksize,blocksize_curr
    
    ! Current subspace dimension
    integer(is), intent(in) :: currdim

    ! Subspace vectors
    real(dp), intent(in)    :: vmat(matdim,maxvec)
    
    ! Rayleigh matrix eigenvectors
    real(dp), intent(in)    :: reigvec(maxvec,blocksize)

    ! Ritz vectors
    real(dp), intent(out)   :: ritzvec(matdim,blocksize)

    call dgemm('N','N',matdim,blocksize_curr,currdim,1.0d0,&
         vmat(1:matdim,1:currdim),matdim,&
         reigvec(1:currdim,1:blocksize_curr),currdim,0.0d0,&
         ritzvec(1:matdim,1:blocksize_curr),matdim)
    
    return
    
  end subroutine calcritzvec

!######################################################################
! calcres: calculation of the residual vectors
!######################################################################
  subroutine calcres(matdim,maxvec,blocksize,blocksize_curr,currdim,&
       nstates_curr,wmat,reigvec,reigval,ritzvec,res,norm,ldeflate,&
       nconv,tol)

    use constants
    
    implicit none

    ! Hamiltonian matrix dimension
    integer(is), intent(in)   :: matdim
    
    ! Maximum subspace dimension
    integer(is), intent(in)   :: maxvec

    ! Block size
    integer(is), intent(in)   :: blocksize,blocksize_curr
    
    ! Current subspace dimension
    integer(is), intent(in)   :: currdim

    ! Current no. roots being worked on
    integer(is), intent(in)   :: nstates_curr
    
    ! Subspace vectors
    real(dp), intent(in)      :: wmat(matdim,maxvec)
    
    ! Rayleigh matrix eigenpairs
    real(dp), intent(in)      :: reigvec(maxvec,blocksize)
    real(dp), intent(in)      :: reigval(blocksize)
    
    ! Ritz vectors
    real(dp), intent(in)      :: ritzvec(matdim,blocksize)

    ! Residual vectors and their norms
    real(dp), intent(out)     :: res(matdim,blocksize)
    real(dp), intent(out)     :: norm(blocksize)

    ! Deflation flag
    logical, intent(in)       :: ldeflate

    ! No. converged roots
    integer(is), intent(inout) :: nconv

    ! Convergence threshold
    real(dp), intent(in)       :: tol
    
    ! Everything else
    integer(is)                :: i

!----------------------------------------------------------------------
! Residual vectors: r_i = lambda_i * x_i - W * y_i
!----------------------------------------------------------------------
! r_i       ith residual vector
!
! lambda_i  ith eigenvalue of the Rayleigh matrix
!
! x_i       ith Ritz vector
!
! W         = H * V (Hamiltonian multiplied against the matrix of
!                   subspace vectors)
!
! y_i      ith eigenvector of the Rayleigh matrix
!----------------------------------------------------------------------
    ! -W * y_i
    call dgemm('N','N',matdim,blocksize_curr,currdim,-1.0d0,&
         wmat(1:matdim,1:currdim),matdim,&
         reigvec(1:currdim,1:blocksize_curr),currdim,0.0d0,&
         res(1:matdim,1:blocksize_curr),matdim)
      
    ! lambda_i * x_i -W * y_i
    do i=1,blocksize_curr
       res(:,i)=res(:,i)+reigval(i)*ritzvec(:,i)
    enddo

!----------------------------------------------------------------------
! Norms of the residual vectors
!----------------------------------------------------------------------
    do i=1,blocksize_curr
       norm(i)=dot_product(res(:,i),res(:,i))
       norm(i)=sqrt(norm(i))
    enddo

!----------------------------------------------------------------------
! Keep track of the no. converged roots
!----------------------------------------------------------------------
    if (.not.ldeflate) nconv=0

    do i=1,nstates_curr
       if (norm(i).lt.tol) nconv=nconv+1
    enddo
    
    return
    
  end subroutine calcres

!######################################################################
! wrtable: outputs the progress of the Davidson iterations
!######################################################################
  subroutine wrtable(k,nstates_curr,blocksize,norm,reigval,tol)

    use constants
    use bitglobal
    
    implicit none

    ! Iteration number
    integer(is), intent(in) :: k

    ! Dimensions
    integer(is), intent(in) :: nstates_curr,blocksize

    ! Residual norms
    real(dp), intent(in)    :: norm(blocksize)

    ! Approximate eigenvalues
    real(dp), intent(in)    :: reigval(blocksize)

    ! Convergence threshold
    real(dp), intent(in)    :: tol
    
    ! Everything else
    integer(is)             :: i,j
    character(len=1)        :: aconv

!----------------------------------------------------------------------
! Table header
!----------------------------------------------------------------------
    if (k.eq.1) then
       write(6,'(53a)') ('*',j=1,53)
       write(6,'(4(a,6x))') &
            'Iteration','Energies','Residuals','Converged'
       write(6,'(53a)') ('*',j=1,53)
    endif

!----------------------------------------------------------------------
! Information from the current iteration
!----------------------------------------------------------------------
    write(6,*)
    do i=1,nstates_curr
       if (norm(i).lt.tol) then
          aconv='y'
       else
          aconv='n'
       endif
       
       if (i.eq.1) then
          write(6,'(i4,10x,F12.7,3x,E13.7,2x,a1)') &
               k,reigval(i)+escf,norm(i),aconv
       else
          write(6,'(14x,F12.7,3x,E13.7,2x,a1)') &
               reigval(i)+escf,norm(i),aconv
       endif
    enddo

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
    flush(6)
    
    return
    
  end subroutine wrtable

!######################################################################
! subspace_expansion: computes the new subspace vectors
!######################################################################
  subroutine subspace_expansion(matdim,maxvec,blocksize,currdim,&
       maxvec_curr,blocksize_curr,nconv,nconv_prev,nstates,&
       nstates_curr,vmat,ritzvec,reigval,res,norm,hdiag,ipre,&
       ldeflate,tol)

    use constants
    use iomod

    implicit none

    ! Dimensions
    integer(is), intent(in)    :: matdim,maxvec,blocksize,&
                                  nconv,nconv_prev,nstates
    integer(is), intent(inout) :: blocksize_curr,maxvec_curr,&
                                  nstates_curr,currdim
    
    ! Subspace vectors
    real(dp), intent(inout)    :: vmat(matdim,maxvec)

    ! Ritz vectors
    real(dp), intent(inout)    :: ritzvec(matdim,blocksize)

    ! Residuals
    real(dp), intent(inout)    :: res(matdim,blocksize)

    ! Eigenvalues of the Rayleigh matrix
    real(dp), intent(inout)    :: reigval(blocksize)

    ! Residual norms
    real(dp), intent(inout)    :: norm(blocksize)
    
    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)       :: hdiag(matdim)

    ! Preconditoner and deflation flags
    integer(is), intent(in)    :: ipre
    logical, intent(in)        :: ldeflate

    ! Convergence threshold
    real(dp), intent(in)       :: tol

    ! Everything else
    logical                    :: lcollapse
    
!----------------------------------------------------------------------
! Removal of converged vectors from the subspace
!----------------------------------------------------------------------
    if (ldeflate .and. nconv > nconv_prev) &
         call deflate(matdim,maxvec,blocksize,ritzvec,reigval,res,&
         norm,nstates,nstates_curr,currdim,maxvec_curr,&
         blocksize_curr,nconv,nconv_prev,tol)

!----------------------------------------------------------------------    
! Determine whether or not we need to collapse the subspace
!----------------------------------------------------------------------
    if (currdim.le.maxvec_curr-blocksize_curr) then
       lcollapse=.false.
    else
       lcollapse=.true.
    endif

    ! Force a collapse of the subspace if a deflation has occurred
    if (ldeflate.and.nconv.gt.nconv_prev) lcollapse=.true.

!----------------------------------------------------------------------
! Calculate the new subspace vectors
!----------------------------------------------------------------------
    if (ipre.eq.1) then
       ! Diagonal preconditioned residue
       call dpr(matdim,maxvec,blocksize,blocksize_curr,currdim,vmat,&
            ritzvec,reigval,res,hdiag,lcollapse)
    else if (ipre.eq.2) then
       ! Olsen's preconditioner
       errmsg='Olsen''s preconditioner is not yet implemented'
       call error_control
    endif

!----------------------------------------------------------------------
! Project the subspace onto the space orthogonal to the that spanned
! by the converged vectors
!----------------------------------------------------------------------
    if (ldeflate.and.nconv.gt.0) &
         call subspace_projection(matdim,maxvec,blocksize,&
         blocksize_curr,currdim,nconv,vmat,ritzvec,lcollapse)
    
!----------------------------------------------------------------------
! Orthogonalise the subspace vectors
!----------------------------------------------------------------------
    call qrortho(matdim,maxvec,blocksize_curr,currdim,vmat,lcollapse)

!----------------------------------------------------------------------    
! Update the dimension of the subspace
!----------------------------------------------------------------------
    if (lcollapse) then
       currdim=2*blocksize_curr
    else
       currdim=currdim+blocksize_curr
    endif
    
    return
    
  end subroutine subspace_expansion

!######################################################################
! deflate: removes converged vectors from the subspace
!######################################################################
  subroutine deflate(matdim,maxvec,blocksize,ritzvec,reigval,res,&
       norm,nstates,nstates_curr,currdim,maxvec_curr,blocksize_curr,&
       nconv,nconv_prev,tol)

    use constants

    implicit none

    ! Dimensions
    integer(is), intent(in)    :: matdim,maxvec,blocksize,nstates,&
                                  nconv,nconv_prev
    integer(is), intent(inout) :: blocksize_curr,maxvec_curr,&
                                  nstates_curr,currdim
    
    ! Ritz vectors
    real(dp), intent(inout)    :: ritzvec(matdim,blocksize)
        
    ! Eigenvalues of the Rayleigh matrix
    real(dp), intent(inout)    :: reigval(blocksize)

    ! Residuals
    real(dp), intent(inout)    :: res(matdim,blocksize)

    ! Residual norms
    real(dp), intent(inout)    :: norm(blocksize)

    ! Convergence threshold
    real(dp), intent(in)       :: tol
        
    ! Everything else
    integer(is)                :: indx,i,j,count
    integer(is)                :: convmap(nstates)
    real(dp), allocatable      :: swapvec(:)
    real(dp)                   :: swapval
      
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(swapvec(matdim))

!----------------------------------------------------------------------
! Rearrange the arrays holding the Ritz vectors, the Ritz values and
! the residual vectors s.t. the elements associated with the converged
! vectors are moved to the ends of the arrays
!----------------------------------------------------------------------
    ! Indices of the converged vectors
    count=0
    convmap=0
    do i=1,blocksize_curr
       if (norm(i).lt.tol) then
          count=count+1
          convmap(count)=i
       endif
    enddo

    ! Rearranegment of the arrays ritzvec, reigval and res
    count=0
    do i=1,nconv-nconv_prev
       
       count=count+1
       
       ! Index of the next column into which the next converged vector
       ! is to be moved
       indx=blocksize+1-nconv_prev-i
       
       ! Ritz vector
       swapvec=ritzvec(:,indx)
       ritzvec(:,indx)=ritzvec(:,convmap(count))
       ritzvec(:,convmap(count))=swapvec

       ! Ritz value
       swapval=reigval(indx)
       reigval(indx)=reigval(convmap(count))
       reigval(convmap(count))=swapval
       
       ! Residual vector
       swapvec=res(:,indx)
       res(:,indx)=res(:,convmap(count))
       res(:,convmap(count))=swapvec
       
    enddo

!----------------------------------------------------------------------
! Update dimensions
!----------------------------------------------------------------------
    blocksize_curr=blocksize_curr-(nconv-nconv_prev)
    maxvec_curr=maxvec_curr-(nconv-nconv_prev)
    nstates_curr=nstates_curr-(nconv-nconv_prev)
    currdim=currdim-(nconv-nconv_prev)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(swapvec)
    
    return
    
  end subroutine deflate
  
!######################################################################
! dpr: application of the diagonal preconditioned residue
!      preconditioner
!######################################################################
  subroutine dpr(matdim,maxvec,blocksize,blocksize_curr,currdim,vmat,&
       ritzvec,reigval,res,hdiag,lcollapse)

    use constants

    implicit none

    ! Dimensions
    integer(is), intent(in) :: matdim,maxvec,blocksize,&
                               blocksize_curr,currdim

    ! Subspace vectors
    real(dp), intent(inout) :: vmat(matdim,maxvec)

    ! Ritz vectors
    real(dp), intent(in)    :: ritzvec(matdim,blocksize)

    ! Residuals
    real(dp), intent(in)    :: res(matdim,blocksize)

    ! Eigenvalues of the Rayleigh matrix
    real(dp), intent(in)    :: reigval(blocksize)

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: hdiag(matdim)

    ! Subspace collapse flag
    logical, intent(in)     :: lcollapse

    ! Everything else
    integer(is)             :: i,j,ilbl,ilast
    real(dp), allocatable   :: tmpvec(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(tmpvec(matdim))

!----------------------------------------------------------------------
! Calculate the new subspace vectors
!----------------------------------------------------------------------
    if (lcollapse) then
       ! Collapse of the subspace
       vmat=0.0d0
       vmat(:,1:blocksize_curr)=ritzvec(:,1:blocksize_curr)
       ilast=blocksize_curr
    else
       ! Expansion of the subspace
       ilast=currdim
    endif

    ! Loop over new subspace vectors
    do i=1,blocksize_curr
       
       ! Index of the next vector
       ilbl=ilast+i
       
       ! Calculate the next vector
       tmpvec=0.0d0
       do j=1,matdim
          tmpvec(j)=1.0d0/(reigval(i)-hdiag(j))
       enddo
       do j=1,matdim
          vmat(j,ilbl)=res(j,i)*tmpvec(j)
       enddo

    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(tmpvec)
    
    return
    
  end subroutine dpr

!######################################################################
! subspace_projection: projection of the subspace vectors onto the
!                      orthogonal complement of the space spanned by
!                      the converged vectors
!######################################################################
  subroutine subspace_projection(matdim,maxvec,blocksize,&
       blocksize_curr,currdim,nconv,vmat,ritzvec,lcollapse)

    use constants
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: matdim,maxvec,blocksize,blocksize_curr,&
                               currdim,nconv
    
    ! Subspace vectors
    real(dp), intent(inout) :: vmat(matdim,maxvec)

    ! Ritz vectors
    real(dp), intent(inout) :: ritzvec(matdim,blocksize)

    ! Subspace collapse flag
    logical, intent(in)     :: lcollapse
    
    ! Everything else
    integer(is)             :: i,j,lower,upper,indx,nsubvec
    real(dp), allocatable   :: overlap(:,:)

!----------------------------------------------------------------------
! Set the lower and upper indices on the unconverged subspace vectors
! that need to be orthogonalised against the converged vectors
!----------------------------------------------------------------------
    if (lcollapse) then
       lower=1
       upper=2*blocksize_curr
    else
       lower=1
       upper=currdim+blocksize_curr
    endif
    
    nsubvec=upper-lower+1

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(overlap(nconv,nsubvec))

!----------------------------------------------------------------------
! Calculate the matrix of overlaps between the converged vectors and
! the unconverged subspace vectors
!----------------------------------------------------------------------
    call dgemm('T','N',nconv,nsubvec,matdim,1.0d0,&
         ritzvec(:,blocksize-nconv+1:blocksize),matdim,&
         vmat(:,lower:upper),matdim,0.0d0,overlap,nconv)
    
!----------------------------------------------------------------------
! Orthogonalise the unconverged subspace vectors againts to the
! converged vectors
!----------------------------------------------------------------------
    do i=lower,upper
       do j=1,nconv
          indx=blocksize-nconv+j
          vmat(:,i)=vmat(:,i)-ritzvec(:,indx)*overlap(j,i)
       enddo
    enddo
      
!---------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(overlap)
    
    return
    
  end subroutine subspace_projection
  
!######################################################################
! qrortho: orthogonalisation of the subspace vectors
!######################################################################
  subroutine qrortho(matdim,maxvec,blocksize_curr,currdim,vmat,&
       lcollapse)

    use constants
    use iomod
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: matdim,maxvec,blocksize_curr,currdim
    
    ! Subspace vectors
    real(dp), intent(inout) :: vmat(matdim,maxvec)

    ! Subspace collapse flag
    logical, intent(in)     :: lcollapse

    ! Everything else
    integer(is)             :: n,info
    real(dp), allocatable   :: tau(:),work(:)

!----------------------------------------------------------------------
! Orthogonalisation of the subspace vectors via a QR factorization
!----------------------------------------------------------------------
    if (lcollapse) then
       n=2*blocksize_curr
    else
       n=currdim+blocksize_curr
    endif

    allocate(tau(n))
    allocate(work(n))

    call dgeqrf(matdim,n,vmat(:,1:n),matdim,tau,work,n,info)
    if (info.ne.0) then
       errmsg='dqerf failed in subroutine qrortho'
       call error_control
    endif
    
    call dorgqr(matdim,n,n,vmat(:,1:n),matdim,tau,work,n,info)
    if (info.ne.0) then
       errmsg='dorgqr failed in subroutine qrortho'
       call error_control
    endif
    
    deallocate(tau)
    deallocate(work)
    
    return
    
  end subroutine qrortho

!######################################################################
! wreigenpairs: writes the eigenpairs to disk
!######################################################################
  subroutine wreigenpairs(irrep,vecscr,vecstem,matdim,blocksize,&
       reigval,ritzvec,nstates,nstates_curr,nconv_prev,ldeflate)

    use constants
    use bitglobal
    use utils
    use iomod
    
    implicit none

    ! Irrep number
    integer(is), intent(in)      :: irrep
    
    ! Eigenpair scratch file number and file stem
    integer(is), intent(out)     :: vecscr
    character(len=*), intent(in) :: vecstem

    ! Dimensions
    integer(is), intent(in)      :: matdim,blocksize
    integer(is), intent(in)      :: nstates,nstates_curr,nconv_prev
    
    ! Eigenpairs
    real(dp), intent(inout)          :: reigval(blocksize)
    real(dp), intent(inout)          :: ritzvec(matdim,blocksize)
                                 
    ! Deflation flag             
    logical, intent(in)          :: ldeflate

    ! Everything else
    integer                      :: i
    integer                      :: indx(nstates)
    integer(is)                  :: iscratch
    real(dp)                     :: esort(nstates)
    character(len=60)            :: vecfile
    character(len=2)             :: amult,airrep
    
!-----------------------------------------------------------------------
! Rearrangement of the ritzvec and reigval array such that the
! converged eigenpairs come first
!
! No. of converged eigenpairs "in storage": nconv_prev
! No. converged eigenpairs not "in storage": nstates_curr      
!
! N.B. This is only relevant if the subspace deflation was used
!----------------------------------------------------------------------
    if (ldeflate) then

       reigval(nstates_curr+1:nstates)=&
            reigval(blocksize-nconv_prev+1:blocksize)
       
       ritzvec(:,nstates_curr+1:nstates)=&
            ritzvec(:,blocksize-nconv_prev+1:blocksize)
       
    endif

!----------------------------------------------------------------------
! Indices of the eigenpairs in order of increasing energy
!
! Note that this is only relevant if subspace deflation was used, but
! makes no difference if not
!----------------------------------------------------------------------
    if (nstates.gt.1) then
       call dsortindxa1('A',nstates,reigval(1:nstates),indx)
    else
       indx(1)=1
    endif

!----------------------------------------------------------------------
! Add on E_SCF
!----------------------------------------------------------------------
    reigval=reigval+escf

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
    do i=1,nstates
       esort(i)=reigval(indx(i))
    enddo
    write(iscratch) esort(1:nstates)
    
    ! Eigenvectors
    do i=1,nstates
       write(iscratch) ritzvec(:,indx(i))
    enddo
     
    ! Close the scratch file
    close(iscratch)
        
    return
    
  end subroutine wreigenpairs
  
!######################################################################
  
end module blockdav
