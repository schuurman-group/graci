!**********************************************************************
! Sundry determinant manipulation routines
!**********************************************************************
module detfuncs

  implicit none

contains

!######################################################################
! truncate_wave_functions: Truncation of the input wave functions
!                          based on a minimum norm threshold
!######################################################################
  subroutine truncate_wave_functions(n_int,ndet,nroots,det,vec,&
       normthrsh,ndet_new,det_new,vec_new)

    use constants
    use global
    use utils
    use timing
    
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: n_int,ndet,nroots

    ! Determinant bit strings
    integer(ib), intent(in)  :: det(n_int,2,ndet)

    ! Eigenvectors
    real(dp), intent(in)     :: vec(ndet,nroots)

    ! Norm-based truncation threshold
    real(dp), intent(in)     :: normthrsh

    ! Truncated determinant and eigenvector arrays
    integer(is), intent(out) :: ndet_new
    integer(ib), allocatable :: det_new(:,:,:)
    real(dp), allocatable    :: vec_new(:,:)
    
    ! Sorting
    integer(is), allocatable :: indx(:)
    real(dp), allocatable    :: cabs(:)

    ! Surviving determinants
    integer(is), allocatable :: idet(:)

    ! Timing variables
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                                twall_end
    
    ! Everything else
    integer(is)              :: i,k,n
    real(dp)                 :: normsq,targ,diff
    real(dp), parameter      :: epsilon=1e-6_dp

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(indx(ndet))
    indx=0

    allocate(cabs(ndet))
    cabs=0.0d0

    allocate(idet(ndet))
    idet=0
    
!----------------------------------------------------------------------
! Determine the indices of the surviving determinants
!----------------------------------------------------------------------
    ! Target squared norm
    targ=normthrsh**2
    
    ! Loop over roots
    do i=1,nroots

       ! Sort the coefficients by absolute value
       cabs=abs(vec(:,i))
       call dsortindxa1('D',ndet,cabs,indx)

       ! Fill in the surviving determinants for this state
       normsq=0.0d0
       do k=1,ndet

          ! Update the squared norm
          normsq=normsq+vec(indx(k),i)**2

          ! Flag the determinant for survival
          idet(indx(k))=1
          
          ! Exit if:
          ! (1) we have hit the target squared norm, and;
          ! (2) this determinant is not degenerate with the next one
          diff = 10.*epsilon
          if(k.lt.ndet) diff=abs(vec(indx(k),i)-vec(indx(k+1),i))
          if (normsq >= targ .and. diff > epsilon) exit

       enddo
       
    enddo

!----------------------------------------------------------------------
! Set up the truncated determinant and eigenvector arrays
!----------------------------------------------------------------------
    !
    ! Number of surviving determinants
    !
    ndet_new=sum(idet)

    !
    ! Allocate arrays
    !
    allocate(det_new(n_int,2,ndet_new))
    det_new=0_ib
    allocate(vec_new(ndet_new,nroots))
    vec_new=0.0d0

    !
    ! Fill in the determinant bit strings
    !
    n=0
    ! Loop over the original set of determinants
    do k=1,ndet

       ! Are we at a surviving determinant?
       if (idet(k) == 1) then
          n=n+1
          det_new(:,:,n)=det(:,:,k)
       endif
       
    enddo

    !
    ! Fill in the eigenvectors
    !
    ! Loop over roots
    do i=1,nroots

       ! Surviving determinant counter
       n=0
       
       ! Loop over the original set of determinants
       do k=1,ndet

          ! Are we at a surviving determinant?
          if (idet(k) == 1) then
             n=n+1
             vec_new(n,i)=vec(k,i)
          endif
             
       enddo

    enddo

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    if (verbose) &
         call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'truncate_wave_functions')
    
    return
    
  end subroutine truncate_wave_functions

!######################################################################
! symm_ortho: Orthonormalisation of a set of wave functions using
!             Lowdin's symmetric orthogonalisation
!######################################################################
  subroutine symm_ortho(n_int,ndet,nroots,vec)

    use constants
    use utils
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: n_int,ndet,nroots

    ! Eigenvectors
    real(dp), intent(inout) :: vec(ndet,nroots)

    ! Everything else
    integer(is)             :: i,j
    real(dp)                :: norm
    real(dp), allocatable   :: Smat(:,:),Sinvsq(:,:),vec_ortho(:,:)
        
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(Smat(nroots,nroots), Sinvsq(nroots,nroots), &
         vec_ortho(ndet,nroots))
    Smat=0.0d0; Sinvsq=0.0d0; vec_ortho=0.0d0
    
!----------------------------------------------------------------------
! Normalisation
!----------------------------------------------------------------------
    do i=1,nroots
       norm=dot_product(vec(:,i),vec(:,i))
       norm=sqrt(norm)
       vec(:,i)=vec(:,i)/norm
    enddo
    
!----------------------------------------------------------------------
! Overlap matrix
!----------------------------------------------------------------------
    call dgemm('T','N',nroots,nroots,ndet,1.0d0,vec,ndet,vec,ndet,&
         0.0d0,Smat,nroots)

!----------------------------------------------------------------------
! Inverse square root of the overlap matrix
!----------------------------------------------------------------------
    call invsqrt_matrix(Smat,Sinvsq,nroots)

!----------------------------------------------------------------------
! Orthogonalisation
!----------------------------------------------------------------------
    call dgemm('N','N',ndet,nroots,nroots,1.0d0,vec,ndet,Sinvsq,&
         nroots,0.0d0,vec_ortho,ndet)

    vec=vec_ortho
    
    return
    
  end subroutine symm_ortho
    
!######################################################################
! get_nel: Determines the number of electrons (total, alpha and beta)
!          in a given determinant d
!######################################################################
  subroutine get_nel(n_int,d,nel,nel_alpha,nel_beta)

    use constants
    
    implicit none

    integer(is), intent(in)  :: n_int
    integer(ib), intent(in)  :: d(n_int,2)
    integer(is), intent(out) :: nel,nel_alpha,nel_beta

    integer(is)              :: k

    !
    ! Number of alpha electrons
    !
    nel_alpha=0
    do k=1,n_int
       nel_alpha=nel_alpha+popcnt(d(k,1))
    enddo
    
    !
    ! Number of beta electrons
    !
    nel_beta=0
    do k=1,n_int
       nel_beta=nel_beta+popcnt(d(k,2))
    enddo

    !
    ! Total number of electrons
    !
    nel=nel_alpha+nel_beta
    
    return
    
  end subroutine get_nel
  
!######################################################################
! mo_occ_string: Given a single alpha/beta string, determines the
!                indices of the occupied MOs
!######################################################################
  subroutine mo_occ_string(n_int,string,dim,nocc,occ)

    use constants
    
    implicit none

    ! Alpha or beta string
    integer(is), intent(in)  :: n_int
    integer(ib), intent(in)  :: string(n_int)

    ! No. occupied orbitals
    integer(is), intent(out) :: nocc

    ! Indices of the occupied orbitals
    integer(is), intent(in)  :: dim
    integer(is), intent(out) :: occ(dim)

    ! Everything else
    integer(is)              :: ic,k,ipos
    integer(ib)              :: h

    !
    ! Initialisation
    !
    occ=0
    nocc=0

    !
    ! Get the indices of the occupied orbitals
    !
    ! Orbital counter
    ic=1

    ! Loop over bit string blocks
    do k=1,n_int

       ! Initialise the work array
       h=string(k)

       ! Get the occupied orbital indices for this block
       do while (h /= 0_ib)

          ! Number of trailing zeros left in h
          ipos=trailz(h)

          ! Index of the next occupied orbital
          occ(ic)=1+ipos+(k-1)*n_bits

          ! Clear the bits up to the occupied orbital in h
          h=ibclr(h,ipos)

          ! Increment the orbital counter
          ic=ic+1
          
       enddo
          
    enddo

    ! No. occupied orbitals
    nocc=ic-1
    
    return
    
  end subroutine mo_occ_string
  
!######################################################################
! anihilate_electron_string: annihilates the electron in an input
!                            alpha/beta string
!######################################################################
  function annihilate_electron_string(n_int,string,imo) &
       result(hole_string)

    use constants

    implicit none

    ! Function result
    integer(is), intent(in) :: n_int
    integer(ib)             :: hole_string(n_int)

    ! Input alpha/beta string
    integer(ib), intent(in) :: string(n_int)

    ! Index of the MO to be annihilated
    integer(is), intent(in) :: imo

    ! Everything else
    integer(is)             :: k,i

    ! Block index
    k=(imo-1)/n_bits+1

    ! Orbital position with the block
    i=imo-1-(k-1)*n_bits

    ! Annihilate the electron
    hole_string=ibclr(string,i)
    
    return
    
  end function annihilate_electron_string
  
!######################################################################
  
end module detfuncs
