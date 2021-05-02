module naf

  use accuracy

contains

!----------------------------------------------------------------------
! Rank reduction of the 3-index density fitting
! integrals using the 'Natural Auxilary Function' procedure detailed
! in J. Chem. Phys. 141, 244113 (2014)
!----------------------------------------------------------------------
  subroutine rank_reduce(n_aux,n_ij,b,thrsh)

    implicit none

    integer(ik), intent(inout) :: n_aux  ! No. aux basis functions
    integer(ik), intent(in)    :: n_ij   ! No. pairs MOs
    real(drk), allocatable     :: b(:,:) ! DF vectors
    real(drk), intent(in)      :: thrsh  ! Truncation threshold

    integer                    :: info,n_new,i,k
    real(drk), allocatable     :: W(:,:),eigval(:),work(:)
    real(drk), allocatable     :: b_new(:,:)
    
    !
    ! Allocate arrays
    !
    allocate(W(n_aux,n_aux))
    allocate(eigval(n_aux))
    allocate(work(3*n_aux))
    
    !
    ! Calculation of W = b b^T
    !
    call dgemm('N','T',n_aux,n_aux,n_ij,1.0d0,b,n_aux,b,n_aux,0.0d0,&
         W,n_aux)

    !
    ! Diagonalisation of W
    !
    call dsyev('V','U',n_aux,W,n_aux,eigval,work,3*n_aux,info)
    if (info /= 0) then
       write(6,'(/,2x,a)') &
            'Diagonalisation of bb^T failed in subroutine rank_reduce'
       stop
    endif

    !
    ! New auxiliary basis size
    !
    n_new=0

    do i=1,n_aux
       if (eigval(i)**2 > thrsh) n_new = n_new+1
    enddo
    
    write(6,'(/,x,a,x,F10.7)') 'DF rank reduction factor:',&
         dble(n_aux)/dble(n_new)
    
    !
    ! Transformation to the truncated NAF basis
    !
    allocate(b_new(n_new,n_ij))

    k = n_aux - n_new + 1

    ! To do: swap this out for a dgemm call
    b_new = matmul(transpose(W(:,k:n_aux)),b)
    
    deallocate(b)
    
    n_aux = n_new
    
    allocate(b(n_aux,n_ij))
    
    b = b_new
    
    return
    
  end subroutine rank_reduce
    
end module naf
