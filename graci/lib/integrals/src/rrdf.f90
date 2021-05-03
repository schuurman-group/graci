module rrdf

  use accuracy

contains

!----------------------------------------------------------------------
! Rank reduction of the 3-index density fitting integrals
!----------------------------------------------------------------------
  subroutine rank_reduce(n_aux,n_ij,b,fac)

    use timing
    
    implicit none

    integer(ik), intent(inout) :: n_aux  ! No. aux basis functions
    integer(ik), intent(in)    :: n_ij   ! No. pairs MOs
    real(drk), allocatable     :: b(:,:) ! DF vectors
    integer(ik), intent(in)    :: fac    ! Compression factor
    
    integer(ik)                :: n_new,workdim,info,i
    integer(ik), allocatable   :: iwork(:)
    real(drk), allocatable     :: U(:,:),VT(:,:),S(:),work(:)
    real(drk)                  :: tmp(1,1),tmp1(1,1)
    real(drk)                  :: tcpu_start,tcpu_end,twall_start,&
                                  twall_end

    !
    ! Output the compression factor
    !
    write(6,'(/,x,a,x,i0)') 'RR-DF compression factor:',fac
    
    !
    ! Allocate arrays
    !
    workdim = 3*n_aux + max(n_ij, 5*n_aux**2 + 4*n_aux)

    allocate(U(n_aux,n_aux), VT(n_aux,n_ij), S(n_aux), work(workdim), &
         iwork(8*n_aux))
    
    !
    ! Start timing
    !
    call get_times(twall_start,tcpu_start)
    
    !
    ! SVD of of the DF 3-index integral tensor
    VT=b

    call dgesdd('O',n_aux,n_ij,VT,n_aux,S,U,n_aux,tmp,1,work,workdim,&
         iwork,info)
    
    if (info /= 0) then
       write(6,'(/,2x,a)') 'SVD of b failed in subroutine rank_reduce'
       stop
    endif

    !
    ! Rank reduced DF integrals 
    !
    n_new=int(real(n_aux)/real(fac))
    n_aux = n_new

    deallocate(b)
    allocate(b(n_aux,n_ij))

    do i=1,n_new
       b(i,:) = VT(i,:) * S(i)
    enddo

    !
    ! Finish timing
    !
    call get_times(twall_end,tcpu_end)
    call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'rank_reduce')
    flush(6)
    
    return
    
  end subroutine rank_reduce
    
end module rrdf
