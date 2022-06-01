!**********************************************************************
! Routines for the sorting of the determinant bit string arrays into
! alpha- and beta-major order, determining the unique alpha and beta
! strings, etc.
!**********************************************************************
module detsort

  implicit none

contains

!######################################################################
! unique_strings: top level routine for the determination of the unique
!                 alpha and beta strings
!######################################################################
  subroutine unique_strings(n_int,ndet,det)

    use constants
    
    implicit none

    ! Determinant bit strings
    integer(is), intent(in)  :: n_int,ndet
    integer(ib), intent(in)  :: det(n_int,2,ndet)

    ! Working arrays
    integer(ib), allocatable :: det_sort(:,:,:),work(:,:,:)

    ! Everything else
    integer(is)              :: i
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(det_sort(n_int,2,ndet))
    det_sort=0_ib

    allocate(work(n_int,2,ndet))
    work=0_ib
    
!----------------------------------------------------------------------
! Determine the uniqe alpha strings
!----------------------------------------------------------------------
    ! Put the determinants into alpha-major order
    det_sort=det
    call radix_sort(n_int, ndet, det_sort, work, 1)


    ! We are going to have to change radix_sort s.t. it returns
    ! a new-to-old index mapping array...
    
    STOP

    
    return
    
  end subroutine unique_strings

!######################################################################
! radix_sort: radix sorting of a given set of determinants into either
!             alpha- or beta-major order
!######################################################################
  subroutine radix_sort(n_int,ndet,det_list,scr,srt_indx)

    use constants
    
    implicit none

    ! Number of integers needed to represent each alpha/beta string
    integer(is), intent(in)            :: n_int
    
    ! Number of dets in the list
    integer(is), intent(in)            :: ndet

    ! Which index to sort, 1=alpha, 2=beta
    integer(is), intent(in)            :: srt_indx

    ! List of dets as bit strings
    integer(ib), target, intent(inout) :: det_list(0:n_int-1,0:1,0:ndet-1)

    ! Currently scratch is same size as det list
    ! This is twice as big as it needs to be, but simplest for now
    integer(ib), target, intent(inout) :: scr(0:n_int-1,0:1,0:ndet-1)

    ! The pointer to the data (or scratch)
    integer(ib), pointer               :: ptra(:,:,:)

    ! The pointer to the scratch (or data)           
    integer(ib), pointer               :: ptrb(:,:,:)

    ! ib integer bit string
    integer(ib)                        :: a,b

    ! hold the array indices for each bit
    integer(is), allocatable           :: indx_table(:,:)

    ! counter variables
    integer(is)                        :: id,ii,ibt,ibi

    ! index variables 
    integer(is)                        :: m,n

    ! how many bytes per ib integer
    integer(is)                        :: n_b

    ! alternate data/scratch pointer target
    logical                            :: det_sorted
    
    det_sorted = .true.

    ! number of bytes per ib int
    n_b = sizeof(a)

    ! we will treat the tuple of integers as a single n_int*nb quantity.
    ! NOTE: future optimization would discard the the unused bits in the
    ! 'last' ib integer in the n_int tuple

    ! table for mapping radix indices
    allocate(indx_table(0:255, 0:n_int*n_b-1))
    
    ptra => det_list
    ptrb => scr
    
    ! First step is to create a histogram of the contents of each
    ! bucket for each element in the list. Once we know what we're
    ! dealing with element-wise, we can create an indexing array for
    ! each of the buckets 
    indx_table = 0
    do id = 0, ndet-1
       ibi = 0
       do ii = 0, n_int-1
          a = ptra(ii, srt_indx-1, id)
          do ibt = 0, n_b-1
             b = iand(a , z'FF')
             indx_table(b, ibi) = indx_table(b, ibi) + 1
             a = shiftr(a , 8_ib)
             ibi = ibi + 1
          enddo
       enddo
    enddo
    
    ! Now step through the each bucket in the histogram to assign a
    ! starting index for the byte/radix value
    ibi = 0
    do ii = 0, n_int-1
       do ibt = 0, n_b-1                                         
          m = 0
          do b = 0, 255
             n = indx_table(b, ibi)
             indx_table(b, ibi) = m
             m = m + n
          enddo
          ibi = ibi + 1
       enddo
    enddo
    
    ! This is the entirety of the radix sort (by LSB)
    ibi = 0
    do ii = 0,n_int-1
       do ibt = 0, n_b-1       
          
          do id = 0, ndet-1
             a = ptra(ii, srt_indx-1, id)
             b = iand(shiftr(a, 8_is*ibt), z'FF')
             ptrb(:, :, indx_table(b, ibi)) = ptra(:, :, id)
             indx_table(b, ibi)             = indx_table(b, ibi) + 1
          enddo
          
          ibi = ibi + 1

          ! Instead of copying back from the temp values swap the array
          ! pointers
          det_sorted = .not.det_sorted
          if( det_sorted ) then
             ! if det_list is current sorted list, set a to det_list
             ptra => det_list                                     
             ptrb => scr                                   
          else
             ! else if scr is the current partially sorted list, set a
             ! to scr
             ptra => scr
             ptrb => det_list
          endif

       enddo
    enddo

    ! if sorted array is scratch, copy to det_list
    if ( .not.det_sorted ) then
       det_list = scr
    endif

    return

  end subroutine radix_sort
  
!######################################################################
  
end module detsort
