!**********************************************************************
! Routines for the sorting of the determinant bit string arrays into
! alpha- and beta-major order, determining the unique alpha and beta
! strings, etc.
!**********************************************************************
module detsort

  implicit none

contains

!######################################################################
! det_sorting: Double sorting of a given set of input determinants
!
!              The determinants are subjected to the following:
!
!              (1) sorting into alpha-major order, followed by;
!
!              (2) sorting of each alpha-string block by beta strings
!
!              During the course of this, the unique alpha and beta
!              strings are also determined, as well as the
!              (sorted) determinant-to-unique beta string mapping
!######################################################################
  subroutine det_sorting(n_int,ndet,nroots,det,vec,nalpha,nbeta,alpha,&
       beta,offset)

    use constants

    implicit none

    ! Determinant bit strings
    integer(is), intent(in)    :: n_int,ndet
    integer(ib), intent(inout) :: det(n_int,2,ndet)

    ! Eigenvectors
    integer(is)                :: nroots
    real(dp), intent(inout)    :: vec(ndet,nroots)
    
    ! No. unique alpha and beta strings
    integer(is), intent(out)   :: nalpha,nbeta

    ! Alpha string offsets
    integer(is), allocatable   :: offset(:)
    
    ! Unique alpha and beta strings
    integer(ib), allocatable   :: alpha(:,:),beta(:,:)
    
    ! Working arrays
    integer(ib), allocatable   :: det_sort(:,:,:),dwork(:,:,:)
    integer(is), allocatable   :: mwork(:)
    real(dp), allocatable      :: vwork(:)
    
    ! Sorted-to-unsorted determinant index mapping array
    integer(is), allocatable   :: imap(:)
    
    ! Everything else
    integer(is)                :: i,j,k,n,istart,iend,nd

!----------------------------------------------------------------------
! (1) Put the determinant and eigenvector arrays into alpha-major
!     order
!----------------------------------------------------------------------
    ! Allocate work arrays
    allocate(dwork(n_int,2,ndet), mwork(ndet), imap(ndet), vwork(ndet))
    dwork=0_ib; mwork=0; imap=0; vwork=0.0d0

    ! Put the determinants into alpha-major order
    call radix_sort(n_int,ndet,det,imap,dwork,mwork,1)

    ! Put the eigenvectors into alpha-major order
    do i=1,nroots
       vwork=vec(:,i)
       do j=1,ndet
          vec(j,i)=vwork(imap(j))
       enddo
    enddo

    ! Deallocate work arrays
    deallocate(dwork,mwork,imap,vwork)
    
!----------------------------------------------------------------------
! (2) Using the alpha-major ordered determinants, determine the unique
!     alpha strings and their offsets
!----------------------------------------------------------------------
    ! Get the no. unique alpha strings
    nalpha=nunique(n_int,ndet,det,1)

    ! Allocate arrays
    allocate(alpha(n_int,nalpha))
    alpha=0_ib
    allocate(offset(nalpha+1))
    offset=0

    ! Fill in the unique alpha strings and their offsets
    call fill_unique(n_int,ndet,nalpha,det,alpha,offset,1)

!----------------------------------------------------------------------
! (3) Sort each alpha-string block of determinants by their beta string
!----------------------------------------------------------------------
    ! Loop over alpha string blocks
    do n=1,nalpha

       ! Start and end of this block
       istart=offset(n)
       iend=offset(n+1)-1
       
       ! Number of determinants in this block
       nd=iend-istart+1

       ! Allocate arrays
       allocate(det_sort(n_int,2,nd), dwork(n_int,2,nd), mwork(nd), &
            imap(nd), vwork(nd))
       dwork=0_ib; mwork=0; imap=0; vwork=0.0d0
       
       ! Sort this block by beta string
       det_sort=det(:,:,istart:iend)
       call radix_sort(n_int,nd,det_sort,imap,dwork,mwork,2)

       ! Rearrange this block of the determinants
       det(:,:,istart:iend)=det_sort
       
       ! Rearrange this block of the eigenvectors
       do i=1,nroots
          vwork=vec(istart:iend,i)
          do j=1,nd
             vec(istart+j-1,i)=vwork(imap(j))
          enddo
       enddo
       
       ! Deallocate arrays
       deallocate(det_sort,dwork,mwork,imap,vwork)
       
    enddo

!----------------------------------------------------------------------
! (4) Using the now double-sorted array of determinants, determine:
!     (i)  the unique beta strings
!     (ii) the determinant-to-unique beta string mapping
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(dwork(n_int,2,ndet), mwork(ndet), imap(ndet), vwork(ndet))
    dwork=0_ib; mwork=0; imap=0; vwork=0.0d0
    
    ! Put a copy of the double-sorted determinants into beta-major
    ! order
    det_sort=det
    call radix_sort(n_int,ndet,det_sort,imap,dwork,mwork,2)

    ! Get the no. unique beta strings
    nbeta=nunique(n_int,ndet,det_sort,2)

    ! Fill in the unique alpha strings (passing a dummy offset array)
    
    
    
    ! Deallocate arrays
    deallocate(det_sort,dwork,mwork,imap,vwork)
    
    return
    
  end subroutine det_sorting
  
!!######################################################################
!! unique_strings: Top level routine for the determination of the unique
!!                 alpha and beta strings
!!######################################################################
!  subroutine unique_strings(n_int,ndet,det,nalpha,nbeta,alpha,beta)
!
!    use constants
!    
!    implicit none
!
!    ! Determinant bit strings
!    integer(is), intent(in)  :: n_int,ndet
!    integer(ib), intent(in)  :: det(n_int,2,ndet)
!
!    ! No. unique alpha and beta strings
!    integer(is), intent(out) :: nalpha,nbeta
!
!    ! Unique alpha and beta strings
!    integer(ib), allocatable :: alpha(:,:),beta(:,:)
!    
!    ! Working arrays
!    integer(ib), allocatable :: det_sort(:,:,:),dwork(:,:,:)
!    integer(is), allocatable :: mwork(:)
!
!    ! Sorted-to-unsorted determinant index mapping array
!    integer(is), allocatable :: imap(:)
!    
!    ! Everything else
!    integer(is)              :: i
!    
!!----------------------------------------------------------------------
!! Allocate arrays
!!----------------------------------------------------------------------
!    allocate(det_sort(n_int,2,ndet))
!    det_sort=0_ib
!
!    allocate(dwork(n_int,2,ndet))
!    dwork=0_ib
!
!    allocate(imap(ndet))
!    imap=0
!
!    allocate(mwork(ndet))
!    mwork=0
!    
!!----------------------------------------------------------------------
!! Determine the uniqe alpha strings
!!----------------------------------------------------------------------
!    ! Put the determinants into alpha-major order
!    det_sort=det
!    call radix_sort(n_int,ndet,det_sort,imap,dwork,mwork,1)
!
!    ! Get the no. unique alpha strings
!    nalpha=nunique(n_int,ndet,det_sort,1)
!
!    ! Allocate the unique alpha string array
!    allocate(alpha(n_int,nalpha))
!    alpha=0_ib
!
!    ! Fill in the unique alpha strings
!    call fill_unique(n_int,ndet,nalpha,det_sort,alpha,1)
!    
!!----------------------------------------------------------------------
!! Determine the uniqe beta strings
!!----------------------------------------------------------------------
!    ! Put the determinants into beta-major order
!    det_sort=det
!    call radix_sort(n_int,ndet,det_sort,imap,dwork,mwork,2)
!
!    ! Get the no. unique beta strings
!    nbeta=nunique(n_int,ndet,det_sort,2)
!
!    ! Allocate the unique beta string array
!    allocate(beta(n_int,nbeta))
!    beta=0_ib
!
!    ! Fill in the unique beta strings
!    call fill_unique(n_int,ndet,nbeta,det_sort,beta,2)
!
!    return
!    
!  end subroutine unique_strings

!######################################################################
! radix_sort: Radix sorting (in base 256) of a given set of
!             determinants into either alpha- or beta-major order, as
!             determined by the value of srt_indx (1<-> alpha-major,
!             2<-> beta-major)
!             Also returns the array imap of sorted-to-unsorted
!             determinant index mappings
!######################################################################
  subroutine radix_sort(n_int,ndet,det_list,imap,dscr,mscr,srt_indx)

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
    integer(ib), target, intent(inout) :: dscr(0:n_int-1,0:1,0:ndet-1)
    
    ! The pointer to the data (or scratch)
    integer(ib), pointer               :: ptrda(:,:,:)

    ! The pointer to the scratch (or data)           
    integer(ib), pointer               :: ptrdb(:,:,:)

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

    ! Sorted-to-unsorted det index mapping
    integer(is), target, intent(out)   :: imap(0:ndet-1)
    integer(is), target                :: mscr(0:ndet-1)
    integer(is), pointer               :: ptrma(:),ptrmb(:)

    det_sorted = .true.

    !
    ! Initialise the sorted-to-unsorted index mapping array
    !
    do id = 0, ndet-1
       imap(id) = id+1
    enddo

    !
    ! Number of bytes per ib int
    !
    n_b = sizeof(a)

    
    ! we will treat the tuple of integers as a single n_int*nb quantity.
    ! NOTE: future optimization would discard the the unused bits in the
    ! 'last' ib integer in the n_int tuple

    !
    ! Table for mapping radix indices
    !
    allocate(indx_table(0:255, 0:n_int*n_b-1))

    !
    ! Pointers to the det list, mapping and scratch arrays
    !
    ptrda => det_list
    ptrdb => dscr
    ptrma => imap
    ptrmb => mscr
    
    !
    ! First step is to create a histogram of the contents of each
    ! bucket for each element in the list. Once we know what we're
    ! dealing with element-wise, we can create an indexing array for
    ! each of the buckets
    !
    indx_table = 0
    ! Loop over determinants
    do id = 0, ndet-1
       ibi = 0
       ! Loop over blocks
       do ii = 0, n_int-1
          a = ptrda(ii, srt_indx-1, id)
          ! Loop over bytes in this block
          do ibt = 0, n_b-1
             b = iand(a , z'FF')
             indx_table(b, ibi) = indx_table(b, ibi) + 1
             a = shiftr(a , 8_ib)
             ibi = ibi + 1
          enddo
       enddo
    enddo
    
    !
    ! Now step through the each bucket in the histogram to assign a
    ! starting index for the byte/radix value
    !
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

    !
    ! This is the entirety of the radix sort (by LSB)
    !
    ibi = 0
    do ii = 0,n_int-1
       do ibt = 0, n_b-1       
          
          do id = 0, ndet-1
             a = ptrda(ii, srt_indx-1, id)
             b = iand(shiftr(a, 8_is*ibt), z'FF')
             ptrdb(:, :, indx_table(b, ibi)) = ptrda(:, :, id)
             ptrmb(indx_table(b, ibi))       = ptrma(id)
             indx_table(b, ibi)              = indx_table(b, ibi) + 1
          enddo
          
          ibi = ibi + 1

          ! Instead of copying back from the temp values swap the array
          ! pointers
          det_sorted = .not.det_sorted
          if( det_sorted ) then
             ! if det_list is current sorted list, set a to det_list
             ptrda => det_list
             ptrdb => dscr
             ptrma => imap
             ptrmb => mscr
          else
             ! else if scr is the current partially sorted list, set a
             ! to scr
             ptrda => dscr
             ptrdb => det_list
             ptrma => mscr
             ptrmb => imap
          endif

       enddo
    enddo

    ! if sorted array is scratch, copy to det_list
    if ( .not.det_sorted ) then
       det_list = dscr
       imap     = mscr
    endif

    return

  end subroutine radix_sort
  
!######################################################################
! nunique: Determines the number of unique alpha (ispin=1) or beta
!          (ispin=2) strings in an array of determinants that have
!          been ***pre-sorted*** into alpha- or beta-major order
!######################################################################
  function nunique(n_int,ndet,det,ispin)

    use constants
    
    implicit none

    integer(is)             :: nunique
    integer(is), intent(in) :: n_int,ndet
    integer(ib), intent(in) :: det(n_int,2,ndet)
    integer(is), intent(in) :: ispin

    integer(is)             :: i
    
    nunique=1

    do i=2,ndet
       if (.not. same_string(n_int,det(:,ispin,i), det(:,ispin,i-1))) &
            nunique=nunique+1
    enddo
    
    return
    
  end function nunique

!######################################################################
! same_string: Returns .true. (.false.) if the two input alpha/beta
!              strings are the same (not the same)
!######################################################################  
  function same_string(n_int,string1,string2)

    use constants

    implicit none

    logical                 :: same_string
    integer(is), intent(in) :: n_int
    integer(ib), intent(in) :: string1(n_int),string2(n_int)
    integer(is)             :: k

    same_string=.true.

    do k=1,n_int
       if (string1(k) /= string2(k)) then
          same_string=.false.
          exit
       endif
    enddo
    
    return
    
  end function same_string
  
!######################################################################
! fill_unique: Returns the unique alpha (ispin=1) or beta (ispin=2)
!              strings given an array of determinants, det, that has
!              been ***pre-sorted*** into alpha- or beta-major order
!
!              Also fills in the alpha string offset array:
!
!              offset(i): starting point in the det array of the i'th
!                         unique ispin-spin string
!######################################################################
  subroutine fill_unique(n_int,ndet,nunique,det,string,offset,ispin)

    use constants

    implicit none

    ! Dimensions
    integer(is), intent(in)  :: n_int,ndet,nunique

    ! Pre-sorted determinant bit strings
    integer(ib), intent(in)  :: det(n_int,2,ndet)

    ! Unique ispin-spin strings
    integer(is), intent(in)  :: ispin
    integer(ib), intent(out) :: string(n_int,nunique)

    ! Offsets
    integer(is), intent(out) :: offset(nunique+1)
    
    ! Everything else
    integer(is)              :: i,n

    !
    ! First string
    !
    string(:,1)=det(:,ispin,1)
    offset(1)=1
    
    !
    ! Remaining strings
    !
    n=1
    do i=2,ndet
       if (.not. &
            same_string(n_int,det(:,ispin,i), det(:,ispin,i-1))) then
          n=n+1
          string(:,n)=det(:,ispin,i)
          offset(n)=i
       endif
    enddo

    offset(nunique+1)=ndet+1
    
    return
    
  end subroutine fill_unique
    
!######################################################################
  
end module detsort
