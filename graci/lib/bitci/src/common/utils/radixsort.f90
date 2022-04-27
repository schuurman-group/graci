module radixsort
   use constants
   use detutils

   implicit none

   contains

   subroutine radix_sort(ndet, det_list , scr, srt_indx)
      integer(is), intent(in)              ::  ndet     ! number of dets in the list
      integer(is), intent(in)              ::  srt_indx ! which index to sort, 1=alpha, 2=beta
      integer(ib), target,  intent(inout)  ::  det_list(0:n_int-1, 0:1, 0:ndet-1) ! list of dets as bit strings
      integer(ib), target,  intent(inout)  ::       scr(0:n_int-1, 0:1, 0:ndet-1) ! currently scratch is same size
                                                                                  ! as det list. This is twice as big
                                                                                  ! as it needs to be, but simplest for now

      ! local variables
      integer(ib), pointer                 ::  ptra(:,:,:)        ! The pointer to the data (or scratch)
      integer(ib), pointer                 ::  ptrb(:,:,:)        ! The pointer to the scratch (or data)           
      integer(ib)                          ::  a, b              ! ib integer bit string
      integer(is)                          ::  id, ii, ib, ibi, i, m, n    ! counter variables 
      integer(is), allocatable             ::  indx_table(:,:)   ! hold the array indices for each bit
      logical                              ::  det_sorted         ! alternate data/scratch pointer target
      integer(is)                          ::  n_b                ! how many bytes per ib integer

      !
      det_sorted = .true.                                     
      n_b        = sizeof(a)                                   ! number of byte per ib int

      ! we will treat the tuple of integers as a single n_int*nb quantity.
      ! NOTE: future optimization would discard the the unused bits in the 'last' ib integer
      ! in the n_int tuple (i.e. n_int*nb % nmo)
      allocate(indx_table(0:255, 0:n_int*n_b-1))           ! table for mapping radix indices

      ptra => det_list
      ptrb => scr
      
      ! First step is to create a histogram of the contents of each bucket for each element
      ! in the list. Once we know what we're dealing with element-wise, we can 
      ! create an indexing array for each of the buckets 
      indx_table = 0
      do id = 0, ndet-1
         ibi = 0
         do ii = 0, n_int-1
             a = ptra(ii, srt_indx-1, id)
             do ib = 0, n_b-1
                b = iand(a , z'FF')
                indx_table(b, ibi) = indx_table(b, ibi) + 1
                a = shiftr(a , 8_ib)
                ibi = ibi + 1
             enddo
         enddo
      enddo

      ! Now step through the each bucket in the histogram to assign a starting index
      ibi = 0
      do ii = 0, n_int-1
         do ib = 0, n_b-1                                         
             m = 0
             do b = 0, 255
                n = indx_table(b, ibi)
                indx_table(b, ibi) = m
                m = m + n   ! this defines the stride through the index array for a single iteration of the radix value
             enddo
             ibi = ibi + 1
         enddo
      enddo

      ! This is the entirety of the radix sort (by LSB)
      ibi = 0
      do ii = 0,n_int-1
         do ib = 0, n_b-1       

            do id = 0, ndet-1
               a = ptra(ii, srt_indx-1, id)
               b = iand(shiftr(a, 8*ib), z'FF')
               ptrb(:, :, indx_table(b, ibi)) = ptra(:, :, id)
               indx_table(b, ibi)             = indx_table(b, ibi) + 1
            enddo

            ibi = ibi + 1

            ! Instead of copying back from the temp values swap the array pointers
            det_sorted = .not.det_sorted
            if( det_sorted ) then  ! if det_list is current sorted list, set a to det_list
               ptra => det_list                                     
               ptrb => scr                                   
            else                   ! else if scr is the current partially sorted list,  set a to scr
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

end module radixsort


