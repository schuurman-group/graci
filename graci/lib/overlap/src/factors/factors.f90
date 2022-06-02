!**********************************************************************
! Routines for the pre-calculation of the unique alpha/beta factors,
! i.e., the unique determinants of the bra/ket spin orbital overlaps
!**********************************************************************
module factors

  implicit none

contains

!######################################################################
! get_factors: Top level routine for the calculation of the complete
!              set of alpha/beta factors given sets of unique bra and
!              ket alpha/beta strings
!######################################################################
  subroutine get_factors(neB,neK,nstringB,nstringK,stringB,stringK,fac)

    use constants
    use global, only: n_intB,n_intK,nmoB,nmoK,smo
    use detfuncs
    use timing
    
    implicit none

    ! Number of electrons in each bra and ket string
    integer(is), intent(in)  :: neB,neK
    
    ! Dimensions
    integer(is), intent(in)  :: nstringB,nstringK

    ! Unique alpha/beta strings
    integer(ib), intent(in)  :: stringB(n_intB,nstringB)
    integer(ib), intent(in)  :: stringK(n_intK,nstringK)    

    ! Alpha/beta factors
    real(dp), intent(out)    :: fac(nstringB,nstringK)

    ! Occupied MOs
    integer(is)              :: noccB,noccK
    integer(is), allocatable :: occB(:),occK(:)

    ! Work arrays
    real(dp), allocatable    :: work(:,:)
    integer(is), allocatable :: ipiv(:)

    ! Timing variables
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
         twall_end
    
    ! Everything else
    integer(is)              :: ibra,iket,m,n,info
    real(dp)                 :: determinant

!*** WE ARE PASSING neB AND neK BUT ALSO ASSUMING THAT THEY ARE EQUAL***

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(occB(neB), occK(neK))
    occB=0; occK=0

    allocate(work(neB,neK))
    work=0.0d0

    allocate(ipiv(neB))
    ipiv=0

!----------------------------------------------------------------------
! Compute the unique factors
!----------------------------------------------------------------------
    ! Loop over ket strings
    do iket=1,nstringK

       ! Get the ket occupied MO indices
       call mo_occ_string(n_intK,stringK(:,iket),neK,noccK,occK)

       ! Loop over bra strings
       do ibra=1,nstringB

          ! Get the bra occupied MO indices
          call mo_occ_string(n_intB,stringB(:,ibra),neB,noccB,occB)
          
          ! Fill in the matrix of occupied bra-ket MO overlaps
          do m=1,neK
             do n=1,neB
                work(n,m)=smo(occB(n),occK(m))
             enddo
          enddo

          ! LU decomposition of the matrix of MO overlaps
          call dgetrf(neB,neB,work,neB,ipiv,info)
          if (info /= 0) then
             write(6,'(/,x,a)') &
                  'Error in get_factors: LU decomposition failed'
             stop
          endif

          ! Determinant of the matrix of MO overlaps
          determinant=1.0d0
          do n=1,neB
             if (ipiv(n) /= n) then
                determinant=-determinant*work(n,n)
             else
                determinant=determinant*work(n,n)
             endif
          enddo

          ! Fill in the factor array
          fac(ibra,iket)=determinant
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'get_factors')
    
    return
    
  end subroutine get_factors
  
!######################################################################
  
end module factors
