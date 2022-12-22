!**********************************************************************
! Routines for the pre-calculation of the unique alpha/beta factors,
! i.e., the unique determinants of the bra/ket spin orbital overlaps
!**********************************************************************
module factors

  implicit none

contains

!######################################################################
! get_all_factors: Calculation of the complete set of alpha/beta
!                  factors given sets of unique bra and ket alpha/beta
!                  strings
!######################################################################
  subroutine get_all_factors(nelX,nstringB,nstringK,stringB,stringK,&
       fac)

    use constants
    use global, only: n_intB,n_intK,nmoB,nmoK,smo,hthrsh,verbose
    use detfuncs
    
    implicit none

    ! Number of electrons in each bra and ket string
    ! (these have to be equal for an overlap calculation)
    integer(is), intent(in)  :: nelX
    
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

    ! Everything else
    integer(is)              :: ibra,iket,m,n
    real(dp)                 :: determinant

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(occB(nelX), occK(nelX))
    occB=0; occK=0

    allocate(work(nelX,nelX))
    work=0.0d0

    allocate(ipiv(nelX))
    ipiv=0

!----------------------------------------------------------------------
! Compute the unique factors
!----------------------------------------------------------------------
    ! Loop over ket strings
    do iket=1,nstringK

       ! Get the ket occupied MO indices
       call mo_occ_string(n_intK,stringK(:,iket),nelX,noccK,occK)

       ! Loop over bra strings
       do ibra=1,nstringB

          ! Get the bra occupied MO indices
          call mo_occ_string(n_intB,stringB(:,ibra),nelX,noccB,occB)
       
          ! Fill in the matrix of occupied bra-ket MO overlaps
          do m=1,nelX
             do n=1,nelX
                work(n,m)=smo(occB(n),occK(m))
             enddo
          enddo

          ! Determinant screening
          if (hadamard_bound(nelX,work) < hthrsh) then
             fac(ibra,iket)=0.0d0
             cycle
          endif
             
          ! Determinant of the matrix of MO overlaps
          fac(ibra,iket)=ludet(nelX,work,ipiv)
          
       enddo
       
    enddo

    return
    
  end subroutine get_all_factors

!######################################################################
! get_one_factor: Calculation of a single alpha/beta factor given a
!                 pair of alpha/beta strings and corresponding MO
!                 occupations  
!######################################################################
  subroutine get_one_factor(nelX,stringB,stringK,occB,occK,work,ipiv,&
       fac)

    use constants
    use global, only: n_intB,n_intK,nmoB,nmoK,smo,hthrsh
    use detfuncs
    use timing

    implicit none
    
    ! Number of electrons in string
    integer(is), intent(in)  :: nelX
    
    ! Unique alpha/beta strings
    integer(ib), intent(in)  :: stringB(n_intB)
    integer(ib), intent(in)  :: stringK(n_intK)    

    ! MO occupations
    integer(is), intent(in)  :: occB(nelX),occK(nelX)

    ! Work arrays
    real(dp), intent(out)    :: work(nelX,nelX)
    integer(is), intent(out) :: ipiv(nelX)
    
    ! Alpha/beta factor
    real(dp), intent(out)    :: fac

    ! Everything else
    integer(is)              :: i,j
    
!----------------------------------------------------------------------
! Fill in the matrix of occupied bra-ket MO overlaps
!----------------------------------------------------------------------    
    do i=1,nelX
       do j=1,nelX
          work(j,i)=smo(occB(j),occK(i))
       enddo
    enddo

!----------------------------------------------------------------------
! Screening based on Hadamard's inequality
!----------------------------------------------------------------------
    if (hadamard_bound(nelX,work) < hthrsh) then
       fac=0.0d0
       return
    endif
    
!----------------------------------------------------------------------    
! Determinant of the matrix of MO overlaps
!----------------------------------------------------------------------    
    fac=ludet(nelX,work,ipiv)

    return
    
  end subroutine get_one_factor
    
!######################################################################
! ludet: Calculation of the determinant of a square matrix via its
!        LU factorisation. The input matrix is overwritten.
!######################################################################
  function ludet(dim,mat,ipiv)

    use constants
    
    implicit none

    real(dp)                 :: ludet
    integer(is), intent(in)  :: dim
    real(dp), intent(inout)  :: mat(dim,dim)
    integer(is), intent(out) :: ipiv(dim)
    
    integer(is)              :: info,i

!----------------------------------------------------------------------
! LU factorisation of the input matrix
! Note that, on exit, mat will contain the L and U factors of the
! LU factorisation
!----------------------------------------------------------------------
    call dgetrf(dim,dim,mat,dim,ipiv,info)

    ! If info > 0, then the factor U is singular and the
    ! determinant is zero
    if (info > 0) then
       ludet=0.0d0
       return
    endif
    
    ! Exit if the call to dgetrf failed
    if (info < 0) then
       write(6,'(/,x,a,i3)') &
            'Error in ludet: LU decomposition failed, info=',info
       stop
    endif
    
!----------------------------------------------------------------------
! Determinant of the input matrix
!----------------------------------------------------------------------
    ludet=1.0d0
    do i=1,dim
       if (ipiv(i) /= i) then
          ludet=-ludet*mat(i,i)
       else
          ludet=ludet*mat(i,i)
       endif
    enddo
    
    return
    
  end function ludet

!######################################################################
! hadamard_bound: Calculation of the upper bound of the determinant
!                 of an input matrix using Hadamard's inequality
!######################################################################
  function hadamard_bound(dim,mat)

    use constants
    
    implicit none

    ! Function result
    real(dp)                :: hadamard_bound

    ! Input matrix
    integer(is), intent(in) :: dim
    real(dp), intent(in)    :: mat(dim,dim)

    ! Everything else
    integer(is)             :: i

    hadamard_bound=1.0d0

    do i=1,dim
       hadamard_bound=hadamard_bound*dot_product(mat(:,i),mat(:,i))
    enddo
    
    return
    
  end function hadamard_bound
    
!######################################################################
  
end module factors
