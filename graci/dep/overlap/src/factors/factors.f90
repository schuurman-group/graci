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
    use timing
    
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

    ! Timing variables
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                                twall_end
    
    ! Everything else
    integer(is)              :: ibra,iket,m,n

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)
    
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

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    if (verbose) &
         call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'get_all_factors')

    return
    
  end subroutine get_all_factors

!######################################################################
! get_all_factors_schur: Calculation of the complete set of alpha/beta
!                        factors given sets of unique bra and ket
!                        alpha/beta strings using Schur's determinant
!                        identity
!######################################################################
  subroutine get_all_factors_schur(nfixedX,nvarX,nstringB,nstringK,&
       stringB,stringK,fac)

    use constants
    use global, only: n_intB,n_intK,nmoB,nmoK,smo,hthrsh,verbose,&
                      invSff,detSff
    use detfuncs
    use timing
    
    implicit none

    ! Number of fixed- and variable-occupation orbitals in each
    ! bra and ket string
    ! (these have to be equal for an overlap calculation)
    integer(is), intent(in)  :: nfixedX,nvarX
    
    ! Dimensions
    integer(is), intent(in)  :: nstringB,nstringK

    ! Unique alpha/beta strings
    integer(ib), intent(in)  :: stringB(n_intB,nstringB)
    integer(ib), intent(in)  :: stringK(n_intK,nstringK)    

    ! Alpha/beta factors
    real(dp), intent(out)    :: fac(nstringB,nstringK)

    ! Occupied MOs
    integer(is)              :: nelX,noccB,noccK
    integer(is), allocatable :: occB(:),occK(:)

    ! Work arrays
    real(dp), allocatable    :: S(:,:),Svv(:,:),Svf(:,:),Sfv(:,:)
    real(dp), allocatable    :: work(:,:),invSffSfv(:,:)
    integer(is), allocatable :: ipiv(:)

    ! Timing variables
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                                twall_end
    
    ! Everything else
    integer(is)              :: ibra,iket,m,n,mm,nn

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Total number of electrons
    nelX=nfixedX+nvarX

    ! Orbital occupations
    allocate(occB(nelX), occK(nelX))
    occB=0; occK=0

    ! Orbital overlap matrix
    allocate(S(nelX,nelX))
    S=0.0d0
    
    ! Blocks of the orbital overlap matrix
    allocate(Svv(nvarX,nvarX), Svf(nvarX,nfixedX), Sfv(nfixedX,nvarX))
    Svv=0.0d0; Svf=0.0d0; Sfv=0.0d0
    
    ! Work arrays
    allocate(ipiv(nvarX), work(nvarX,nvarX), invSffSfv(nfixedX,nvarX))
    ipiv=0; work=0.0d0; invSffSfv=0.0d0

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

          ! Fill in the complete bra-ket MO overlap matrix
          do m=1,nelX
             do n=1,nelX
                S(n,m)=smo(occB(n),occK(m))
             enddo
          enddo

          ! Determinant screening
          if (hadamard_bound(nelX,S) < hthrsh) then
             fac(ibra,iket)=0.0d0
             cycle
          endif
          
          ! Fill in the var-var block of the bra-ket MO overlap matrix,
          ! S^(v,v)
          do m=1,nvarX
             mm=m+nfixedX
             do n=1,nvarX
                nn=n+nfixedX
                Svv(n,m)=S(nn,mm)
             enddo
          enddo

          ! Fill in the var-fixed block of the bra-ket MO overlap matrix,
          ! S^(v,f)
          do m=1,nfixedX
             do n=1,nvarX
                nn=n+nfixedX
                Svf(n,m)=S(nn,m)
             enddo
          enddo

          ! Fill in the fixed-var block of the bra-ket MO overlap matrix,
          ! S^(f,v)
          do m=1,nvarX
             mm=m+nfixedX
             do n=1,nfixedX
                Sfv(n,m)=S(n,mm)
             enddo
          enddo

          ! [S^(f,f)]^-1 S^(f,v)
          call dgemm('N','N',nfixedX,nvarX,nfixedX,1.0d0,invSff,&
               nfixedX,Sfv,nfixedX,0.0d0,invSffSfv,nfixedX)
          
          ! S^(v,v) - S^(v,f) [S^(f,f)]^-1 S^(f,v)
          work=Svv
          call dgemm('N','N',nvarX,nvarX,nfixedX,-1.0d0,Svf,nvarX,&
               invSffSfv,nfixedX,1.0d0,work,nvarX)

          ! Determinant of the matrix of MO overlaps
          fac(ibra,iket)=detSff*ludet(nvarX,work,ipiv)
          
       enddo

    enddo

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    if (verbose) &
         call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'get_all_factors_schur')
    
    return
    
  end subroutine get_all_factors_schur
    
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
    
    ! Alpha/beta strings
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
! get_one_factor_schur: Calculation of a single alpha/beta factor given
!                       a pair of alpha/beta strings and corresponding
!                       MO occupations using Schur's determinant
!                       identity
!######################################################################
  subroutine get_one_factor_schur(nelX,nfixedX,nvarX,stringB,stringK,&
       occB,occK,S,Svv,Svf,Sfv,invSffSfv,work,ipiv,fac)

    use constants
    use global, only: n_intB,n_intK,nmoB,nmoK,smo,hthrsh,verbose,&
                      invSff,detSff
    use detfuncs
    
    implicit none

    ! Number of electrons
    integer(is), intent(in) :: nelX
    
    ! Number of fixed- and variable-occupation orbitals in each
    ! bra and ket string
    ! (these have to be equal for an overlap calculation)
    integer(is), intent(in)  :: nfixedX,nvarX
    
    ! Alpha/beta strings
    integer(ib), intent(in)  :: stringB(n_intB)
    integer(ib), intent(in)  :: stringK(n_intK)    

    ! MO occupations
    integer(is), intent(in)  :: occB(nelX),occK(nelX)
    
    ! Work arrays
    real(dp), intent(out)    :: S(nelX,nelX),Svv(nvarX,nvarX),&
                                Svf(nvarX,nfixedX),&
                                Sfv(nfixedX,nvarX),&
                                invSffSfv(nfixedX,nvarX)
    real(dp), intent(out)    :: work(nvarX,nvarX)
    integer(is), intent(out) :: ipiv(nvarX)
    
    ! Alpha/beta factor
    real(dp), intent(out)    :: fac

    ! Everything else
    integer(is)              :: m,n,mm,nn

!----------------------------------------------------------------------
! Fill in the complete bra-ket MO overlap matrix
!----------------------------------------------------------------------
    do m=1,nelX
       do n=1,nelX
          S(n,m)=smo(occB(n),occK(m))
       enddo
    enddo

!----------------------------------------------------------------------
! Determinant screening
!----------------------------------------------------------------------
    if (hadamard_bound(nelX,S) < hthrsh) then
       fac=0.0d0
       return
    endif    
    
!----------------------------------------------------------------------
! Fill in the var-var block of the bra-ket MO overlap matrix, S^(v,v)
!----------------------------------------------------------------------
    do m=1,nvarX
       mm=m+nfixedX
       do n=1,nvarX
          nn=n+nfixedX
          Svv(n,m)=S(nn,mm)
       enddo
    enddo

!----------------------------------------------------------------------
! Fill in the var-fixed block of the bra-ket MO overlap matrix, S^(v,f)
!----------------------------------------------------------------------
    do m=1,nfixedX
       do n=1,nvarX
          nn=n+nfixedX
          Svf(n,m)=S(nn,m)
       enddo
    enddo

!----------------------------------------------------------------------
! Fill in the fixed-var block of the bra-ket MO overlap matrix, S^(f,v)
!----------------------------------------------------------------------
    do m=1,nvarX
       mm=m+nfixedX
       do n=1,nfixedX
          Sfv(n,m)=S(n,mm)
       enddo
    enddo

!----------------------------------------------------------------------
! [S^(f,f)]^-1 S^(f,v)
!----------------------------------------------------------------------
    call dgemm('N','N',nfixedX,nvarX,nfixedX,1.0d0,invSff,nfixedX,Sfv,&
         nfixedX,0.0d0,invSffSfv,nfixedX)
    
!----------------------------------------------------------------------
! S^(v,v) - S^(v,f) [S^(f,f)]^-1 S^(f,v)
!----------------------------------------------------------------------
    work=Svv
    call dgemm('N','N',nvarX,nvarX,nfixedX,-1.0d0,Svf,nvarX,invSffSfv,&
         nfixedX,1.0d0,work,nvarX)
    
!----------------------------------------------------------------------
! Determinant of the matrix of MO overlaps
!----------------------------------------------------------------------
    fac=detSff*ludet(nvarX,work,ipiv)
    
    return
    
  end subroutine get_one_factor_schur
    
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
  
end module factors
