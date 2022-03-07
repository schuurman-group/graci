module clebsch_gordan

  implicit none

contains

!######################################################################
! cgcoeff_eigen: For a given pair of subsystem angular momentum quantum
!                numbers j1 and j2, computes the Clebsh-Gordan
!                coefficients for total angular momentum J and all
!                possible (m1,m2,M) tuples
!                
!                Phases are fixed according the Condon-Shortley
!                convention
!
!                The output array cg is assumed to be dimensioned as
!                ((2*j1+1)*(2*j2+1), 2*J+1), and the direct product
!                basis |j1m1> x |j2m2> is assumed to be traversed in
!                lexicographical order
!
!                Follows the algorithm given in the Appendix of
!                J. Chem. Phys., 151, 034106 (2019)
!######################################################################
  subroutine cgcoeff(j1,j2,J,cg)

    use constants
    use iomod
    
    implicit none

    ! Angular momentum quantum numbers
    real(dp), intent(in)     :: j1,j2,J

    ! Output array of Clebsch-Gordan coefficients
    real(dp)                 :: cg(:,:)
    
    ! Direct product basis representation of J^2, J_z and J^-, etc
    integer(is)              :: dim1,dim2,dim
    real(dp), allocatable    :: Jsq(:,:),Jz(:,:),Jminus(:,:)
    real(dp), allocatable    :: j1z(:,:),j2z(:,:)
    real(dp), allocatable    :: mval(:,:)

    ! Number of eigenvectors of J^2 with eigenvalues J(J+1)
    integer(is)              :: dimJ

    ! Projector onto the subspace spanned by the eigenvectors
    ! of J^2 with eigenvalues J(J+1)
    real(dp), allocatable    :: P(:,:)

    ! Projected J_z matrix
    real(dp), allocatable    :: PJzP(:,:)
    real(dp), allocatable    :: Jzvec(:,:)
    
    ! Work arrays
    real(dp), allocatable    :: eigvec(:,:),eigval(:)
    real(dp), allocatable    :: work(:)
        
    ! Evertything else
    integer(is)              :: ibra,iket,i,k,l,i1,i2,n,info
    integer(is), allocatable :: indx(:)
    real(dp)                 :: M,m1val,m2val
    real(dp)                 :: fac
    real(dp), parameter      :: shift=100.0d0
    real(dp), parameter      :: thrsh=1e-10_dp
    
!----------------------------------------------------------------------
! Check on the input angular momentum quantum numbers
!----------------------------------------------------------------------
    if (j1 < 0.0d0 .or. mod(j1,0.5d0) /= 0.0d0) then
       errmsg='Error in cgcoeff: illegal input value of j1'
       call error_control
    endif

    if (j2 < 0.0d0 .or. mod(j2,0.5d0) /= 0.0d0) then
       errmsg='Error in cgcoeff: illegal input value of j2'
       call error_control
    endif

    if (J < 0.0d0 .or. mod(J,0.5d0)  /= 0.0d0) then
       errmsg='Error in cgcoeff: illegal input value of J'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Size of the direct product basis |j1,m1> x |j2,m2>
!----------------------------------------------------------------------
    dim1=int(2*j1)+1
    dim2=int(2*j2)+1
    dim=dim1*dim2
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(Jsq(dim,dim))
    Jsq=0.0d0

    allocate(Jz(dim,dim))
    Jz=0.0d0

    allocate(j1z(dim,dim))
    j1z=0.0d0

    allocate(j2z(dim,dim))
    j2z=0.0d0
    
    allocate(Jminus(dim,dim))
    Jminus=0.0d0
    
    allocate(mval(2,dim))
    mval=0.0d0

    allocate(eigvec(dim,dim))
    eigvec=0.0d0
    
    allocate(eigval(dim))
    eigval=0.0d0

    allocate(work(3*dim))
    work=0.0d0

    dimJ=int(2*J)+1
    allocate(indx(dimJ))
    indx=0

    allocate(P(dim,dim))
    P=0.0d0

    allocate(PJzP(dim,dim))
    PJzP=0.0d0

    allocate(Jzvec(dim,dimJ))
    Jzvec=0.0d0
    
!----------------------------------------------------------------------
! Get the values of (m1,m2) for each element of the direct product
! basis
!----------------------------------------------------------------------
    ! Loop over direct product basis functions |j1,m1> x |j2,m2>
    n=0
    do i1=1,dim1
       m1val=-j1-1+i1
       do i2=1,dim2
          m2val=-j2-1+i2
          n=n+1

          ! (m1,m2) value
          mval(1,n)=m1val
          mval(2,n)=m2val
          
       enddo
    enddo

!----------------------------------------------------------------------
! Construct the direct product basis representation of J^2
!----------------------------------------------------------------------
    ! Loop over kets
    do iket=1,dim

       ! Loop over bras
       do ibra=iket,dim

          ! Diagonal elements: contribution from
          ! j1^2 + j2^2 + 2 j1_z j2_z
          if (iket == ibra) then
             Jsq(ibra,iket)=j1*(j1+1.0d0)+j2*(j2+1.0d0) &
                  +2.0d0*mval(1,iket)*mval(2,iket)
          endif
          
          ! Off-diagonal elements: contribution from
          ! j1^+ j2^- + 2 j1^- j2^+
          if (iket /= ibra) then
             
             ! Contribution (i): j1^+ j2^-
             if (abs(mval(1,iket)+1.0d0) <= j1 .and. &
                  abs(mval(2,iket)-1.0d0) <= j2) then

                if (mval(1,ibra) == mval(1,iket)+1 .and. &
                     mval(2,ibra) == mval(2,iket)-1) then

                   Jsq(ibra,iket)= &
                        sqrt(j1*(j1+1.0d0)-mval(1,iket)*(mval(1,iket)+1)) &
                        *sqrt(j2*(j2+1.0d0)-mval(2,iket)*(mval(2,iket)-1))
                   
                endif
                                   
             endif

             ! Contribution (ii): j1^- j2^+
             if (abs(mval(1,iket)-1.0d0) <= j1 .and. &
                  abs(mval(2,iket)+1.0d0) <= j2) then

                if (mval(1,ibra) == mval(1,iket)-1 .and. &
                     mval(2,ibra) == mval(2,iket)+1) then

                   Jsq(ibra,iket)=Jsq(ibra,iket)+ &
                        sqrt(j1*(j1+1.0d0)-mval(1,iket)*(mval(1,iket)-1)) &
                        *sqrt(j2*(j2+1.0d0)-mval(2,iket)*(mval(2,iket)+1))
                   
                endif
                   
             endif

             ! Fill in the upper triangle
             Jsq(iket,ibra)=Jsq(ibra,iket)
             
          endif
          
       enddo

    enddo

!----------------------------------------------------------------------
! Construct the direct product representations of J_z, j1_z and j2_z
! Note that the j1_z and j2_z matrices will be used to check that the
! generated Clebsch-Gordan coefficients satisfy the Condon-Shortley
! phase convention
!----------------------------------------------------------------------
    ! Loop over kets
    do i=1,dim

       ! < j1 m1 j2 m2 | J_z | j1 m1 j2 m2 >
       Jz(i,i)=mval(1,i)+mval(2,i)

       ! < j1 m1 j2 m2 | j1_z | j1 m1 j2 m2 >
       j1z(i,i)=mval(1,i)

       ! < j1 m1 j2 m2 | j2_z | j1 m1 j2 m2 >
       j2z(i,i)=mval(1,i)
       
    enddo

!----------------------------------------------------------------------
! Construct the direct product representation of J^- = j1^- + j2^-
!----------------------------------------------------------------------
    ! Loop over kets
    do iket=1,dim
       
       ! Loop over bras
       do ibra=1,dim

          ! (a) < j1 m1 j2 m2 | j1^- | j1 m1' j2 m2' >
          if (abs(mval(1,iket)-1.0d0) <= j1 .and. &
               mval(1,ibra) == mval(1,iket)-1) then
             Jminus(ibra,iket)=&
                  sqrt(j1*(j1+1.0d0)-mval(1,iket)*(mval(1,iket)-1))
          endif
             
          ! (b) < j1 m1 j2 m2 | j2^- | j1 m1' j2 m2' >
          if (abs(mval(2,iket)-1.0d0) <= j2 .and. &
               mval(2,ibra) == mval(2,iket)-1) then
             Jminus(ibra,iket)=Jminus(ibra,iket) &
                  +sqrt(j2*(j2+1.0d0)-mval(2,iket)*(mval(2,iket)-1))
          endif
          
       enddo

    enddo
       
!----------------------------------------------------------------------
! Diagonalise the direct product representation of J^2
!----------------------------------------------------------------------
    eigvec=Jsq
    call dsyev('V','U',dim,eigvec,dim,eigval,work,3*dim,info)

    if (info /= 0) then
       errmsg='Error in cgcoeff: diagonalisation of J^2 failed'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Get the indices of the eigenvectors with eigenvalues J(J+1)
!----------------------------------------------------------------------
    n=0
    do i=1,dim
       if (abs(eigval(i)-J*(J+1)) < thrsh) then
          n=n+1
          indx(n)=i
       endif
    enddo

    ! Sanity check
    if (n /= dimJ) then
       errmsg='Error in cgcoeff: incorrect eigenvalues of J^2 found'
       call error_control
    endif

!----------------------------------------------------------------------
! Form the projector P onto the subspace spanned by the eigenvectors
! of J^2 with eigenvalues J(J+1)
!----------------------------------------------------------------------
    do i=1,dimJ       
       do k=1,dim
          do l=1,dim
             P(l,k)=P(l,k)+eigvec(l,indx(i))*eigvec(k,indx(i))
          enddo
       enddo
    enddo
        
!----------------------------------------------------------------------
! Project form the direct product representation of
! P (J_z + shift*1) P
! Here, the shift is used to enable us to distinguish the M=0
! eigenvectors from those spanning the null space
!----------------------------------------------------------------------
    do i=1,dim
       Jz(i,i)=Jz(i,i)+shift
    enddo

    PJzP=matmul(P,matmul(Jz,P))
    
!----------------------------------------------------------------------
! Diagonalise the projected J_z matrix
!----------------------------------------------------------------------
    eigvec=PJzP

    call dsyev('V','U',dim,eigvec,dim,eigval,work,3*dim,info)

    if (info /= 0) then
       errmsg='Error in cgcoeff: diagonalisation of P J_z P failed'
       call error_control
    endif

!----------------------------------------------------------------------
! Fill in the eigenvectors P(Jz + shift*1)P orthogonal to the null
! space
!----------------------------------------------------------------------
    n=0
    do i=1,dim
       if (abs(eigval(i)) > thrsh) then
          n=n+1
          Jzvec(:,n)=eigvec(:,i)
       endif
    enddo

    if (n /= dimJ) then
       errmsg='Error in cgcoeff: incorrect eigenvalues of ' &
            //'P(Jz+shift*1)P found'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Fix the phases of the Clebsch-Gordan coefficients s.t. they obey
! the Condon-Shortley convention
!----------------------------------------------------------------------
    !
    ! Anchor point: ensure that <j1 m1=j1; j2 m2=(J-j1) | J M=J>
    ! is positive
    !
    ! Index of the coefficient corresponding to
    ! <j1 m1=j1; j2 m2=(J-j1) | J M=J>
    do i=1,dim
       if (abs(mval(1,i)-j1) < thrsh &
            .and. abs(mval(2,i)-(J-j1)) < thrsh) then
          n=i
          exit
       endif
    enddo

    ! Fix the phase of the eigenvector for |J M=J>
    fac=sign(1.0d0,Jzvec(n,dimJ))
    Jzvec(:,dimJ)=fac*Jzvec(:,dimJ)

    !
    ! Remaining eigenvectors
    !
    M=J+1
    do i=dimJ,2,-1
       M=M-1

       ! Compute < J M-1 | J^- | J M >
       fac=dot_product(Jzvec(:,i-1),matmul(Jminus,Jzvec(:,i)))

       ! Sanity check
       if (abs(abs(fac)-sqrt(J*(J+1)-M*(M-1))) > thrsh) then
          errmsg='Error in the calculation of  < J M-1 | J^- | J M >'
          call error_control
       endif

       ! Sign of < J M-1 | J^- | J M >
       fac=sign(1.0d0,fac)

       ! Fix the phase of the current eigenvector
       Jzvec(:,i-1)=fac*Jzvec(:,i-1)
       
    enddo

!----------------------------------------------------------------------
! Package up the Clebsch-Gordan coefficients
!----------------------------------------------------------------------
    cg=Jzvec
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(Jsq)
    deallocate(Jz)
    deallocate(j1z)
    deallocate(j2z)
    deallocate(Jminus)
    deallocate(mval)
    deallocate(eigvec)
    deallocate(eigval)
    deallocate(work)
    deallocate(indx)
    deallocate(P)
    deallocate(PJzP)
    deallocate(Jzvec)
        
    return

  end subroutine cgcoeff

!######################################################################
  
end module clebsch_gordan
