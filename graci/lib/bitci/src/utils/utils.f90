module utils

  implicit none
  
contains

!######################################################################
! dsortindxa1: heap sort for 8-byte real arrays
!######################################################################
  subroutine dsortindxa1(order,ndim,arrin,indx)
  
    use constants
    implicit none

    character(1), intent(in) :: order
    integer, intent(in) :: ndim
    real(dp), dimension(ndim), intent(in)   :: arrin
    integer, dimension(ndim), intent(inout) :: indx
    
    integer :: i,l,ir,indxt,j
    real(dp) :: q
!!$ The subroutine is taken from the NR p233, employs heapsort.

    if (ndim.eq.1) then
       indx(1)=1
       return
    endif

       
    do i= 1,ndim
       indx(i)=i
    end do
    
    l=ndim/2+1
    ir=ndim
    
    if(order .eq. 'D') then
       
10     continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
20     if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .gt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .gt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 20
       end if
       indx(i)=indxt
       go to 10
       
    elseif(order .eq. 'A') then
       
100    continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
200    if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .lt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .lt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 200
       end if
       indx(i)=indxt
       go to 100
       
    end if
       
  end subroutine dsortindxa1
  
!######################################################################
! i4sortindx: heap sort for 4-byte integer arrays
!######################################################################
  subroutine i4sortindx(order,ndim,arrin,indx)
  
    use constants
  
    implicit none
  
    character(1), intent(in)                    :: order
    integer(is), intent(in)                     :: ndim
    integer(is), dimension(ndim), intent(in)    :: arrin
    integer(is), dimension(ndim), intent(inout) :: indx
    integer(is)                                 :: i,l,ir,indxt,j
    integer(is)                                 :: q
  
!!$ The subroutine is adapted from the NR p233, employs heapsort.
    
    if (ndim.eq.1) then
       indx(1)=1
       return
    endif
       
    do i= 1,ndim
       indx(i)=i
    end do
    
    l=ndim/2+1
    ir=ndim
    
    if(order .eq. 'D') then
       
10     continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
20     if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .gt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .gt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 20
       end if
       indx(i)=indxt
       go to 10
       
    elseif(order .eq. 'A') then
       
100    continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
200    if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .lt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .lt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 200
       end if
       indx(i)=indxt
       go to 100
       
    end if
    
  end subroutine i4sortindx

!######################################################################
! i8sortindx: heap sort for 8-byte integer arrays
!######################################################################
  subroutine i8sortindx(order,ndim,arrin,indx)
  
    use constants
  
    implicit none
  
    character(1), intent(in)                    :: order
    integer(is), intent(in)                     :: ndim
    integer(il), dimension(ndim), intent(in)    :: arrin
    integer(is), dimension(ndim), intent(inout) :: indx
    integer(il)                                 :: i,l,ir,indxt,j
    integer(il)                                 :: q
  
!!$ The subroutine is adapted from the NR p233, employs heapsort.
    
    if (ndim.eq.1) return

    do i= 1,ndim
       indx(i)=i
    end do
    
    l=ndim/2+1
    ir=ndim
    
    if(order .eq. 'D') then
       
10     continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
20     if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .gt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .gt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 20
       end if
       indx(i)=indxt
       go to 10
       
    elseif(order .eq. 'A') then
       
100    continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
200    if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .lt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .lt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 200
       end if
       indx(i)=indxt
       go to 100
       
    end if
    
  end subroutine i8sortindx
  
!######################################################################
! uppercase: replaces all lowercase letters in a given character
!            string with the corresponding uppercase letters
!######################################################################
  subroutine uppercase(string)

    use constants
    
    implicit none
      
    integer(is)  :: i
    character(*) :: string
      
    do i=1,len(string)
       if(string(i:i).ge."a".and.string(i:i).le."z")&
            string(i:i)=achar(iachar(string(i:i))-32)
    enddo
    
    return

  end subroutine uppercase

!######################################################################
! lowercase: replaces all uppercase letters in a given character
!            string with the corresponding lowercase letters
!######################################################################
  subroutine lowercase(string)

    use constants
    
    implicit none
      
    integer(is)  ::  i,j
    character(*) :: string
    
    do i=1,len(string)
       do j=65,90
          if (ichar(string(i:i)).eq.j) string(i:i)=char(j+32)
       enddo
    enddo
    
    return
    
  end subroutine lowercase

!######################################################################
! diag_matrix_real: computes the eigenpairs of a real, double
!                   precision matrix
!######################################################################
  subroutine diag_matrix_real(mat,eigval,eigvec,dim)

    use constants
    
    implicit none

    integer(is), intent(in)      :: dim
    integer(is)                  :: error
    real(dp), dimension(dim,dim) :: mat
    real(dp), dimension(dim)     :: eigval
    real(dp), dimension(dim,dim) :: eigvec
    real(dp), dimension(3*dim)   :: work

    eigvec=mat
    
    call dsyev('V','U',dim,eigvec,dim,eigval,work,3*dim,error)
    
    if (error.ne.0) then
       write(6,'(/,a,/)') 'Diagonalisation failed in subroutine '&
            //'diag_matrix_real'
       stop
    endif
    
    return
    
  end subroutine diag_matrix_real
  
!######################################################################
    
end module utils
