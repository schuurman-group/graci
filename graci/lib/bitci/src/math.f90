module math

  implicit none

contains

!######################################################################
! bico: numerical recipes binomial coefficient function
!       *note that this returns the binomial coefficient as a
!       floating-point number*
!######################################################################

  function bico(n,k)

    use constants
    
    implicit none

    integer(is), intent(in) :: k,n
    real(dp)                :: bico

    bico=nint(exp(factln(n)-factln(k)-factln(n-k)))
    
    return
    
  end function bico

!######################################################################
! factln: numerical recipes function to compute ln(n!)
!######################################################################

  function factln(n)

    use constants

    implicit none

    integer(is), intent(in) :: n
    real(dp)                :: factln
    real(dp)                :: a(100)

    save a
    data a/100*-1./
    
    if (n.lt.0) then
       write(6,'(/,a,/)') 'negative factorial in factln'
       stop
    endif
       
    if (n.le.99) then
       if (a(n+1).lt.0._dp) a(n+1)=gammln(n+1._dp)
       factln=a(n+1)
    else
       factln=gammln(n+1._dp)
    endif
    
    return
    
  end function factln

!######################################################################
! gammln: numerical recipes routine to calculate ln(Gamma(xx)) for
!         xx > 0
!######################################################################
  function gammln(xx)

    use constants
    
    implicit none

    real(dp), intent(in) :: xx
    real(dp)             :: gammln
    integer(is)          :: i,j
    real(dp)             :: ser,stp,tmp,x,y,cof(6)

    save cof,stp

    data cof,stp/76.18009172947146d0,-86.50532032941677d0,&
         24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
         -.5395239384953d-5,2.5066282746310005d0/

    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
       y=y+1.d0
       ser=ser+cof(j)/y
    enddo
    gammln=tmp+log(stp*ser/x)
    
    return
    
  end function gammln
    
!######################################################################
  
end module math
