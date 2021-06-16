!**********************************************************************
! Routines for the calculation of the determinant expansions of MRCI
! wave functions
!**********************************************************************
module csf2det

  implicit none

contains
  
!######################################################################
! csf2det: master routine controling the calculation and dumping to
!          file of a set of MRCI wave functions in the determinant
!          basis
!######################################################################
  subroutine csf2det_mrci(cfg,nroots,csfdim,vec,ndet)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! Dimensions
    integer(is), intent(in)  :: nroots,csfdim

    ! Eigenvectors
    real(dp), intent(in)     :: vec(csfdim,nroots)

    ! No. determinants
    integer(is), intent(out) :: ndet
    
    print*,'here?'
    stop
    
    return
    
  end subroutine csf2det_mrci

!######################################################################
  
end module csf2det
