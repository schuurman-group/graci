!**********************************************************************
! Routines for the pruning of the reference space configurations based
! upon the reference space eigenvectors
!**********************************************************************
module ref_prune

  implicit none

contains

!######################################################################
! prune_ref_space: Pruning of the reference space based on the
!                  maximum coefficient values across all reference
!                  space eigenvectors
!######################################################################
  subroutine prune_ref_space(nconf,hdim,offset,conf,sop,n_int_I)

    use constants
    use bitglobal
    use iomod
    
    implicit none

    ! No. configurations & CSFs
    integer(is), intent(in) :: nconf,hdim

    ! CSF offsets
    integer(is), intent(in) :: offset(nconf+1)

    ! Reference configurations and SOPs
    integer(is), intent(in) :: n_int_I
    integer(ib), intent(in) :: conf(n_int_I,2,nconf)
    integer(ib), intent(in) :: sop(n_int_I,2,nconf)

    print*,"Here"
    stop
    
    return
    
  end subroutine prune_ref_space
    
!######################################################################
  
end module ref_prune
