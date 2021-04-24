!**********************************************************************
! Routines for the calculation of 1-RDMs for MRCI wavefunctions
!**********************************************************************
module rdm

  implicit none

contains

!######################################################################
! 1rdm_mrci: Master routine for the calculation of MRCI 1-RDMs for
!            a specified set of states of a single irrep
!######################################################################
  subroutine rdm_mrci(cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    use iomod
    use timing
    
    implicit none

    ! No. CSFs
    integer(is), intent(in) :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in) :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(out)   :: dmat(nmo,nmo,nroots)
   
    ! Timing variables
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                                twall_end 
    
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)  
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    dmat=0.0d0

!----------------------------------------------------------------------
! (1) On-diagonal elements
!----------------------------------------------------------------------
    call rdm_mrci_diag(cfg,csfdim,nroots,vec,dmat)
    
!----------------------------------------------------------------------
! (2) Off-diagonal elements of the Ref-Ref block
!----------------------------------------------------------------------
    call rdm_mrci_0h_0h(cfg,csfdim,nroots,vec,dmat)
    
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'rdm_mrci')
    
    return
    
  end subroutine rdm_mrci

!######################################################################
! rdm_mrci_ondiag: Calculation of the on-diagonal elements of the MRCI
!                  1-RDMs
!######################################################################
  subroutine rdm_mrci_diag(cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! No. CSFs
    integer(is), intent(in) :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in) :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)

!----------------------------------------------------------------------
! Ref space configurations
!----------------------------------------------------------------------
    
    
    return
    
  end subroutine rdm_mrci_diag
    
!######################################################################
! rdm_mrci_0h_0h: Calculation of the off-diagonal elements of the
!                 Ref-Ref block of the MRCI 1-RDMs
!######################################################################  
  subroutine rdm_mrci_0h_0h(cfg,csfdim,nroots,vec,dmat)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! No. CSFs
    integer(is), intent(in) :: csfdim
    
    ! No. roots for which the RDMs will be computed
    integer(is), intent(in) :: nroots
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Eigenvectors
    real(dp), intent(in)    :: vec(csfdim,nroots)
    
    ! Density matrix
    real(dp), intent(inout) :: dmat(nmo,nmo,nroots)

    print*,'>Here?'
    stop
    
    return
    
  end subroutine rdm_mrci_0h_0h
    
!######################################################################
  
end module rdm
