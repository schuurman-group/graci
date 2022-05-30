!**********************************************************************
! Routines for the calculation of expectation values in the determinant
! basis
!**********************************************************************
module det_expec

  implicit none

contains

!#######################################################################
! s2expec_detbas: Computes the <S^2> expectation values for a set of
!                 MRCI wave functions expanded in the determinant basis
!#######################################################################
  subroutine s2expec_detbas(detdim,nroots,det,vec,s2)

    use constants
    use bitglobal
    use bitutils
    use slater_condon, only: exc_degree_det,exc,s2ii,s2ij,&
         phasemask,phase_pure_exc
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: detdim,nroots

    ! Determinant bit strings
    integer(ib), intent(in) :: det(n_int,2,detdim)

    ! Eigenvectors in the determinant basis
    real(dp), intent(in)    :: vec(detdim,nroots)

    ! S^2 expectation values
    real(dp), intent(out)   :: s2(nroots)
    
    ! S^2 build variables
    integer(ib)              :: phase_mask(n_int,2)
    integer(is)              :: nexci
    integer(ib)              :: p(n_int,2),h(n_int,2)
    integer(is), parameter   :: maxex=2
    integer(is)              :: plist(maxex,2),hlist(maxex,2)
    integer(is)              :: phase

    ! S^2 - eigenvector product
    real(dp), allocatable    :: s2vec(:,:)
    
    ! Everything else
    integer(is)             :: idet,iket,ibra,ista
    real(dp)                :: val

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(s2vec(detdim,nroots))
    s2vec=0.0d0
    
!----------------------------------------------------------------------
! Contribution from the on-diagonal elements
!----------------------------------------------------------------------
    do idet=1,detdim
       s2vec(idet,:)=s2ii(det(:,:,idet))*vec(idet,:)
    enddo
    
!----------------------------------------------------------------------
! Off-diagonal elements
!----------------------------------------------------------------------
    ! Loop over ket determinants
    do iket=1,detdim-1

       ! Compute the phase mask
       call phasemask(det(:,:,iket),phase_mask)
       
       ! Loop over bra determinants
       do ibra=iket+1,detdim

          ! Cycle if the no. excitations connecting the bra and
          ! ket determinants is greater than 2
          nexci=exc_degree_det(det(:,:,ibra),det(:,:,iket))
          if (nexci > 2) cycle

          ! Get the indices of the spin-orbitals involved in the
          ! excitations linking the bra and ket determinants
          call exc(det(:,:,iket),det(:,:,ibra),p,h)
          call list_from_bitstring(p(:,ialpha),plist(:,ialpha),maxex)
          call list_from_bitstring(h(:,ialpha),hlist(:,ialpha),maxex)
          call list_from_bitstring(p(:,ibeta),plist(:,ibeta),maxex)
          call list_from_bitstring(h(:,ibeta),hlist(:,ibeta),maxex)

          ! Compute the phase factor
          phase=phase_pure_exc(phase_mask(:,ialpha),hlist(:,ialpha), &
               plist(:,ialpha),maxex,nexci) &
               *phase_pure_exc(phase_mask(:,ibeta),hlist(:,ibeta), &
               plist(:,ibeta),maxex,nexci)

          ! Compute the S^2 matrix element
          val=s2ij(phase,hlist,plist)
          
          ! Contributions to the S^2 expectation values
          s2vec(iket,:)=s2vec(iket,:)+val*vec(ibra,:)
          s2vec(ibra,:)=s2vec(ibra,:)+val*vec(iket,:)
             
       enddo

    enddo
    
!----------------------------------------------------------------------
! Compute the <S^2> values
!----------------------------------------------------------------------
    do ista=1,nroots
       s2(ista)=dot_product(vec(:,ista),s2vec(:,ista))
    enddo
    
    return
    
  end subroutine s2expec_detbas
  
!#######################################################################
  
end module det_expec
