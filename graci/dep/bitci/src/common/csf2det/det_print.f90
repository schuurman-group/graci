module det_print

  implicit none

contains
  
!######################################################################
! print_eigvec_det: Outputs a set of eigenvectors in the determinant
!                   representation. For the sake of brevity, the
!                   dominant determinants are output as their
!                   difference wrt the base determinant.
!######################################################################
  subroutine print_eigvec_det(ndet,nroots,lddet,vec,det)
    
    use constants
    use bitglobal
    use utils
    
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: ndet,nroots,lddet

    ! Eigenvectors
    real(dp), intent(in)     :: vec(ndet,nroots)

    ! Determinant bit strings
    integer(ib), intent(in)  :: det(lddet,2,ndet)

    ! Sorting
    integer(is), allocatable :: indx(:)
    real(dp), allocatable    :: absvec(:)
    
    ! Everything else
    integer(is)              :: i,j,k
    real(dp), parameter      :: pthrsh=0.03d0
    real(dp)                 :: coe
    character(len=60)        :: string
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(indx(ndet))
    allocate(absvec(ndet))
    indx=0
    absvec=0.0d0

!----------------------------------------------------------------------
! Output the dominant determinants and coefficients
!----------------------------------------------------------------------
    ! Loop over roots
    do i=1,nroots

       ! Sort the eigenvector by absolute coefficient value
       absvec=abs(vec(:,i))
       call dsortindxa1('D',ndet,absvec,indx)

       ! Header
       write(6,'(/,x,50a)') ('-',k=1,50)
       write(6,'(x,a,x,i0)') 'State',i
       write(6,'(x,50a)') ('-',k=1,50)
       
       ! Output the dominant determinants
       do j=1,ndet

          ! Coefficient
          coe=vec(indx(j),i)

          ! Exit if we have reached the end of the dominant
          ! coefficients
          if (abs(coe) < pthrsh) exit

          ! Get the determinant difference character string
          call get_diff_string(lddet,det(:,:,indx(j)),string)

          ! Output the coefficient, difference string pair
          write(6,'(x,F10.7,x,a,x,a)') coe,'|',trim(string)
          
       enddo

       ! Footer
       write(6,'(x,50a)') ('-',k=1,50)
       
    enddo
    
    return
    
  end subroutine print_eigvec_det

!######################################################################
! diff_string: Given a determinant bit string, det, returns a
!              character string corresponding to the excitation wrt
!              the base determinant
!######################################################################
  subroutine get_diff_string(lddet,det_in,string)

    use constants
    use bitglobal
    use slater_condon
    use bitutils
    
    implicit none

    ! Input determinant bit string
    integer(is), intent(in)        :: lddet
    integer(ib), intent(in)        :: det_in(lddet,2)

    ! Output character string
    character(len=60), intent(out) :: string

    ! Particle and hole indices
    integer(ib)                    :: p(n_int,2),h(n_int,2)
    integer(is)                    :: plist(nexmax,2),hlist(nexmax,2)
    integer(is)                    :: nexci(2)

    ! Everything else
    integer(is)                    :: i,i1,dim
    integer(ib)                    :: det(n_int,2)

!----------------------------------------------------------------------
! Working determinant bit string
!
! This is necessary because we may have been passed a set of
! determinants from a previous geometry in, e.g., a root-following
! calculation
!
! In such a case, it is possible that the leading dimension does not
! equal n_int, but possibly n_int+-1 due to a variation of the no.
! MOs within the (truncated) virtual space
!
! For now, we will just assume that the first n_int blocks are all
! that is needed, as should be the case for the dominant determinants
!
! This is _slightly_ risky, but it avoids the alteration of the
! Slater-Condon and bitutils routines, which assume bit strings of
! leading dimension n_int
!----------------------------------------------------------------------
    if (n_int == lddet) then
       det=det_in
    else if (n_int < lddet) then
       det=det_in(1:n_int,:)
    else if (n_int > lddet) then
       det=0_ib
       det(1:lddet,:)=det_in
    endif
    
!----------------------------------------------------------------------
! Get the bit string encodings of the particle and hole MOs relative
! base determinant
!----------------------------------------------------------------------
    call exc(det0,det,p,h)

!----------------------------------------------------------------------
! Get the alpha and beta excitation degrees
!----------------------------------------------------------------------
    nexci=0
    nexci(1)=exc_degree_string(det0(:,1),det(:,1))
    nexci(2)=exc_degree_string(det0(:,2),det(:,2))
    
!----------------------------------------------------------------------
! Get the lists of alpha and beta particle/hole MO indices
!----------------------------------------------------------------------
    ! alpha particles
    call list_from_bitstring(p(:,1),plist(:,1),nexmax)

    ! alpha holes
    call list_from_bitstring(h(:,1),hlist(:,1),nexmax)

    ! beta particles
    call list_from_bitstring(p(:,2),plist(:,2),nexmax)

    ! beta holes
    call list_from_bitstring(h(:,2),hlist(:,2),nexmax)

!----------------------------------------------------------------------
! Write the excitation character string
!----------------------------------------------------------------------
    ! Is this the base determinant?
    if (sum(nexci) == 0) then
       string=' base det'
       return
    endif

    ! Initialisation
    string=''

    ! Aplha holes
    do i=1,nexci(1)
       i1=len_trim(string)
       write(string(i1+1:),'(x,i0,a1)') hlist(i,1),'a'
    enddo

    ! Beta holes
    do i=1,nexci(2)
       i1=len_trim(string)
       write(string(i1+1:),'(x,i0,a1)') hlist(i,2),'b'
    enddo

    ! Hole-particle delimiter
    i1=len_trim(string)
    write(string(i1+1:),'(x,a)') '->'

    ! Aplha particles
    do i=1,nexci(1)
       i1=len_trim(string)
       write(string(i1+1:),'(x,i0,a1)') plist(i,1),'a'
    enddo

    ! Beta particles
    do i=1,nexci(2)
       i1=len_trim(string)
       write(string(i1+1:),'(x,i0,a1)') plist(i,2),'b'
    enddo
    
    return
    
  end subroutine get_diff_string
  
!######################################################################
  
end module det_print
