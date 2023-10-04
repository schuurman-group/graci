!**********************************************************************
! Sundry determinant manipulation routines
!**********************************************************************
module detfuncs

  implicit none

contains

!######################################################################
! truncate_wave_functions: Truncation of the input wave functions
!                          based on a minimum norm threshold
!######################################################################
  subroutine truncate_wave_functions(n_int,ndet,nroots,det,vec,&
       normthrsh,ndet_new,det_new,vec_new)

    use constants
    use global
    use utils
    
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: n_int,ndet,nroots

    ! Determinant bit strings
    integer(ib), intent(in)  :: det(n_int,2,ndet)

    ! Eigenvectors
    real(dp), intent(in)     :: vec(ndet,nroots)

    ! Norm-based truncation threshold
    real(dp), intent(in)     :: normthrsh

    ! Truncated determinant and eigenvector arrays
    integer(is), intent(out) :: ndet_new
    integer(ib), allocatable :: det_new(:,:,:)
    real(dp), allocatable    :: vec_new(:,:)
    
    ! Sorting
    integer(is), allocatable :: indx(:)
    integer(is), allocatable :: iwork(:)

    ! Surviving determinants
    integer(is), allocatable :: idet(:)

    ! Everything else
    integer(is)              :: i,k,n,ntrim
    real(dp)                 :: normsq,targ,diff
    real(dp), parameter      :: epsilon=1e-8_dp
    real(dp)                 :: trim_thrsh
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(iwork(ndet))
    iwork=0
    
    allocate(indx(ndet))
    indx=0

    allocate(idet(ndet))
    idet=0
    
!----------------------------------------------------------------------
! Determine the indices of the surviving determinants
!----------------------------------------------------------------------
    ! Target squared norm
    targ=normthrsh**2

    ! Wave function trimming threshold
    trim_thrsh=sqrt((1.0d0-normthrsh)**2/ndet)
    
    ! Loop over roots
    do i=1,nroots

       ! Trim and sort the eigenvector for this state
       call trim_and_sort(ndet,vec(:,i),indx,ntrim,trim_thrsh,iwork)
       
       ! Fill in the surviving determinants for this state
       normsq=0.0d0
       do k=1,ntrim
          
          ! Update the squared norm
          normsq=normsq+vec(indx(k),i)**2

          ! Flag the determinant for survival
          idet(indx(k))=1
          
          ! Exit if:
          ! (1) we have hit the target squared norm, and;
          ! (2) this determinant is not degenerate with the next one
          if (k < ntrim) then
             diff=abs(vec(indx(k),i)-vec(indx(k+1),i))
          else
             diff=10*epsilon
          endif
          if (normsq >= targ .and. diff > epsilon) exit
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Set up the truncated determinant and eigenvector arrays
!----------------------------------------------------------------------
    !
    ! Number of surviving determinants
    !
    ndet_new=sum(idet)

    !
    ! Allocate arrays
    !
    allocate(det_new(n_int,2,ndet_new))
    det_new=0_ib
    allocate(vec_new(ndet_new,nroots))
    vec_new=0.0d0

    !
    ! Fill in the determinant bit strings
    !
    n=0
    ! Loop over the original set of determinants
    do k=1,ndet

       ! Are we at a surviving determinant?
       if (idet(k) == 1) then
          n=n+1
          det_new(:,:,n)=det(:,:,k)
       endif
       
    enddo

    !
    ! Fill in the eigenvectors
    !
    ! Loop over roots
    do i=1,nroots

       ! Surviving determinant counter
       n=0
       
       ! Loop over the original set of determinants
       do k=1,ndet

          ! Are we at a surviving determinant?
          if (idet(k) == 1) then
             n=n+1
             vec_new(n,i)=vec(k,i)
          endif
             
       enddo

    enddo

    return
    
  end subroutine truncate_wave_functions

!######################################################################
! trim_and_sort: for a given eigenvector, sorts a trimmed subset where
!                the trimming is based on an absolute value threshold,
!                thrsh
!######################################################################
  subroutine trim_and_sort(ndet,vec,indx,ntrim,thrsh,iwork)

    use constants
    use utils
    
    implicit none
    
    ! Dimensions
    integer(is), intent(in)    :: ndet

    ! Input eigenvector
    real(dp), intent(in)       :: vec(ndet)

    ! Indices of the determinants in order of increasing absolute
    ! coefficient value
    integer(is), intent(out)   :: indx(ndet)

    ! Size of the trimmed determinant basis
    integer(is), intent(out)   :: ntrim
    
    ! Trimming threshold
    real(dp), intent(in)       :: thrsh

    ! Work array
    integer(is), intent(inout) :: iwork(ndet)
    
    ! Everything else
    integer(is)                :: i,n
    integer(is), allocatable   :: imap(:)
    real(dp), allocatable      :: cabs(:)

!----------------------------------------------------------------------
! Determine the indices of the surviving determinants
!----------------------------------------------------------------------
    ! Initialisation
    iwork=0
        
    ! Loop over determinants
    do i=1,ndet
       
       if (abs(vec(i)) > thrsh) iwork(i)=1
       
    enddo

    ! Number of surviving determinants
    ntrim=sum(iwork)

!----------------------------------------------------------------------
! Fill in the truncated vector of absolute coefficient values
!----------------------------------------------------------------------
    allocate(cabs(ntrim))
    allocate(imap(ntrim))
    
    n=0
    do i=1,ndet
       if (iwork(i) == 1) then
          n=n+1
          cabs(n)=abs(vec(i))
          imap(n)=i
       endif
    enddo

!----------------------------------------------------------------------
! Sort the truncated vector of absolute coefficient values
!----------------------------------------------------------------------
    iwork=0
    call dsortindxa1('D',ntrim,cabs,iwork(1:ntrim))

!----------------------------------------------------------------------
! Fill in the output array of sorted determinant indices
!----------------------------------------------------------------------
    indx=0
    do i=1,ntrim
       indx(i)=imap(iwork(i))
    enddo
    
    return
    
  end subroutine trim_and_sort
  
!######################################################################
! symm_ortho: Orthonormalisation of a set of wave functions using
!             Lowdin's symmetric orthogonalisation
!######################################################################
  subroutine symm_ortho(n_int,ndet,nroots,vec)

    use constants
    use utils
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: n_int,ndet,nroots

    ! Eigenvectors
    real(dp), intent(inout) :: vec(ndet,nroots)

    ! Everything else
    integer(is)             :: i,j
    real(dp)                :: norm
    real(dp), allocatable   :: Smat(:,:),Sinvsq(:,:),vec_ortho(:,:)
        
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(Smat(nroots,nroots), Sinvsq(nroots,nroots), &
         vec_ortho(ndet,nroots))
    Smat=0.0d0; Sinvsq=0.0d0; vec_ortho=0.0d0
    
!----------------------------------------------------------------------
! Normalisation
!----------------------------------------------------------------------
    do i=1,nroots
       norm=dot_product(vec(:,i),vec(:,i))
       norm=sqrt(norm)
       vec(:,i)=vec(:,i)/norm
    enddo
    
!----------------------------------------------------------------------
! Overlap matrix
!----------------------------------------------------------------------
    call dgemm('T','N',nroots,nroots,ndet,1.0d0,vec,ndet,vec,ndet,&
         0.0d0,Smat,nroots)

!----------------------------------------------------------------------
! Inverse square root of the overlap matrix
!----------------------------------------------------------------------
    call invsqrt_matrix(Smat,Sinvsq,nroots)

!----------------------------------------------------------------------
! Orthogonalisation
!----------------------------------------------------------------------
    call dgemm('N','N',ndet,nroots,nroots,1.0d0,vec,ndet,Sinvsq,&
         nroots,0.0d0,vec_ortho,ndet)

    vec=vec_ortho
    
    return
    
  end subroutine symm_ortho
    
!######################################################################
! get_nel: Determines the number of electrons (total, alpha and beta)
!          in a given determinant d
!######################################################################
  subroutine get_nel(n_int,d,nel,nel_alpha,nel_beta)

    use constants
    
    implicit none

    integer(is), intent(in)  :: n_int
    integer(ib), intent(in)  :: d(n_int,2)
    integer(is), intent(out) :: nel,nel_alpha,nel_beta

    integer(is)              :: k

    !
    ! Number of alpha electrons
    !
    nel_alpha=0
    do k=1,n_int
       nel_alpha=nel_alpha+popcnt(d(k,1))
    enddo
    
    !
    ! Number of beta electrons
    !
    nel_beta=0
    do k=1,n_int
       nel_beta=nel_beta+popcnt(d(k,2))
    enddo

    !
    ! Total number of electrons
    !
    nel=nel_alpha+nel_beta
    
    return
    
  end subroutine get_nel

!######################################################################
! get_nel_string: Determines the number of electrons 
!                 in a given alpha or beta string 'string'
!######################################################################
  subroutine get_nel_string(n_int,string,nel)

    use constants
    
    implicit none

    integer(is), intent(in)  :: n_int
    integer(ib), intent(in)  :: string(n_int)
    integer(is), intent(out) :: nel

    integer(is)              :: k

    !
    ! Number of electrons in the string
    !
    nel=0
    do k=1,n_int
       nel=nel+popcnt(string(k))
    enddo
    
    return
    
  end subroutine get_nel_string
  
!######################################################################
! mo_occ_string: Given a single alpha/beta string, determines the
!                indices of the occupied MOs
!######################################################################
  subroutine mo_occ_string(n_int,string,dim,nocc,occ)

    use constants
    
    implicit none

    ! Alpha or beta string
    integer(is), intent(in)  :: n_int
    integer(ib), intent(in)  :: string(n_int)

    ! No. occupied orbitals
    integer(is), intent(out) :: nocc

    ! Indices of the occupied orbitals
    integer(is), intent(in)  :: dim
    integer(is), intent(out) :: occ(dim)

    ! Everything else
    integer(is)              :: ic,k,ipos
    integer(ib)              :: h

    !
    ! Initialisation
    !
    occ=0
    nocc=0

    !
    ! Get the indices of the occupied orbitals
    !
    ! Orbital counter
    ic=1

    ! Loop over bit string blocks
    do k=1,n_int

       ! Initialise the work array
       h=string(k)

       ! Get the occupied orbital indices for this block
       do while (h /= 0_ib)

          ! Number of trailing zeros left in h
          ipos=trailz(h)

          ! Index of the next occupied orbital
          occ(ic)=1+ipos+(k-1)*n_bits

          ! Clear the bits up to the occupied orbital in h
          h=ibclr(h,ipos)

          ! Increment the orbital counter
          ic=ic+1
          
       enddo
          
    enddo

    ! No. occupied orbitals
    nocc=ic-1
    
    return
    
  end subroutine mo_occ_string
  
!######################################################################
! anihilate_electron_string: annihilates the electron in an input
!                            alpha/beta string
!######################################################################
  function annihilate_electron_string(n_int,string,imo) &
       result(hole_string)

    use constants

    implicit none

    ! Function result
    integer(is), intent(in) :: n_int
    integer(ib)             :: hole_string(n_int)

    ! Input alpha/beta string
    integer(ib), intent(in) :: string(n_int)

    ! Index of the MO to be annihilated
    integer(is), intent(in) :: imo

    ! Everything else
    integer(is)             :: k,i

    ! Block index
    k=(imo-1)/n_bits+1

    ! Orbital position with the block
    i=imo-1-(k-1)*n_bits

    ! Annihilate the electron
    hole_string=string
    hole_string(k)=ibclr(string(k),i)
    
    return
    
  end function annihilate_electron_string

!######################################################################
! print_eigvec_det: Outputs a set of eigenvectors in the determinant
!                   representation. For the sake of brevity, the
!                   dominant determinants are output as their
!                   difference wrt the base determinant.
!######################################################################
  subroutine print_eigvec_det(ndet,nroots,n_int,vec,det)

    use constants
    use utils
    
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: ndet,nroots,n_int

    ! Eigenvectors
    real(dp), intent(in)     :: vec(ndet,nroots)

    ! Determinant bit strings
    integer(ib), intent(in)  :: det(n_int,2,ndet)

    ! Sorting
    integer(is), allocatable :: indx(:)
    real(dp), allocatable    :: absvec(:)

    ! Base determinant
    integer(ib)              :: det0(n_int,2)
    
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
! Set up the base determinant
! For now, we will take this as the dominant determinant in first state
!----------------------------------------------------------------------
    absvec=abs(vec(:,1))
    call dsortindxa1('D',ndet,absvec,indx)
    det0=det(:,:,indx(1))
        
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
          call get_diff_string(n_int,det0,det(:,:,indx(j)),string)          
          
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
  subroutine get_diff_string(n_int,det0,det,string)

    use constants

    implicit none

    ! Input determinant bit string
    integer(is), intent(in)        :: n_int
    integer(ib), intent(in)        :: det0(n_int,2)
    integer(ib), intent(in)        :: det(n_int,2)

    ! Output character string
    character(len=60), intent(out) :: string

    ! Particle and hole indices
    integer(ib)                    :: p(n_int,2),h(n_int,2)
    integer(is), parameter         :: nexmax=10
    integer(is)                    :: plist(nexmax,2),hlist(nexmax,2)
    integer(is)                    :: nexci(2)

    ! Everything else
    integer(is)                    :: i,i1,dim
        
!----------------------------------------------------------------------
! Get the bit string encodings of the particle and hole MOs relative
! base determinant
!----------------------------------------------------------------------
    call exc(n_int,det0,det,p,h)

!----------------------------------------------------------------------
! Get the alpha and beta excitation degrees
!----------------------------------------------------------------------
    nexci=0
    nexci(1)=exc_degree_string(n_int,det0(:,1),det(:,1))
    nexci(2)=exc_degree_string(n_int,det0(:,2),det(:,2))

!----------------------------------------------------------------------
! Get the lists of alpha and beta particle/hole MO indices
!----------------------------------------------------------------------
    ! alpha particles
    call list_from_bitstring(n_int,p(:,1),plist(:,1),nexmax)

    ! alpha holes
    call list_from_bitstring(n_int,h(:,1),hlist(:,1),nexmax)

    ! beta particles
    call list_from_bitstring(n_int,p(:,2),plist(:,2),nexmax)

    ! beta holes
    call list_from_bitstring(n_int,h(:,2),hlist(:,2),nexmax)

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
! exc: Given two determinants d1 and d2, returns bitstrings p and h
!      encoding, respectively, the indices of the particles and holes
!      created by the excitation linking d1 and d2
!######################################################################
  subroutine exc(n_int,d1,d2,p,h)

    use constants
    
    implicit none

    integer(is), intent(in)  :: n_int
    integer(ib), intent(in)  :: d1(n_int,2),d2(n_int,2)
    integer(ib), intent(out) :: p(n_int,2),h(n_int,2)
    integer(ib)              :: c
    integer(is)              :: ispin,k
    
    ! Loop over spins
    do ispin=1,2
       ! Loop over blocks
       do k=1,n_int
          c=ieor(d1(k,ispin),d2(k,ispin))
          p(k,ispin)=iand(c,d2(k,ispin))
          h(k,ispin)=iand(c,d1(k,ispin))
       enddo
    enddo
    
    return
    
  end subroutine exc

!######################################################################
! exc_degree_string: Computes the excitation degree between two alpha
!                    or beta strings s1 and s2
!######################################################################  
  function exc_degree_string(n_int,s1,s2) result(nexci)

    use constants
    
    implicit none

    integer(is)             :: nexci

    integer(is), intent(in) :: n_int
    integer(ib), intent(in) :: s1(n_int),s2(n_int)
    integer(is)             :: ispin,k

    ! Compute the number of creation and annihilation operators
    ! connecting s1 and s2
    nexci=popcnt(ieor(s1(1),s2(1)))
    do k=2,n_int
       nexci=nexci+popcnt(ieor(s1(k),s2(k)))
    enddo

    ! Divide by 2 to get the no. excitations connecting s1 and s2
    ! Note that this is equivalent to performing a right bit shift
    ! by 1 place
    nexci=shiftr(nexci,1)

    return
    
  end function exc_degree_string

!######################################################################
! list_from_bitstring: Given a bitstring I, fills in the array list
!                      with the indices of the set bits in I. Note
!                      that I is destroyed on output.
!######################################################################
  subroutine list_from_bitstring(n_int,I,list,listdim)

    use constants
    
    implicit none

    integer(is), intent(in)    :: n_int
    integer(is), intent(in)    :: listdim
    integer(ib), intent(inout) :: I(n_int)
    integer(is), intent(out)   :: list(listdim)
    integer(is)                :: ispin,k,n
    integer(is)                :: e
    
    ! Initialisation
    list=0
    n=1
    
    ! Loop over blocks
    do k=1,n_int
       
       ! Determine the indices of any set bits
       do while (I(k)/=0_ib)
          e=trailz(I(k))
          I(k)=ibclr(I(k),e)
          list(n)=e+(k-1)*n_bits+1
          n=n+1
       enddo
       
    enddo

    return
    
  end subroutine list_from_bitstring
  
!######################################################################
  
end module detfuncs
