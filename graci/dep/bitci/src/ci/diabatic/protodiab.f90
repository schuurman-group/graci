!**********************************************************************
! Routines for the calculation of prototype diabatic and related states
!**********************************************************************
module protodiab

  implicit none

contains

!######################################################################
! get_pds_basis: Given a set of diabatic states from a previous
!                geometry, computes a set of 'prototype diabatic
!                states' via:
!
!                (1) the projection of the previous geometry diabatic
!                    states onto the current geometry reference space
!                    states to form 'precursor states', followed by;
!
!                (2) Lowdin's symmetric orthonormalisation of the
!                    precursor states to yield the prototype diabatic
!                    states
!######################################################################
  subroutine get_pds_basis(cfg,refdim,nvec,vec0,nmoR0,n_intR0,ndetR0,&
       nrootsR0,detR0,vecR0,smoR0,ncore,icore,lfrzcore,adtR0,vec_pds)
    
    use constants
    use bitglobal
    use conftype
    use csf2det
    use utils
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg
    
    ! Dimension of the ref CSF space
    integer(is), intent(in)  :: refdim
    
    ! No. ref space states
    integer(is), intent(in)  :: nvec
    
    ! Ref space wave functions
    real(dp), intent(in)     :: vec0(refdim,nvec)

    ! Previous geometry eigensates expressed in a Slater determinant
    ! basis
    integer(is), intent(in)  :: n_intR0,ndetR0,nrootsR0
    integer(ib), intent(in)  :: detR0(n_intR0,2,ndetR0)
    real(dp), intent(in)     :: vecR0(ndetR0,nrootsR0)
    
    ! MO overlaps
    integer(is), intent(in)  :: nmoR0
    real(dp), intent(in)     :: smoR0(nmoR0,nmo)

    ! Frozen core orbital info
    integer(is), intent(in)  :: ncore
    integer(is), intent(in)  :: icore(ncore)
    logical(is), intent(in)  :: lfrzcore

    ! ADT matrix of the previous geometry
    real(dp), intent(in)     :: adtR0(nrootsR0,nrootsR0)
    
    ! Protoype diabatic states in the ref CSF basis
    real(dp), intent(out)    :: vec_pds(refdim,nrootsR0)

    ! Dummy ref space configuration derived type
    type(mrcfg)              :: cfg_ref

    ! Precursor states in the ref state basis
    real(dp), allocatable    :: precoe(:,:)
    
    ! Prototype diabatic states in the ref state basis
    real(dp), allocatable    :: pdscoe(:,:)

    ! Ref space wave functions in the determinant basis
    integer(is)              :: ndet_ref
    integer(ib), allocatable :: det_ref(:,:,:)
    real(dp), allocatable    :: vec0_det(:,:)

    ! Ref-state - previous-geometry-diabatic-state overlaps
    integer(is)              :: npairs
    integer(is), allocatable :: ipairs(:,:)
    real(dp), allocatable    :: Sij(:),Smat(:,:)
    real(dp), allocatable    :: smoT(:,:)
    real(dp)                 :: normthrsh
    logical                  :: lprint
    
    ! Everything else
    integer(is)              :: i,j,n
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(precoe(nvec,nrootsR0))
    precoe=0.0d0
    
    allocate(pdscoe(nvec,nrootsR0))
    pdscoe=0.0d0

    npairs=nvec*nrootsR0
    allocate(Sij(npairs))
    allocate(Smat(nvec,nrootsR0))
    allocate(smoT(nmo,nmoR0))
    allocate(ipairs(npairs,2))
    Sij=0.0d0
    Smat=0.0d0
    smoT=0.0d0
    ipairs=0
    
!----------------------------------------------------------------------
! Get the determinant representation of the zeroth-order eigenfunctions
!----------------------------------------------------------------------
    ! Set up a dummy MRCI configuration derived data type to hold the
    ! ref conf info
    allocate(cfg_ref%conf0h(cfg%n_int_I,2,cfg%n0h))
    allocate(cfg_ref%sop0h(cfg%n_int_I,2,cfg%n0h))
    cfg_ref%n0h=cfg%n0h
    cfg_ref%n_int_I=cfg%n_int_I
    cfg_ref%n1I=0
    cfg_ref%n2I=0
    cfg_ref%n1E=0
    cfg_ref%n2E=0
    cfg_ref%n1I1E=0
    cfg_ref%conf0h=cfg%conf0h
    cfg_ref%sop0h=cfg%sop0h
    
    ! Determine the dimension of the reference space determinant basis
    call get_detdim(cfg_ref,ndet_ref)
    
    ! Allocate arrays
    allocate(det_ref(n_int,2,ndet_ref))
    allocate(vec0_det(ndet_ref,nvec))
    det_ref=0_ib
    vec0_det=0.0d0
    
    ! Compute the determinant representation of the reference space
    ! wave functions
    call det_trans(cfg_ref,cfg%m2c,nvec,refdim,ndet_ref,vec0,vec0_det,&
         det_ref)

!----------------------------------------------------------------------
! Compute the overlaps of the reference space wave functions with
! the wave functions of the previous geometry
!----------------------------------------------------------------------
! In the following,
!
! Smat(i,j) = <psi_i^(0)(R_n)|psi_j(R_n-1)>
!
! i.e., Smat(:,j) is the projection of the previous geometry
!       wave functions in terms onto the space spanned by the current
!       geometry reference space wave functions
!----------------------------------------------------------------------
    ! Truncation threshold
    normthrsh=0.999d0

    ! Fill in the array of bra-ket overlaps required
    n=0
    do i=1,nvec
       do j=1,nrootsR0
          n=n+1
          ipairs(n,1)=i
          ipairs(n,2)=j
       enddo
    enddo

    ! MO overlaps
    smoT=transpose(smoR0)
    
    ! Compute the overlaps
    lprint=.true.
    call overlap(nmo,nmoR0,n_int,n_intR0,ndet_ref,ndetR0,nvec,&
         nrootsR0,det_ref,detR0,vec0_det,vecR0,smoT,normthrsh,&
         ncore,icore,lfrzcore,npairs,Sij,ipairs,lprint)
  
    ! Put the overlaps into a more useful form
    do n=1,npairs
       i=ipairs(n,1)
       j=ipairs(n,2)
       Smat(i,j)=Sij(n)
    enddo

!----------------------------------------------------------------------
! Transform the overlaps using the previous geometry ADT matrix
! to yield the precursor states
!----------------------------------------------------------------------
    precoe=matmul(Smat,adtR0)
  
!----------------------------------------------------------------------
! Output the squared norms of projections of the previous geometry
! diabatic wave functions onto the space spanned by the current
! geometry reference space wave functions
!----------------------------------------------------------------------
    ! Table header
    write(6,'(/,x,21a)') ('-', i=1,21)
    write(6,'(4x,a)') 'Precursor state'
    write(6,'(5x,a)') 'squared norms'
    write(6,'(x,21a)') ('-', i=1,21)
    write(6,'(2x,a)') 'I  | <phi#_I|phi#_I>'
    write(6,'(x,21a)') ('-', i=1,21)

    ! Squared norms of the precursor states
    do i=1,nrootsR0
       write(6,'(x,i3,x,a,x,F10.7)') &
            i,'|',dot_product(precoe(:,i),precoe(:,i))
    enddo

    ! Table footer
    write(6,'(x,21a)') ('-', i=1,21)

!----------------------------------------------------------------------
! Lowdin's symmetric orthonormalisation of the precursor states to
! yield the prototype diabatic states in the ref space state basis
!----------------------------------------------------------------------
    pdscoe=precoe
    call symm_ortho(nvec,nrootsR0,pdscoe)

!----------------------------------------------------------------------
! Prototype diabatic states in the ref space CSF basis
!----------------------------------------------------------------------
    vec_pds(1:refdim,1:nrootsR0)=matmul(&
         vec0(1:refdim,1:nvec),&
         pdscoe(1:nvec,1:nrootsR0))

    return
    
  end subroutine get_pds_basis
  
!######################################################################
! get_complememt_basis: Given the sets of prototype diabatic states
!                       and ref space states, computes the 'complement'
!                       states via:
!
!                       (1) projection of the reference space states
!                           onto the space orthogonal to the prototype
!                           diabatic states, followed by;
!
!                       (2) the dropping of the null space
!######################################################################
  subroutine get_complement_basis(refdim,nrootsR0,nvec,ncomp,vec0,&
       vec_pds,vec_comp)

    use constants
    use bitglobal
    use utils
    use iomod

    implicit none

    ! Dimensions
    integer(is), intent(in) :: refdim,nrootsR0,nvec,ncomp

    ! Ref space eigenstates
    real(dp), intent(in)    :: vec0(refdim,nvec)

    ! Prototype diabatic states in the ref space CSF basis
    real(dp), intent(in)    :: vec_pds(refdim,nrootsR0)

    ! Complement states in the ref space CSF basis
    real(dp), intent(out)   :: vec_comp(refdim,ncomp)
    
    ! Projectors
    real(dp), allocatable   :: Pmat(:,:),Qmat(:,:)

    ! Q-projected ref space states
    real(dp), allocatable   :: Qvec0(:,:)

    ! Gram matrix and its eigenpairs
    real(dp), allocatable   :: gram(:,:)
    real(dp), allocatable   :: eigvec(:,:),eigval(:)
    
    ! Everything else
    integer(is)             :: i,j
    integer(is)             :: nulldim
    real(dp), parameter     :: nullthrsh=1e-10_dp
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(Pmat(refdim,refdim))
    allocate(Qmat(refdim,refdim))
    Pmat=0.0d0
    Qmat=0.0d0

    allocate(Qvec0(refdim,nvec))
    Qvec0=0.0d0

    allocate(gram(nvec,nvec))
    allocate(eigvec(nvec,nvec))
    allocate(eigval(nvec))
    gram=0.0d0
    eigvec=0.0d0
    eigval=0.0d0
        
!----------------------------------------------------------------------
! Set up the projector onto the orthogonal complement of the prototype
! diabatic states
!----------------------------------------------------------------------
    ! Projector, P, onto the space spanned by the prototype diabatic
    ! states
    Pmat=matmul(vec_pds,transpose(vec_pds))

    ! Q = 1 - P
    Qmat=-Pmat
    do i=1,refdim
       Qmat(i,i)=Qmat(i,i)+1.0d0
    enddo

!----------------------------------------------------------------------
! Compute the projection of the ref space states onto the orthogonal
! complement of the prototype diabatic states
!----------------------------------------------------------------------
    Qvec0=matmul(Qmat,vec0)

!----------------------------------------------------------------------
! Compute the overlaps between the projected ref space states
! (the Gram matrix)
!----------------------------------------------------------------------
    gram=matmul(transpose(Qvec0),Qvec0)

!----------------------------------------------------------------------
! Diagonalise the Gram matrix
!----------------------------------------------------------------------
    call diag_matrix_real(gram,eigval,eigvec,nvec)

!----------------------------------------------------------------------
! Check on the dimension of the null space
!----------------------------------------------------------------------
    nulldim=0
    do i=1,nvec
       if (abs(eigval(i)) < nullthrsh) nulldim=nulldim+1
    enddo

    if (nulldim /= nrootsR0) then
       errmsg='Error in get_complement_basis: '&
            //'incorrect null space dimension'
       call error_control
    endif

!----------------------------------------------------------------------
! Form the complement states by dropping the null space vectors
!----------------------------------------------------------------------
    vec_comp=matmul(Qvec0,eigvec(:,nulldim+1:nvec))

    return
    
  end subroutine get_complement_basis

!######################################################################
! get_model_basis: Block diagonalisation of the Hamiltonian matrix
!                  between the prototype diabatic state and complement
!                  state blocks
!
!                  The transformed states in the prototype diabatic
!                  state block will be taken as the model space states
!######################################################################
  subroutine get_model_basis(nvec,nrootsR0,ncomp,refdim,E0,vec_pds,&
       vec_comp,vec_mod)

    use constants
    use bitglobal
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: nvec,nrootsR0,ncomp
    integer(is), intent(in) :: refdim

    ! Reference space eigenvalues
    real(dp), intent(in)    :: E0(nvec)

    ! Prototype diabatic states
    real(dp), intent(in)    :: vec_pds(refdim,nrootsR0)

    ! Complement states
    real(dp), intent(in)    :: vec_comp(refdim,ncomp)

    ! Model space states
    real(dp), intent(out)   :: vec_mod(refdim,nrootsR0)
    
    print*,''
    print*,'We need the (rectangular) ref-to-PDS and ref-to-comp'&
         //' transformations...'
    print*,''
    stop
    
    return
    
  end subroutine get_model_basis
    
!######################################################################
  
end module protodiab
