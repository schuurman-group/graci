!**********************************************************************
! Routines for the calculation of quasi-diabatic model states
!**********************************************************************
module model_states

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
       nrootsR0,detR0,vecR0,smoR0,ncore,icore,lfrzcore,vec_pds,&
       normthrsh)
    
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

    ! Previous geometry diabatic states expressed in a Slater
    ! determinant basis
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

    ! Protoype diabatic states in the ref CSF basis
    real(dp), intent(out)    :: vec_pds(refdim,nrootsR0)

    ! Norm-based wave function truncation threshold
    real(dp), intent(in)     :: normthrsh

    ! Determinant screening threshold
    real(dp)                 :: hthrsh
    
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
! the diabatic wave functions of the previous geometry
!----------------------------------------------------------------------
! In the following,
!
! Smat(i,j) = <psi_i^(0)(R_n)|psi_j(R_n-1)>
!
! i.e., Smat(:,j) is the projection of the previous geometry
!       diabatic wave functions in terms onto the space spanned by the
!       current geometry reference space wave functions
!----------------------------------------------------------------------
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
    hthrsh=1e-6_dp
    call overlap(nmo,nmoR0,n_int,n_intR0,ndet_ref,ndetR0,nvec,&
         nrootsR0,det_ref,detR0,vec0_det,vecR0,smoT,normthrsh,&
         hthrsh,ncore,icore,lfrzcore,npairs,Sij,ipairs,lprint)
    
!----------------------------------------------------------------------
! Fill in the precursor state coefficients
!----------------------------------------------------------------------
    do n=1,npairs
       i=ipairs(n,1)
       j=ipairs(n,2)
       precoe(i,j)=Sij(n)
    enddo
    
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
    call dgemm('N','N',refdim,nrootsR0,nvec,1.0d0,vec0,refdim,&
         pdscoe,nvec,0.0d0,vec_pds,refdim)
    
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
  subroutine get_model_basis(nvec,nrootsR0,ncomp,refdim,E0,vec0,&
       vec_pds,vec_comp,vec_mod,ref2mod)

    use constants
    use bitglobal
    use utils
    
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: nvec,nrootsR0,ncomp
    integer(is), intent(in)  :: refdim

    ! Reference space eigenpairs
    real(dp), intent(in)     :: E0(nvec)
    real(dp), intent(in)     :: vec0(refdim,nvec)
    
    ! Prototype diabatic states
    real(dp), intent(in)     :: vec_pds(refdim,nrootsR0)

    ! Complement states
    real(dp), intent(in)     :: vec_comp(refdim,ncomp)

    ! Model space states
    real(dp), intent(out)    :: vec_mod(refdim,nrootsR0)

    ! Ref-state-to-model-state transformation
    real(dp), intent(out)    :: ref2mod(nvec,nrootsR0)
    
    ! Transformation matrices
    real(dp), allocatable    :: ref2pds(:,:)
    real(dp), allocatable    :: ref2comp(:,:)
    real(dp), allocatable    :: ref2tot(:,:)

    ! Hamiltonian matrix
    real(dp), allocatable    :: H0(:,:),Htot(:,:)

    ! Block diagonalisation transformation
    integer(is), allocatable :: indx(:)
    real(dp), allocatable    :: norm_pds(:),Vtmp(:,:)
    real(dp), allocatable    :: eigval(:),V(:,:),VBDD(:,:)
    real(dp), allocatable    :: X(:,:),Y(:,:),invsqrtY(:,:),T(:,:)
    real(dp), allocatable    :: vec_tot(:,:)    
        
    ! Everything else
    integer(is)              :: i,j
    
!----------------------------------------------------------------------
! Set up the transformations from the ref state to prototype diabatic
! and complement bases
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(ref2pds(nvec,nrootsR0))
    allocate(ref2comp(nvec,ncomp))
    allocate(ref2tot(nvec,nvec))
    
    ! Ref-to-prototype-diabatic transformation
    ref2pds=matmul(transpose(vec0),vec_pds)
    
    ! Ref-to-complement transformation
    ref2comp=matmul(transpose(vec0),vec_comp)

    ! Transformation for the union of the two bases
    ref2tot(:,1:nrootsR0)=ref2pds
    ref2tot(:,nrootsR0+1:nvec)=ref2comp
    
!----------------------------------------------------------------------
! Hamiltonian matrix in the PDS U COMP basis
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(H0(nvec,nvec))
    allocate(Htot(nvec,nvec))
    
    ! Hamiltonian matrix in the ref state basis
    ! (subtraction of E_SCF gives the true eigenvalues)
    H0=0.0d0
    do i=1,nvec
       H0(i,i)=E0(i)-escf
    enddo

    ! Hamiltonian matrix in the PDS U COMP basis
    Htot=matmul(H0,ref2tot)
    Htot=matmul(transpose(ref2tot),Htot)

!----------------------------------------------------------------------
! Eigenvectors of the Hamiltonian matrix in the PDS U COMP basis
!----------------------------------------------------------------------
    allocate(eigval(nvec))
    allocate(V(nvec,nvec))

    call diag_matrix_real(Htot,eigval,V,nvec)

!----------------------------------------------------------------------
! Re-order the eigenvectors s.t. the ones with maximum projection
! onto the PDS space come first
!----------------------------------------------------------------------
    allocate(indx(nvec))
    allocate(norm_pds(nvec))
    allocate(Vtmp(nvec,nvec))
    
    ! Norms of the eigenvectors projected onto the PDS space
    do i=1,nvec
       norm_pds(i)=sqrt(dot_product(V(1:nrootsR0,i),V(1:nrootsR0,i)))
    enddo

    ! Sort the eigenvectors by the norm of their projection onto
    ! the PDS space
    call dsortindxa1('D',nvec,norm_pds,indx)
    
    ! Rearrange the eigenvectors
    do i=1,nvec
       Vtmp(:,i)=V(:,indx(i))
    enddo
    
    V=Vtmp
    
!----------------------------------------------------------------------
! Block diagonalisation transformation T satisfying ||T-1|| = min,
!
! T = V V_BDD^T (V_BDD V_BDD^T)^-1/2,
!
! where V is the matrix of eigenvectors of the Hamiltonian matrix
! and V_BDD its block diagonal part
!----------------------------------------------------------------------
    allocate(VBDD(nvec,nvec))
    allocate(X(nvec,nvec))
    allocate(Y(nvec,nvec))
    allocate(invsqrtY(nvec,nvec))
    allocate(T(nvec,nvec))

    ! V_BDD
    VBDD=0.0d0
    do i=1,nrootsR0
       do j=1,nrootsR0
          VBDD(i,j)=V(i,j)
       enddo
    enddo
    do i=nrootsR0+1,nvec
       do j=nrootsR0+1,nvec
          VBDD(i,j)=V(i,j)
       enddo
    enddo
    
    ! X = V V_BDD^T
    X=matmul(V,transpose(VBDD))
    
    ! Y = V_BDD V_BDD^T
    Y=matmul(VBDD,transpose(VBDD))

    ! Y^-1/2
    call invsqrt_matrix(Y,invsqrtY,nvec)
    
    ! T = X Y^-1/2
    T=matmul(X,invsqrtY)

!----------------------------------------------------------------------
! Compute the model space states
!----------------------------------------------------------------------
    ! PDS + COMP states
    allocate(vec_tot(refdim,nvec))
    vec_tot(:,1:nrootsR0)=vec_pds
    vec_tot(:,nrootsR0+1:nvec)=vec_comp

    ! Apply the block diagonalisation transformation
    vec_tot=matmul(vec_tot,T)

    ! Extract the model space states
    vec_mod=vec_tot(:,1:nrootsR0)

!----------------------------------------------------------------------
! Set up the ref-state-to-model-state transformation
!----------------------------------------------------------------------
    ref2mod(1:nvec,1:nrootsR0)=matmul(transpose(vec0),vec_mod)

    return
    
  end subroutine get_model_basis
    
!######################################################################
  
end module model_states
