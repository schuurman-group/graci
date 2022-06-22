!**********************************************************************
! Routines for the full diagonalisation of the small Hamiltonian
! matrices. Used for both reference space diagonalisation and the
! generation of spin-pure guess vectors for use in Davidson
! diagonalisation
!**********************************************************************

!######################################################################
! ref_diag: Diagonalisation of a reference space Hamiltonian
!######################################################################
#ifdef CBINDING
subroutine ref_diag(irrep,detscr,vecscr,nroots) &
     bind(c,name="ref_diag")
#else
subroutine ref_diag(irrep,detscr,vecscr,nroots)
#endif

  use constants
  use bitglobal
  use utils
  use iomod
  use detsort
  use slater_condon
  use timing

  implicit none

  ! Input irrep and requested number of roots
  integer(is), intent(in)    :: irrep
  integer(is), intent(inout) :: nroots
  
  ! Determinant scratch file number
  integer(is), intent(in)    :: detscr

  ! Determinant arrays
  integer(is)                :: ndet,offdima,offdimb
  integer(ib), allocatable   :: da(:,:,:),db(:,:,:)
  integer(is)                :: nsym(0:nirrep-1)
  integer(is)                :: nsym_sum(0:nirrep-1)
  integer(is)                :: nuniquea(0:nirrep-1),&
                                nuniqueb(0:nirrep-1)
  integer(is), allocatable   :: offseta(:,:),offsetb(:,:)
  integer(is), allocatable   :: mapab(:)

  ! Hamiltonian and S2 matrices
  real(dp), allocatable      :: hmat(:,:),s2mat(:,:)
  
  ! Spin purifcation variables
  integer(is)                :: info
  integer(is)                :: ncsf,counter
  integer(is), allocatable   :: csfindx(:)
  real(dp), allocatable      :: s2vec(:,:),s2val(:)
  real(dp)                   :: s2targ
  real(dp), parameter        :: s2thrsh=1e-10
  real(dp), allocatable      :: hmat_csf(:,:),det2csf(:,:),&
                                hval_csf(:),hvec_csf(:,:),&
                                hvec_det(:,:)
  
  ! I/O variables
  integer(is), intent(out)   :: vecscr
  integer(is)                :: iscratch
  character(len=60)          :: vecfile
  character(len=2)           :: amult,airrep

  ! Timing variables
  real(dp)                   :: tcpu_start,tcpu_end,twall_start,&
                                twall_end

  ! Everything else
  integer(is)                :: i

  integer(is)                :: j
  real(dp), allocatable      :: hvec(:,:),hval(:)
  real(dp)                   :: rand

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)
  
!----------------------------------------------------------------------
! Read in the determinants and associated information from disk
!----------------------------------------------------------------------
  call read_det_file_sorted(detscr,ndet,offdima,offdimb,nsym,da,db,&
       nuniquea,nuniqueb,offseta,offsetb,mapab)

!----------------------------------------------------------------------
! Return if there are no determinants for the current irrep
!----------------------------------------------------------------------
  if (nsym(irrep) == 0) then
     if (verbose) &
          write(6,'(/,x,a)') 'No reference space determinants of '&
          //trim(irreplbl(irrep,ipg))//' symmetry'
     nroots=0
     vecscr=-1
     deallocate(da)
     deallocate(db)
     deallocate(offseta)
     deallocate(offsetb)
     return
  endif
  
!----------------------------------------------------------------------
! Allocate and initialise the Hamiltonian and S2 matrix arrays
!----------------------------------------------------------------------
  allocate(hmat(nsym(irrep),nsym(irrep)))
  allocate(s2mat(nsym(irrep),nsym(irrep)))
  allocate(s2vec(nsym(irrep),nsym(irrep)))
  allocate(s2val(nsym(irrep)))
  hmat=0.0d0
  s2mat=0.0d0
  s2vec=0.0d0
  s2val=0.0d0

!----------------------------------------------------------------------
! Build the Hamiltonian and S2 matrices in the determinant basis
!----------------------------------------------------------------------
  call build_h_s2(irrep,hmat,s2mat,ndet,offdima,offdimb,nsym,da,db,&
       nuniquea,nuniqueb,offseta,offsetb,mapab)

!----------------------------------------------------------------------
! Diagonalise the S2 matrix
!----------------------------------------------------------------------
  call diag_matrix_real(s2mat,s2val,s2vec,nsym(irrep))

!----------------------------------------------------------------------
! Transform the Hamiltonian matrix to a CSF basis
!----------------------------------------------------------------------
  ! Determine how many eigenfunctions of the S2 matrix have the
  ! desired spin multiplicity
  ncsf=0
  s2targ=(dble(imult)-1.0d0)/2.0d0
  s2targ=s2targ*(s2targ+1.0d0)
  do i=1,nsym(irrep)
     if (abs(s2val(i)-s2targ) < s2thrsh) ncsf=ncsf+1
  enddo

  ! Allocate the CSF-related arrays
  allocate(csfindx(ncsf))
  allocate(hmat_csf(ncsf,ncsf))
  allocate(det2csf(nsym(irrep),ncsf))
  
  ! Determine the indices of the eigenfunctions of the S2 matrix with
  ! the desired spin multiplicity
  counter=0
  do i=1,nsym(irrep)
     if (abs(s2val(i)-s2targ) < s2thrsh) then
        counter=counter+1
        csfindx(counter)=i
     endif
  enddo

  ! Construct the determinant-to-CSF transformation matrix
  do i=1,ncsf
     det2csf(:,i)=s2vec(:,csfindx(i))
  enddo

  ! Transform the Hamiltonian matrix to the CSF basis
  hmat_csf=matmul(transpose(det2csf),matmul(hmat,det2csf))
  
!----------------------------------------------------------------------
! Diagonalise the Hamiltonian represented in the CSF basis
!----------------------------------------------------------------------
  ! Allocate the eigenpair arrays
  allocate(hval_csf(ncsf))
  allocate(hvec_csf(ncsf,ncsf))

  ! Diagonalise the CSF representation of the Hamiltonian
  call diag_matrix_real(hmat_csf,hval_csf,hvec_csf,ncsf)

!----------------------------------------------------------------------
! Transform the eigenvectors back to the determinant basis
!----------------------------------------------------------------------
  ! Allocate the determinant basis eigenvector array
  allocate(hvec_det(nsym(irrep),ncsf))

  ! Transform the eigenvectors
  hvec_det=matmul(det2csf,hvec_csf)

!!----------------------------------------------------------------------
!! S2 test
!!----------------------------------------------------------------------
!  allocate(hvec(nsym(irrep),nsym(irrep)),hval(nsym(irrep)))
!  
!  call diag_matrix_real(hmat,hval,hvec,nsym(irrep))
!
!  do i=1,nsym(irrep)
!     write(6,'(F10.7)') dot_product(hvec(:,i),matmul(s2mat,hvec(:,i)))
!  enddo
!
!  stop
  
!----------------------------------------------------------------------
! Write the eigenpairs to disk. Note that we will only save the
! requested number, nroots, of eigenpairs.
!----------------------------------------------------------------------
  ! Decrease nroots if it is greater than the number of CSFs
  if (nroots > ncsf) then
     if (verbose) &
          write(6,'(/,x,a,x,i0)') 'Adjusting the number of reference '&
          //'space roots of '//trim(irreplbl(irrep,ipg))&
          //' symmetry to',ncsf
     nroots=ncsf
  endif
     
  ! Register a new scratch file
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name('refvec.mult'//trim(amult)//&
       '.irrep'//trim(airrep),vecfile)
  call register_scratch_file(vecscr,vecfile)
  
  ! Open the scratch file
  iscratch=scrunit(vecscr)
  open(iscratch,file=scrname(vecscr),form='unformatted',&
       status='unknown')

  ! Dimensions
  write(iscratch) nsym(irrep)
  write(iscratch) nroots

  ! Eigenvalues
  write(iscratch) hval_csf(1:nroots)

  ! Eigenvectors
  write(iscratch) hvec_det(:,1:nroots)
  
  ! Close the scratch file
  close(iscratch)
  
!----------------------------------------------------------------------
! Deallocate arrays  
!----------------------------------------------------------------------
  deallocate(da)
  deallocate(db)
  deallocate(offseta)
  deallocate(offsetb)
  deallocate(mapab)
  deallocate(hmat)
  deallocate(s2mat)
  deallocate(s2vec)
  deallocate(s2val)
  deallocate(csfindx)
  deallocate(hmat_csf)
  deallocate(det2csf)
  deallocate(hval_csf)
  deallocate(hvec_csf)
  deallocate(hvec_det)
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
!  call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
!       'full_diag')

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine ref_diag

!######################################################################
! build_hmat: Builds the full Hamiltonian and S2 matrices
!######################################################################
subroutine build_h_s2(irrep,hmat,s2mat,ndet,offdima,offdimb,nsym,&
     da,db,nuniquea,nuniqueb,offseta,offsetb,mapab)

  use constants
  use bitglobal
  use bitutils
  use slater_condon
  use utils
  
  implicit none

  ! Determinant arrays
  integer(is), intent(in) :: irrep
  integer(is), intent(in) :: ndet,offdima,offdimb
  integer(is), intent(in) :: nsym(0:nirrep-1)
  integer(ib), intent(in) :: da(n_int,2,ndet),db(n_int,2,ndet)
  integer(is), intent(in) :: nuniquea(0:nirrep-1),&
                             nuniqueb(0:nirrep-1)
  integer(is), intent(in) :: offseta(offdima,0:nirrep-1),&
                             offsetb(offdimb,0:nirrep-1)
  integer(is), intent(in) :: mapab(ndet)

  ! Hamiltonian and S2 matrices
  real(dp), intent(out)   :: hmat(nsym(irrep),nsym(irrep))
  real(dp), intent(out)   :: s2mat(nsym(irrep),nsym(irrep))
  
  ! Hamiltonian and S2 matrix build variables
  integer(is)              :: i,i1,i2,j,n
  integer(is)              :: a,b1,b2,b,a1,a2
  integer(is)              :: ibra,iket,idet
  integer(is)              :: nexci
  integer(ib)              :: p(n_int,2),h(n_int,2)
  integer(is), parameter   :: maxex=2
  integer(is)              :: plist(maxex,2),hlist(maxex,2)
  integer(ib)              :: phase_mask(n_int,2)
  integer(is)              :: phase
  integer(is)              :: occlist(nmo,2)
  integer(is)              :: nsym_sum(0:nirrep-1)

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
  hmat=0.0d0
  s2mat=0.0d0

!----------------------------------------------------------------------
! Cumulative numbers of determinants in the symmetry blocks
!----------------------------------------------------------------------
  nsym_sum=0
  nsym_sum(1:nirrep-1)=nsym(0:nirrep-2)
  do i=1,nirrep-1
     nsym_sum(i)=nsym_sum(i)+nsym_sum(i-1)
  enddo
  
!----------------------------------------------------------------------
! (1) On-diagonal elements
!----------------------------------------------------------------------
  ! Loop over unique alpha strings for the current irrep
  do a=1,nuniquea(irrep)

     ! Loop over determinants with the current alpha
     ! string
     do b=offseta(a,irrep),offseta(a+1,irrep)-1

        ! Determinant index
        idet=b-nsym_sum(irrep)

        ! Get the indices occupied alpha and beta spin orbitals in
        ! the ket determinant
        call occ_orbs(da(:,:,b),occlist)
        
        ! On-diagonal element of H
        hmat(idet,idet)=hii(occlist)

        ! On-diagonal element of S2
        s2mat(idet,idet)=s2ii(da(:,:,b))

     enddo

  enddo

!----------------------------------------------------------------------
! (2) Pairs of determinants linked by pure beta-spin single and
!     double excitations
!----------------------------------------------------------------------
  ! Initialisation (important for phase factor testing)
  hlist=0
  plist=0

  ! Loop over unique alpha strings for the current irrep
  do a=1,nuniquea(irrep)

     ! Loop over pairs of determinants with the current alpha
     ! string
     !
     ! ket determinants
     do b1=offseta(a,irrep),offseta(a+1,irrep)-1
        
        ! Compute the phase mask for the ket determinant
        call phasemask(da(:,:,b1),phase_mask)

        ! Get the indices occupied alpha and beta spin orbitals in
        ! the ket determinant
        call occ_orbs(da(:,:,b1),occlist)

        ! ket index
        iket=b1-nsym_sum(irrep)
        
        ! bra determinants
        do b2=b1+1,offseta(a+1,irrep)-1
           
           ! Compute the excitation degree between the beta parts
           ! of two determinants
           nexci=exc_degree_string(da(:,ibeta,b2),da(:,ibeta,b1))
           
           ! Cycle if the excitation degree is greater than 2
           if (nexci>2) cycle

           ! bra index
           ibra=b2-nsym_sum(irrep)
           
           !********************************************************
           ! Note that if we are here, then we are considering
           ! an excitation between determinants
           ! |b1> -> |b2> = O T |b1> (not the otherway around!)
           !********************************************************
           
           ! Determine the indices of the beta orbitals involved
           ! in the excitation           
           call exc(da(:,:,b1),da(:,:,b2),p,h)
           call list_from_bitstring(p(:,ibeta),plist(:,ibeta),maxex)
           call list_from_bitstring(h(:,ibeta),hlist(:,ibeta),maxex)
           
           ! Compute the phase factor
           phase=phase_pure_exc(phase_mask(:,ibeta),hlist(:,ibeta),&
                plist(:,ibeta),maxex,nexci)
           
           ! Off-diagonal Hamiltonian matrix element
           select case(nexci)
           case(1)
              hmat(ibra,iket)=hij_single(ibeta,phase,occlist,&
                   hlist,plist,maxex)
           case(2)
              hmat(ibra,iket)=hij_double_same_spin(ibeta,phase,&
                   occlist,hlist,plist,maxex)
           end select
           hmat(iket,ibra)=hmat(ibra,iket)

        enddo
     enddo

  enddo

!----------------------------------------------------------------------
! (3) Pairs of determinants linked by pure alpha-spin single and
!     double excitations
!----------------------------------------------------------------------
  ! Initialisation (important for phase factor testing)
  hlist=0
  plist=0
  
  ! Loop over unique beta strings for the current irrep
  do b=1,nuniqueb(irrep)

     ! Loop over pairs of determinants with the current beta
     ! string
     !
     ! ket determinants
     do a1=offsetb(b,irrep),offsetb(b+1,irrep)-1
        
        ! Compute the phase mask for the ket determinant
        call phasemask(db(:,:,a1),phase_mask)

        ! Get the indices occupied alpha and beta spin orbitals in
        ! the ket determinant
        call occ_orbs(db(:,:,a1),occlist)

        ! ket index
        iket=mapab(a1)-nsym_sum(irrep)
                
        ! bra determinants
        do a2=a1+1,offsetb(b+1,irrep)-1

           ! Compute the excitation degree between the alpha parts
           ! of two determinants
           nexci=exc_degree_string(db(:,ialpha,a2),db(:,ialpha,a1))
           
           ! Cycle if the excitation degree is greater than 2
           if (nexci>2) cycle

           ! bra index
           ibra=mapab(a2)-nsym_sum(irrep)

           !********************************************************
           ! Note that if we are here, then we are considering
           ! an excitation between determinants
           ! |a1> -> |a2> = O T |a1> (not the otherway around!)
           !********************************************************

           ! Determine the indices of the alpha orbitals involved
           ! in the excitation
           call exc(db(:,:,a1),db(:,:,a2),p,h)
           call list_from_bitstring(p(:,ialpha),plist(:,ialpha),maxex)
           call list_from_bitstring(h(:,ialpha),hlist(:,ialpha),maxex)

           ! Compute the phase factor
           phase=phase_pure_exc(phase_mask(:,ialpha),hlist(:,ialpha),&
                plist(:,ialpha),maxex,nexci)
                      
           ! Off-diagonal Hamiltonian matrix element
           select case(nexci)
           case(1)
              hmat(ibra,iket)=hij_single(ialpha,phase,occlist,&
                   hlist,plist,maxex)
           case(2)
              hmat(ibra,iket)=hij_double_same_spin(ialpha,phase,&
                   occlist,hlist,plist,maxex)
           end select
           hmat(iket,ibra)=hmat(ibra,iket)
           
        enddo
        
     enddo
     
  enddo

!----------------------------------------------------------------------
! (4) Pairs of determinants linked by alpha,beta double excitations
!----------------------------------------------------------------------
  ! Initialisation (important for phase factor testing)
  hlist=0
  plist=0

  ! Loop over pairs of unique alpha strings for the current irrep
  do a1=1,nuniquea(irrep)
     do a2=a1+1,nuniquea(irrep)

        ! Compute the excitation degree between the alpha strings
        ! and cycle if not equal to 1
        nexci=exc_degree_string(da(:,ialpha,offseta(a2,irrep)),&
             da(:,ialpha,offseta(a1,irrep)))
        if (nexci /= 1) cycle

        ! Loop over pairs of determinants with the current alpha
        ! strings
        !
        ! ket determinants
        do b1=offseta(a1,irrep),offseta(a1+1,irrep)-1

           ! Compute the phase mask for the ket determinant
           call phasemask(da(:,:,b1),phase_mask)

           ! Get the indices occupied alpha and beta spin orbitals in
           ! the ket determinant
           call occ_orbs(da(:,:,b1),occlist)

           ! ket index
           iket=b1-nsym_sum(irrep)
                      
           ! bra determinants
           do b2=offseta(a2,irrep),offseta(a2+1,irrep)-1
              
              ! Compute the excitation degree between the beta
              ! parts of the determinants
              nexci=exc_degree_string(da(:,ibeta,b2),da(:,ibeta,b1))

              ! Cycle if the excitation degree is not equal to 1
              if (nexci /= 1) cycle

              ! bra index
              ibra=b2-nsym_sum(irrep)
              
              !********************************************************
              ! Note that if we are here, then we are considering
              ! an excitation between determinants
              ! |b1> -> |b2> = O T |b1> (not the otherway around!)
              !********************************************************

              ! Determine the indices of the alpha beta orbitals involved
              ! in the excitation
              call exc(da(:,:,b1),da(:,:,b2),p,h)
              call list_from_bitstring(p(:,ialpha),plist(:,ialpha),maxex)
              call list_from_bitstring(h(:,ialpha),hlist(:,ialpha),maxex)
              call list_from_bitstring(p(:,ibeta),plist(:,ibeta),maxex)
              call list_from_bitstring(h(:,ibeta),hlist(:,ibeta),maxex)

              ! Compute the phase factor
              phase=phase_pure_exc(phase_mask(:,ialpha),hlist(:,ialpha),&
                   plist(:,ialpha),maxex,nexci)&
                   *phase_pure_exc(phase_mask(:,ibeta),hlist(:,ibeta),&
                   plist(:,ibeta),maxex,nexci)

              ! Off-diagonal Hamiltonian matrix element
              hmat(ibra,iket)=hij_double_diff_spin(phase,occlist,&
                   hlist,plist,maxex)
              hmat(iket,ibra)=hmat(ibra,iket)
              
              ! Compute the S2 matrix element
              s2mat(ibra,iket)=s2ij(phase,hlist,plist)
              s2mat(iket,ibra)=s2mat(ibra,iket)
              
           enddo

        enddo
           
     enddo
  enddo
  
  return
  
end subroutine build_h_s2
