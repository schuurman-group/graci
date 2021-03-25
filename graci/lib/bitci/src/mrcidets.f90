!**********************************************************************
! Routines for the generation of the Slater determinants of an MRCI
! wavefunction
!**********************************************************************

!######################################################################
! generate_mrci_dets: generates single or single+double excitations
!                     out of some set of a given set of reference
!                     space determinants
!######################################################################
#ifdef CBINDING
subroutine generate_mrci_dets(dref,nref,ndet,class_type,selection_type) &
     bind(c,name="generate_mrci_dets")
#else
subroutine generate_mrci_dets(dref,nref,ndet,class_type,selection_type)
#endif

  use constants
  use bitglobal
  use bitutils
  use utils
  use dethash
  use iomod
  use timing
  
  implicit none

  integer(is), intent(in)      :: nref
  integer(ib), intent(in)      :: dref(n_int,2,nref)
  integer(is), intent(out)     :: ndet
  character(len=*), intent(in) :: class_type
  character(len=*), intent(in) :: selection_type

  integer(is), allocatable     :: occ(:,:,:),unocc(:,:,:)
  integer(is)                  :: nea,neb,nha,nhb
  integer(is)                  :: i,j,k,iref
  integer(is)                  :: ia1,iap1,ia,iap,ika,ikap
  integer(is)                  :: ic1,icp1,ic,icp,ikc,ikcp
  integer(is)                  :: iexc1,iexc2
  real(dp)                     :: tcpu_start,tcpu_end,twall_start,&
                                  twall_end
  character(len=60)            :: mrcifile
                               
  type(dhtbl)                  :: h
  integer(ib)                  :: key(n_int,2)
  integer(is)                  :: initial_size,ntot

!----------------------------------------------------------------------
! Check the excitation class and determinant selector types
!----------------------------------------------------------------------
  ! Excitation class
  if (class_type.ne.'single'.and.class_type.ne.'single+double') then
     errmsg='Error in generate_mrci_dets: unknown excitation class '&
          //'type: '//trim(class_type)
     call error_control
  endif

  ! Determinant selector
  if (selection_type.ne.'dftmrci'.and.selection_type.ne.'abinitio') then
     errmsg='Error in generate_mrci_dets: unknown determinant '&
          //'selector type: '//trim(selection_type)
     call error_control
  endif

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)
  
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  ! Indices of the occupied and unoccupied orbitals in the reference
  ! space alpha and beta bit strings
  allocate(occ(nmo,2,nref))
  allocate(unocc(nmo,2,nref))
  occ=0
  unocc=0

!----------------------------------------------------------------------
! Get the indices of the occupied orbitals in the alpha and beta bit
! strings
!----------------------------------------------------------------------
  ! Loop over reference determinants
  do i=1,nref
     
     ! Get the indices of the occupied orbitals
     call get_occ(dref(:,:,i),occ(:,:,i))
     
     ! Get the indices of the unoccupied orbitals
     call get_unocc(dref(:,:,i),unocc(:,:,i))
     
  enddo

!----------------------------------------------------------------------
! Determine the number of occupied and unoccupied alpha and beta
! orbitals
!----------------------------------------------------------------------
  ! Number of occupied alpha and beta orbitals
  nea=n_alpha_electrons(dref(:,:,1))
  neb=n_beta_electrons(dref(:,:,1))

  ! Number of unoccupied alpha and beta orbitals
  nha=nmo-nea
  nhb=nmo-neb
  
!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------  
  initial_size=1000
  call h%initialise_table(initial_size)

!----------------------------------------------------------------------
! Insert the reference space determinants into the hash table
!----------------------------------------------------------------------
  do iref=1,nref
     call h%insert_key(dref(:,:,iref))
  enddo

!----------------------------------------------------------------------
! Generate all the singly-excited determinants
!----------------------------------------------------------------------
  ! Initialise the determinant counter
  iexc1=0

  !
  ! Alpha spin excitations
  !
  ! Loop over reference space determinants
  do iref=1,nref

     ! Loop over annihilators
     do ia1=1,nea

        ! Annihilator index
        ia=occ(ia1,1,iref)

        ! Block that the annihilated electron sits in
        ika=(ia-1)/64+1
        
        ! Loop over creators
        do ic1=1,nha

           ! Creator index
           ic=unocc(ic1,1,iref)

           ! Block that the created electron sits in
           ikc=(ic-1)/64+1
           
           ! Increment the determinant counter
           iexc1=iexc1+1
           
           ! Form the singly-excited determinant
           key=dref(:,:,iref)
           key(ika,1)=ibclr(key(ika,1),ia-1)
           key(ikc,1)=ibset(key(ikc,1),ic-1)

           ! Insert the singly-excited determinant into the hash table
           call h%insert_key(key)
           
        enddo
        
     enddo
     
  enddo
  
  !
  ! Beta spin excitations
  !
  ! Loop over reference space determinants
  do iref=1,nref

     ! Loop over annihilators
     do ia1=1,neb

        ! Annihilator index
        ia=occ(ia1,2,iref)

        ! Block that the annihilated electron sits in
        ika=(ia-1)/64+1
        
        ! Loop over creators
        do ic1=1,nhb

           ! Creator index
           ic=unocc(ic1,2,iref)

           ! Block that the created electron sits in
           ikc=(ic-1)/64+1
           
           ! Increment the determinant counter
           iexc1=iexc1+1
           
           ! Form the singly-excited determinant
           key=dref(:,:,iref)
           key(ika,2)=ibclr(key(ika,2),ia-1)
           key(ikc,2)=ibset(key(ikc,2),ic-1)

           ! Insert the singly-excited determinant into the hash table
           call h%insert_key(key)
           
        enddo
        
     enddo
     
  enddo

!----------------------------------------------------------------------
! Generate all the doubly-excited determinants
!----------------------------------------------------------------------
  if (class_type=='single+double') then

     ! Initialise the determinant counter
     iexc2=0

     !
     ! alpha alpha -> alpha alpha excitations
     !
     ! Loop over reference space determinants
     do iref=1,nref
        
        ! Loop over the first annihilator
        do ia1=1,nea-1
           
           ! First annihilator index
           ia=occ(ia1,1,iref)
           
           ! Block index for the first annihilator
           ika=(ia-1)/64+1
           
           ! Loop over the second annihilator
           do iap1=ia1+1,nea

              ! Second annihilator index
              iap=occ(iap1,1,iref)

              ! Block index for the first annihilator
              ikap=(iap-1)/64+1
           
              ! Loop over the first creator
              do ic1=1,nha-1
                 
                 ! First creator index
                 ic=unocc(ic1,1,iref)
                 
                 ! Block index for the first creator
                 ikc=(ic-1)/64+1
                 
                 ! Loop over the second creator
                 do icp1=ic1+1,nha
                    
                    ! Second creator index
                    icp=unocc(icp1,1,iref)
                    
                    ! Block index for the second creator
                    ikcp=(icp-1)/64+1
                    
                    ! Increment the determinant counter
                    iexc2=iexc2+1
                    
                    ! Form the doubly-excited determinant
                    key=dref(:,:,iref)
                    key(ika,1)=ibclr(key(ika,1),ia-1)
                    key(ikap,1)=ibclr(key(ikap,1),iap-1)
                    key(ikc,1)=ibset(key(ikc,1),ic-1)
                    key(ikcp,1)=ibset(key(ikcp,1),icp-1)
                    
                    ! Insert the doubly-excited determinant into the
                    ! hash table
                    call h%insert_key(key)
                    
                 enddo
                 
              enddo
              
           enddo
           
        enddo
        
     enddo
     
     !
     ! beta beta -> beta beta excitations
     !
     ! Loop over reference space determinants
     do iref=1,nref
        
        ! Loop over the first annihilator
        do ia1=1,neb-1
           
           ! First annihilator index
           ia=occ(ia1,2,iref)
           
           ! Block index for the first annihilator
           ika=(ia-1)/64+1
           
           ! Loop over the second annihilator
           do iap1=ia1+1,neb
              
              ! Second annihilator index
              iap=occ(iap1,2,iref)
              
              ! Block index for the first annihilator
              ikap=(iap-1)/64+1
              
              ! Loop over the first creator
              do ic1=1,nhb-1
                 
                 ! First creator index
                 ic=unocc(ic1,2,iref)
                 
                 ! Block index for the first creator
                 ikc=(ic-1)/64+1
                 
                 ! Loop over the second creator
                 do icp1=ic1+1,nhb
                    
                    ! Second creator index
                    icp=unocc(icp1,2,iref)
                    
                    ! Block index for the second creator
                    ikcp=(icp-1)/64+1
                 
                    ! Increment the determinant counter
                    iexc2=iexc2+1

                    ! Form the doubly-excited determinant
                    key=dref(:,:,iref)
                    key(ika,2)=ibclr(key(ika,2),ia-1)
                    key(ikap,2)=ibclr(key(ikap,2),iap-1)
                    key(ikc,2)=ibset(key(ikc,2),ic-1)
                    key(ikcp,2)=ibset(key(ikcp,2),icp-1)

                    ! Insert the doubly-excited determinant into the
                    ! hash table
                    call h%insert_key(key)
                    
                 enddo
                 
              enddo
              
           enddo
           
        enddo
        
     enddo
     
     !
     ! alpha beta -> alpha beta excitations
     !
     ! Loop over reference space determinants
     do iref=1,nref
        
        ! Loop over the first annihilator (alpha)
        do ia1=1,nea
           
           ! First annihilator index
           ia=occ(ia1,1,iref)
           
           ! Block index for the first annihilator
           ika=(ia-1)/64+1
        
           ! Loop over the second annihilator (beta)
           do iap1=1,neb

              ! Second annihilator index
              iap=occ(iap1,2,iref)

              ! Block index for the first annihilator
              ikap=(iap-1)/64+1
              
              ! Loop over the first creator (alpha)
              do ic1=1,nha
                 
                 ! First creator index
                 ic=unocc(ic1,1,iref)

                 ! Block index for the first creator
                 ikc=(ic-1)/64+1
              
                 ! Loop over the second creator (beta)
                 do icp1=1,nhb

                    ! Second creator index
                    icp=unocc(icp1,2,iref)

                    ! Block index for the second creator
                    ikcp=(icp-1)/64+1
                 
                    ! Increment the determinant counter
                    iexc2=iexc2+1
                    
                    ! Form the doubly-excited determinant
                    key=dref(:,:,iref)
                    key(ika,1)=ibclr(key(ika,1),ia-1)    ! alpha annihilator
                    key(ikap,2)=ibclr(key(ikap,2),iap-1) ! beta annihilator
                    key(ikc,1)=ibset(key(ikc,1),ic-1)    ! alpha creator
                    key(ikcp,2)=ibset(key(ikcp,2),icp-1) ! beta creator

                    ! Insert the doubly-excited determinant into the
                    ! hash table
                    call h%insert_key(key)
                    
                 enddo
                 
              enddo
              
           enddo
           
        enddo
        
     enddo

  endif
     
!----------------------------------------------------------------------
! Number of unique singly-excited determinants 
!----------------------------------------------------------------------
  ndet=h%n_keys_stored
  
!----------------------------------------------------------------------
! Write the unique singly-excited determinants to disk
!----------------------------------------------------------------------
  call scratch_name('detmrci',mrcifile)
  call h%dump_keys(mrcifile)

!----------------------------------------------------------------------
! Delete the hash table
!----------------------------------------------------------------------
  call h%delete_table
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(occ,unocc)
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'generate_mrci_excitations')
  
  return
  
end subroutine generate_mrci_dets
  
