!**********************************************************************
! Routines for the SCI diagonalisation of the MRCI reference space
! Hamiltonian
!**********************************************************************

module ref_sci

contains
  
!######################################################################
! ref_diag_mrci_sci: Diagonalisation of the reference space Hamiltonian
!                    using an SCI algorithm
!######################################################################
#ifdef CBINDING
  subroutine ref_diag_mrci_sci(irrep,nroots,confscr,nconf,vecscr) &
       bind(c,name="ref_diag_mrci_sci")
#else
  subroutine ref_diag_mrci_sci(irrep,nroots,confscr,nconf,vecscr)
#endif
      
    use constants
    use bitglobal
    use hbuild_double
    use iomod
  
    implicit none

    ! Irrep number and the requested number of roots
    integer(is), intent(in)    :: irrep
    integer(is), intent(inout) :: nroots
    
    ! Array of reference configuration scratch file numbers
    integer(is), intent(in)    :: confscr(0:nirrep-1)

    ! Array of numbers of referecne configurations
    integer(is), intent(in)    :: nconf(0:nirrep-1)

    ! Eigenvector scratch file index
    integer(is), intent(out)   :: vecscr

    ! Number of configurations found in the scratch file
    integer(is)                :: nconf1

    ! Number of (n_bits)-bit integers required to represent
    ! the configurations
    integer(is)                :: n_int_I
      
    ! Number of internal and external MOs
    integer(is)                :: nmoI,nmoE

    ! Reference configurations and SOPs
    integer(ib), allocatable   :: conf(:,:,:),sop(:,:,:)
  
    ! MO mapping arrays
    integer(is)                :: m2c(nmo),c2m(nmo)
      
    ! Numbers of CSFs
    integer(is)                :: csfdim
    integer(is), allocatable   :: offset(:)
      
    ! Dimension of the CSF basis past which full diagonalisation
    ! will not be used
    integer(is), parameter     :: full_lim=1200
    
    ! On-diagonal Hamiltonian matrix elements
    real(dp), allocatable      :: hii(:)
    
    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    real(dp), allocatable      :: averageii(:)

    ! Hamiltonian scratch file number
    integer(is)                :: hscr
    
    ! Number of Hamiltonian scratch file records
    integer(is)                :: nrec

    ! On-diagonal Hamiltonian matrix elements for the P and Q
    ! subspaces
    integer(is)                :: csfdimP
    real(dp), allocatable      :: hiiP(:),averageiiP(:)
    real(dp), allocatable      :: hiiQ(:)
    
    ! P and Q space information
    integer(is)                :: nP,nQ,csfdimQ
    integer(is), allocatable   :: iP(:),iQ(:)
    integer(is), allocatable   :: offsetP(:),offsetQ(:)
    integer(is), allocatable   :: confmap(:)
    
    ! Eigenpairs of the projected reference space Hamiltonian
    real(dp), allocatable      :: vecP(:,:),EP(:)

    ! 1st-order corrected wave functions
    real(dp), allocatable      :: Avec(:,:)

    ! 2nd-order energy corrections
    real(dp), allocatable      :: E2(:)

    ! P space weights
    real(dp), allocatable      :: WP(:)
    
    ! Temporary hard-wiring of the maximum number of iterations
    integer(is), parameter     :: maxiter=15

    ! Temporary hard-wiring of the P space weight convergence
    ! threshold
    real(dp), parameter        :: WP_thrsh=0.95_dp
    
    ! Work arrays
    integer(is)                :: harr2dim
    real(dp), allocatable      :: harr2(:)
    
    ! Everything else
    integer(is)                :: i,iroot
    logical                    :: converged
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    harr2dim=maxval(ncsfs(0:nomax))**2
    allocate(harr2(harr2dim))
    harr2=0.0d0
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    if (verbose) then
       write(6,'(/,72a)') ('-',i=1,52)
       write(6,'(2(x,a),/,x,a)') &
            'Reference space SCI diagonalisation in the',&
            trim(irreplbl(irrep,ipg)),'subspace'
       write(6,'(72a)') ('-',i=1,52)
    endif

!----------------------------------------------------------------------
! Return if there are no configurations for the current irrep
!----------------------------------------------------------------------
    if (nconf(irrep) == 0) then
       if (verbose) &
            write(6,'(/,x,a)') 'No reference space configurations of '&
            //trim(irreplbl(irrep,ipg))//' symmetry'
       nroots=0
       return
    endif

!----------------------------------------------------------------------
! Read the configurations from file
!----------------------------------------------------------------------
    call read_ref_confs(confscr(irrep),nconf1,n_int_I,nmoI,nmoE,conf,&
         sop,m2c,c2m)

!----------------------------------------------------------------------
! Sanity check on the number of configurations
!----------------------------------------------------------------------
    if (nconf(irrep) /= nconf1) then
       errmsg='Error in ref_diag_mrci: '&
            //'inconsistent configuration numbers'
       call error_control
    endif
  
!----------------------------------------------------------------------
! Determine the total number of CSFs and the offsets for each
! configuration
!----------------------------------------------------------------------
    allocate(offset(nconf1+1))
    offset=0

    call basis_dimensions(csfdim,offset,sop,n_int_I,nconf1)

!----------------------------------------------------------------------
! Sanity check on the number of roots
!----------------------------------------------------------------------
    if (csfdim < nroots) then
       errmsg='Error in ref_diag_mrci: N_roots > N_CSF'
       call error_control
    endif

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(hii(csfdim))
    allocate(averageii(nconf1))
    allocate(E2(nroots))
    allocate(Avec(csfdim,nroots))
    allocate(confmap(csfdim))
    allocate(iP(nconf1))
    allocate(iQ(nconf1))
    allocate(hiiP(csfdim))
    allocate(hiiQ(csfdim))
    allocate(averageiiP(nconf1))
    allocate(offsetP(nconf1+1))
    allocate(offsetQ(nconf1+1))
    allocate(EP(nroots))
    allocate(WP(nroots))
    hii=0.0d0
    averageii=0.0d0
    E2=0.0d0
    Avec=0.0d0
    confmap=0
    iP=0
    iQ=0
    hiiP=0.0d0
    hiiQ=0.0d0
    averageiiP=0.0d0
    offsetP=0
    offsetQ=0
    EP=0.0d0
    WP=0.0d0
    
!----------------------------------------------------------------------
! Compute the on-diagonal Hamiltonian matrix elements and the
! their spin-coupling averaged values
!----------------------------------------------------------------------
    call hii_double(nconf1,csfdim,offset,conf,sop,n_int_I,m2c,irrep,&
         hii,averageii)

!----------------------------------------------------------------------
! Initialialise the P and Q spaces
!----------------------------------------------------------------------
    ! Determine the P and Q spaces
    call partition_confs(nroots,nconf1,offset,csfdim,hii,nP,iP,nQ,iQ,&
         offsetP,offsetQ,confmap)

    ! Fill in the P and Q space on-diagonal Hamiltonian matrix
    ! elements
    call partition_hii(csfdim,hii,csfdimP,csfdimQ,hiiP,hiiQ,nconf1,&
         offset,nP,iP,nQ,iQ,averageii,averageiiP)
    
!----------------------------------------------------------------------
! Iterative improvement of the P space
!----------------------------------------------------------------------
    ! Table header
    if (verbose) then
       write(6,'(/,x,29a)') ('*', i=1,29)
       write(6,'(2x,a)') 'Iteration  Nconf   min W(P)'
       write(6,'(x,29a)') ('*', i=1,29)
    endif
       
    ! Perform the iterations
    do i=1,maxiter

       ! Allocate the P space eigenvector array
       if (allocated(vecP)) deallocate(vecP)
       allocate(vecP(csfdimP,nroots))
       
       ! Diagonalise the P space Hamiltonian
       call diag_pspace(nconf1,nP,csfdim,csfdimP,nroots,iP,offset,&
            offsetP,averageii,hiiP,conf,sop,n_int_I,m2c,irrep,&
            vecP,EP,vecscr)

       ! Compute the ENPT2 wave function and energy corrections
       call pt2_corrections(csfdim,csfdimP,csfdimQ,nroots,vecP,EP,&
            hiiQ,Avec,E2,n_int_I,nconf1,conf,sop,nP,nQ,iP,iQ,offset,&
            averageii,harr2dim,harr2,m2c,offsetP,offsetQ,WP)

       ! Update the P and Q spaces
       call update_partitioning(nconf1,csfdim,csfdimP,csfdimQ,nroots,&
            Avec,confmap,nP,nQ,iP,iQ,offset,offsetP,offsetQ)

       ! Fill in the P and Q space on-diagonal Hamiltonian matrix
       ! elements
       call partition_hii(csfdim,hii,csfdimP,csfdimQ,hiiP,hiiQ,nconf1,&
            offset,nP,iP,nQ,iQ,averageii,averageiiP)
       
       ! Output our progress
       if (verbose) write(6,'(4x,i5,4x,i5,4x,F6.4)') i,nP,minval(WP)
       
       ! Exit if we have reached convergence
       converged=check_conv(nroots,WP,WP_thrsh)
       if (converged) exit
       
    enddo

    ! Table footer
    if (verbose) write(6,'(x,29a)') ('*', i=1,29)

!----------------------------------------------------------------------
! Prune the reference space
!----------------------------------------------------------------------
    
    STOP
    
    return
  
  end subroutine ref_diag_mrci_sci

!######################################################################

  subroutine partition_confs(nroots,nconf,offset,csfdim,hii,nP,iP,nQ,iQ,&
       offsetP,offsetQ,confmap)

    use constants
    use bitglobal
    use utils
  
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: nroots,nconf
  
    ! CSF offsets
    integer(is), intent(in)  :: offset(nconf+1)

    ! On-diagonal Hamiltonian matrix elements
    integer(is), intent(in)  :: csfdim
    real(dp), intent(in)     :: hii(csfdim)

    ! P and Q space configurations
    integer(is), intent(out) :: nP,nQ
    integer(is), intent(out) :: iP(nconf),iQ(nconf)
    integer(is), intent(out) :: offsetP(nconf+1),offsetQ(nconf+1)
    integer(is), intent(out) :: confmap(csfdim)
    
    ! Sorting arrays
    integer(is), allocatable :: indx(:)
    
    ! Everything else
    integer(is)              :: iconf,icsf,i
    integer(is)              :: countP,countQ
    integer(is)              :: totalP,totalQ
    integer(is)              :: nlow
    integer(is), allocatable :: isurvive(:)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(isurvive(nconf))
    isurvive=0

    allocate(indx(csfdim))
    indx=0
    
!----------------------------------------------------------------------
! Sort the on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
    call dsortindxa1('A',csfdim,hii,indx)

!----------------------------------------------------------------------  
! Determine the configurations from which each CSF is generated
!----------------------------------------------------------------------
    ! Loop over configurations
    do iconf=1,nconf
       
       ! Loop over the CSFs generated by this configuration
       do icsf=offset(iconf),offset(iconf+1)-1

          ! Fill in the CSF-to-conf mapping array
          confmap(icsf)=iconf
          
       enddo
       
    enddo

!----------------------------------------------------------------------  
! Determine the configurations corresponding to the P space CSFs
!----------------------------------------------------------------------
    ! Initialisation
    isurvive=0

    ! For now we shall hardwire the number of initially selected CSFs
    nlow=2*nroots
        
    ! Loop over the lowest energy CSFs
    do i=1,nlow

       ! i'th lowest energy CSF
       icsf=indx(i)

       ! Flag the configuration that generates this CSF for
       ! survival
       iconf=confmap(icsf)
       isurvive(iconf)=1
       
    enddo

!----------------------------------------------------------------------
! Number of P and Q space configurations
!----------------------------------------------------------------------
    nP=sum(isurvive)
    nQ=nconf-nP
    
!----------------------------------------------------------------------
! Fill in the arrays of P and Q space configurations
!----------------------------------------------------------------------
    countP=0
    countQ=0
  
    do i=1,nconf
       if (isurvive(i) == 1) then
          countP=countP+1
          iP(countP)=i
       else
          countQ=countQ+1
          iQ(countQ)=i
       endif
    enddo
    
!----------------------------------------------------------------------
! Fill in the CSF offsets for the P and Q space configurations
!----------------------------------------------------------------------
    totalP=1
    totalQ=1
    countP=0
    countQ=0
        
    do iconf=1,nconf
       if (isurvive(iconf) == 1) then
          countP=countP+1
          offsetP(countP)=totalP
          totalP=totalP+offset(iconf+1)-offset(iconf)
       else
          countQ=countQ+1
          offsetQ(countQ)=totalQ
          totalQ=totalQ+offset(iconf+1)-offset(iconf)
       endif
    enddo

    offsetP(nP+1)=totalP
    offsetQ(nQ+1)=totalQ
    
    return
  
  end subroutine partition_confs

!######################################################################

  subroutine partition_hii(csfdim,hii,csfdimP,csfdimQ,hiiP,hiiQ,&
       nconf,offset,nP,iP,nQ,iQ,averageii,averageiiP)

    use constants
    use bitglobal

    implicit none
    
    ! On-diagonal Hamiltonian matrix elements for the full space
    integer(is), intent(in)  :: csfdim
    real(dp), intent(in)     :: hii(csfdim)
 
    ! On-diagonal Hamiltonian matrix elements for the P and Q
    ! spaces
    integer(is), intent(out) :: csfdimP,csfdimQ
    real(dp), intent(out)    :: hiiP(csfdim),hiiQ(csfdim)

    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: averageii(nconf)

    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    ! for the P space
    real(dp), intent(out)    :: averageiiP(nconf+1)
    
    ! CSF offsets
    integer(is), intent(in)  :: nconf
    integer(is), intent(in)  :: offset(nconf+1)

    ! P and Q space configurations
    integer(is), intent(in)  :: nP,nQ
    integer(is), intent(in)  :: iP(nconf),iQ(nconf)
    
    ! Everything else
    integer(is)              :: i,iconf,icsf,ncsf,count

!----------------------------------------------------------------------
! Number of P and Q space CSFs
!----------------------------------------------------------------------
    ! Initialisation
    csfdimP=0

    ! Loop over the P space configurations
    do i=1,nP
       iconf=iP(i)

       ! Number of CSFs generated by this configuration
       ncsf=offset(iconf+1)-offset(iconf)

       ! Running total number of P space CSFs
       csfdimP=csfdimP+ncsf
       
    enddo

    ! Number of Q space CSFs
    csfdimQ=csfdim-csfdimP

!----------------------------------------------------------------------    
! Fill in the P space on-diagonal matrix element array
!----------------------------------------------------------------------
    ! P space CSF counter
    count=0
    
    ! Loop over P space configurations
    do i=1,nP
       iconf=iP(i)

       ! Loop over the CSFs generated by this configuration
       do icsf=offset(iconf),offset(iconf+1)-1

          ! Update the P space CSF counter
          count=count+1

          ! Save the on-diagonal Hamiltonian matrix element
          ! for this CSF
          hiiP(count)=hii(icsf)
          
       enddo

       ! Save the spin-coupling averaged on-diagonal Hamiltonian
       ! matrix element value for this configuration
       averageiiP(i)=averageii(iconf)
              
    enddo

!----------------------------------------------------------------------    
! Fill in the Q space on-diagonal matrix element array
!----------------------------------------------------------------------
    ! Q space CSF counter
    count=0

    ! Loop over Q space configurations
    do i=1,nQ
       iconf=iQ(i)

       ! Loop over the CSFs generated by this configuration
       do icsf=offset(iconf),offset(iconf+1)-1

          ! Update the Q space CSF counter
          count=count+1

          ! Save the on-diagonal Hamiltonian matrix element
          ! for this CSF
          hiiQ(count)=hii(icsf)
          
       enddo
       
    enddo
    
    return
    
  end subroutine partition_hii

!######################################################################  

  subroutine diag_pspace(nconf,nP,csfdim,csfdimP,nroots,iP,offset,&
       offsetP,averageii,hiiP,conf,sop,n_int_I,m2c,irrep,vecP,EP,&
       vecscr)

    use constants
    use bitglobal
    use hbuild_double
    use full_diag
    use iomod
    
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: nconf,nP,csfdim,csfdimP,nroots

    ! P space configurations
    integer(is), intent(in)  :: iP(nconf)

    ! CSF offsets
    integer(is), intent(in)  :: offset(nconf+1)
    integer(is), intent(in)  :: offsetP(nconf+1)

    ! P space on-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: hiiP(csfdim)
    
    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: averageii(nconf)

    ! Configurations and SOPs
    integer(is), intent(in)  :: n_int_I
    integer(ib), intent(in)  :: conf(n_int_I,2,nconf)
    integer(ib), intent(in)  :: sop(n_int_I,2,nconf)

    ! MO mapping array
    integer(is), intent(in)  :: m2c(nmo)

    ! Irrep number
    integer(is), intent(in)  :: irrep

    ! P space eigenpairs
    real(dp), intent(out)    :: vecP(csfdimP,nroots),EP(nroots)
    
    ! P space eigenvector scratch file number
    integer(is), intent(out) :: vecscr

    ! Everything else
    integer(is)              :: hscr,nrec
    logical                  :: verbose_save
    
!----------------------------------------------------------------------
! Make the output non-verbose
!----------------------------------------------------------------------
    verbose_save=verbose
    verbose=.false.
    
!----------------------------------------------------------------------
! Save to disk the non-zero P space off-diagonal Hamiltonian matrix
! elements
!----------------------------------------------------------------------
    call save_hij_double_selected(nP,iP,offsetP,nconf,csfdim,offset,&
         averageii,conf,sop,n_int_I,m2c,irrep,hscr,nrec,'hij_ref')

!----------------------------------------------------------------------
! Diagonalise the P space Hamiltonian
!----------------------------------------------------------------------
    call diag_full(hscr,nrec,csfdimP,hiiP,irrep,nroots,vecscr,&
         'refvec')

!----------------------------------------------------------------------
! Retrieve the P space eigenpairs
!----------------------------------------------------------------------
    ! Read in the eigenpairs
    call read_all_eigenpairs(vecscr,vecP,EP,csfdimP,nroots)

    ! Subtract off E_SCF from the energies to get the true eigenvalues
    EP=EP-escf

!----------------------------------------------------------------------
! Restore the original verbose flag value
!----------------------------------------------------------------------
    verbose=verbose_save
    
    return
    
  end subroutine diag_pspace
    
!######################################################################
! pt2_corrections: Calculation of the ENPT2 corrections
!######################################################################
  subroutine pt2_corrections(csfdim,csfdimP,csfdimQ,nroots,vecP,EP,&
       hiiQ,Avec,E2,n_int_I,nconf,conf,sop,nP,nQ,iP,iQ,offset,&
       averageii,harr2dim,harr2,m2c,offsetP,offsetQ,WP)

    use constants
    use bitglobal
    use conftype
    use mrci_integrals
    use mrciutils
    use hbuild_mrci
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: csfdim,csfdimP,csfdimQ,nroots

    ! Eigenpairs of the reference space Hamiltonian projected
    ! onto the P space
    real(dp), intent(in)    :: vecP(csfdimP,nroots)
    real(dp), intent(in)    :: EP(nroots)

    ! On-diagonal Hamiltonian matrix elements for the Q space
    ! CSFs
    real(dp), intent(in)    :: hiiQ(csfdim)
    
    ! 1st-order corrected wave functions
    real(dp), intent(out)   :: Avec(csfdim,nroots)

    ! 2nd-order energy corrections
    real(dp), intent(out)   :: E2(nroots)
    
    ! Configurations and SOPs
    integer(is), intent(in) :: nconf,n_int_I
    integer(ib), intent(in) :: conf(n_int_I,2,nconf)
    integer(ib), intent(in) :: sop(n_int_I,2,nconf)

    ! P and Q space configurations
    integer(is), intent(in) :: nP,nQ
    integer(is), intent(in) :: iP(nconf),iQ(nconf)

    ! CSF offsets
    integer(is), intent(in) :: offset(nconf+1)

    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: averageii(nconf)

    ! Work arrays
    integer(is), intent(in) :: harr2dim
    real(dp)                :: harr2(harr2dim)

    ! MO mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! CSF offsets for the P and Q space configurations
    integer(is), intent(in) :: offsetP(nconf+1)
    integer(is), intent(in) :: offsetQ(nconf+1)

    ! P space weights
    real(dp), intent(out)   :: WP(nroots)
    
    ! Difference configuration information
    integer(is)             :: ndiff
    integer(is)             :: Dw(nmo,2)
    
    ! Number of open shells preceding each MO
    integer(is)             :: nbefore(nmo)

    ! Indices of the singly-occupied ket MOs
    integer(is)             :: socc(nmo)
    integer(is)             :: nsocc
    
    ! Working arrays
    integer(is), parameter  :: maxexci=2
    integer(is)             :: hlist(maxexci),plist(maxexci)
    integer(ib)             :: kconf_full(n_int,2),bconf_full(n_int,2)
    integer(ib)             :: ksop_full(n_int,2),bsop_full(n_int,2)
    
    ! Everything else
    integer(is)             :: iket,ibra,ikcsf,ibcsf,icsf
    integer(is)             :: kconf,bconf
    integer(is)             :: knopen,knsp,bnopen,bnsp
    integer(is)             :: nexci
    integer(is)             :: iroot,counter
    real(dp)                :: ediff,norm

!----------------------------------------------------------------------
! Initialise the A-vectors to the P-space eigenvectors
!----------------------------------------------------------------------
    Avec=0.0d0

    do iroot=1,nroots
       Avec(1:csfdimP,iroot)=vecP(:,iroot)
    enddo
    
!----------------------------------------------------------------------
! Compute the matrix elements <w omega|H|Psi_I^0>
!----------------------------------------------------------------------
    ! Loop over the ket (P space) configurations
    do iket=1,nP
       kconf=iP(iket)
       
       ! Number of open shells in the ket configuration
       knopen=sop_nopen(sop(:,:,kconf),n_int_I)

       ! Number of ket CSFs
       knsp=ncsfs(knopen)

       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_I,:)=conf(:,:,kconf)
       ksop_full(1:n_int_I,:)=sop(:,:,kconf)

       ! Package the ket configuration information
       call package_confinfo_offdiag(ksop_full,kconf_full,socc,nsocc,&
            Dw,ndiff,nbefore)

       ! Loop over the bra (Q space) configurations
       do ibra=1,nQ
          bconf=iQ(ibra)

          ! Compute the excitation degree between the two
          ! configurations
          nexci=exc_degree_conf(conf(:,:,kconf),conf(:,:,bconf),&
               n_int_I)

          ! Cycle if the excitation degree is greater than 2
          if (nexci > 2) cycle

          ! Number of open shells in the bra configuration
          bnopen=sop_nopen(sop(:,:,bconf),n_int_I)

          ! Number of bra CSFs
          bnsp=ncsfs(bnopen)

          ! Bra configuration and SOP in the full MO space
          bconf_full=0_ib
          bsop_full=0_ib
          bconf_full(1:n_int_I,:)=conf(:,:,bconf)
          bsop_full(1:n_int_I,:)=sop(:,:,bconf)

          ! Get the indices of the MOs involved in the excitation
          hlist=0
          plist=0
          call get_exci_indices(conf(:,:,kconf),conf(:,:,bconf),&
               n_int_I,hlist(1:nexci),plist(1:nexci),nexci)

          ! Compute the matrix elements between the CSFs generated
          ! by the bra and ket configurations
          call hij_mrci(harr2,harr2dim,nexci,bconf,kconf,&
               bsop_full,ksop_full,bnsp,knsp,bnopen,knopen,&
               hlist,plist,m2c,socc,nsocc,nbefore,Dw,ndiff,&
               offset,offset,nconf+1,nconf+1,averageii(bconf),&
               averageii(kconf))
          
          ! Loop over roots
          do iroot=1,nroots

             counter=0

             ! Loop over the ket (P space) CSFs
             do ikcsf=offsetP(iket),offsetP(iket+1)-1

                ! Cycle if the ket CSF coefficient is tiny
                if (abs(vecP(ikcsf,iroot)) < epshij) cycle
                
                ! Loop over the bra (Q space) CSFs
                do ibcsf=offsetQ(ibra),offsetQ(ibra+1)-1
                   counter=counter+1

                   Avec(csfdimP+ibcsf,iroot)=&
                        Avec(csfdimP+ibcsf,iroot)&
                        +harr2(counter)*vecP(ikcsf,iroot)
                   
                enddo
                   
             enddo
                
          enddo
                    
       enddo
       
    enddo

!----------------------------------------------------------------------
! Apply the denominators and compute the energy corrections
!----------------------------------------------------------------------
    ! Initialisation
    E2=0.0d0

    ! Loop over roots
    do iroot=1,nroots

       ! Loop over the Q space CSFs
       do icsf=1,csfdimQ

          ! E^(0) - H_ii
          ediff=EP(iroot)-hiiQ(icsf)

          ! Energy correction
          E2(iroot)=E2(iroot)+Avec(csfdimP+icsf,iroot)**2/ediff

          ! A-vector element
          Avec(csfdimP+icsf,iroot)=Avec(csfdimP+icsf,iroot)/ediff

       enddo

    enddo

!----------------------------------------------------------------------
! Normalisation of the 1st-order wave function corrections
!----------------------------------------------------------------------
    ! Loop over roots
    do iroot=1,nroots

       ! Norm of the 1st-order corrected wave function
       norm=dot_product(Avec(:,iroot),Avec(:,iroot))
       norm=sqrt(norm)

       ! P space weight
       WP(iroot)=1.0d0/norm
       
       ! Normalisation
       Avec(:,iroot)=Avec(:,iroot)/norm
       
    enddo
    
    return
    
  end subroutine pt2_corrections

!######################################################################

  subroutine update_partitioning(nconf,csfdim,csfdimP,csfdimQ,nroots,&
       Avec,confmap,nP,nQ,iP,iQ,offset,offsetP,offsetQ)

    use constants
    use bitglobal
    use utils
    
    implicit none

    ! Dimensions
    integer(is), intent(in)    :: nconf,csfdim,nroots
    integer(is), intent(inout) :: csfdimP,csfdimQ
    
    ! 1st-order corrected wave functions
    real(dp), intent(in)       :: Avec(csfdim,nroots)

    ! CSF-to-conf map
    integer(is), intent(in)    :: confmap(csfdim)

    ! P and Q space configurations
    integer(is), intent(inout) :: nP,nQ
    integer(is), intent(inout) :: iP(nconf),iQ(nconf)
    
    ! CSF offsets
    integer(is), intent(in)    :: offset(nconf+1)
    integer(is), intent(inout) :: offsetP(nconf+1)
    integer(is), intent(inout) :: offsetQ(nconf+1)
    
    ! Configuration selection threshold (hard-wired for now)
    real(dp), parameter        :: thrsh=0.03_dp

    ! Selected CSFs and confs
    integer(is), allocatable   :: isel_csf(:),isel_conf(:)
    
    ! Everything else
    integer(is)                :: i,iroot,icsf,icsfP,icsfQ,iconf
    integer(is)                :: countP,countQ,totalP,totalQ

!----------------------------------------------------------------------    
! Allocate arrays
!----------------------------------------------------------------------
    allocate(isel_csf(csfdim))
    allocate(isel_conf(nconf))

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    isel_csf=0
    isel_conf=0
    
!----------------------------------------------------------------------
! Determine the indices of the dominant CSFs
!----------------------------------------------------------------------
    ! Loop over roots
    do iroot=1,nroots

       ! Loop over CSFs
       do icsf=1,csfdim

          if (abs(Avec(icsf,iroot)) > thrsh) isel_csf(icsf)=1
          
       enddo
              
    enddo

!----------------------------------------------------------------------  
! Determine the P space configurations that generate dominant CSFs
!----------------------------------------------------------------------
    ! Loop over P space configurations
    do i=1,nP
       iconf=iP(i)
       
       ! Loop over the P space CSFs generated by this configuration
       do icsfP=offsetP(i),offsetP(i+1)-1
       
          ! Are we at a dominant CSF?
          if (isel_csf(icsfP) == 1) then
       
             isel_conf(iconf)=1
             
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------  
! Determine the Q space configurations that generate dominant CSFs
!----------------------------------------------------------------------
    ! Loop over Q space configurations
    do i=1,nQ
       iconf=iQ(i)

       ! Loop over the Q space CSFs generated by this configuration
       do icsfQ=offsetQ(i),offsetQ(i+1)-1
       
          ! Are we at a dominant CSF?
          if (isel_csf(csfdimP+icsfQ) == 1) then
       
             isel_conf(iconf)=1
             
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Number of P and Q space configurations
!----------------------------------------------------------------------
    nP=sum(isel_conf)
    nQ=nconf-nP

!----------------------------------------------------------------------
! Fill in the arraya of P and Q space configurations
!----------------------------------------------------------------------
    iP=0
    iQ=0

    countP=0
    countQ=0
    
    do iconf=1,nconf
       if (isel_conf(iconf) == 1) then
          countP=countP+1
          iP(countP)=iconf
       else
          countQ=countQ+1
          iQ(countQ)=iconf
       endif
    enddo

!----------------------------------------------------------------------
! Fill in the CSF offsets for the P and Q space configurations
!----------------------------------------------------------------------
    offsetP=0
    offsetQ=0

    totalP=1
    totalQ=1
    countP=0
    countQ=0
    
    do iconf=1,nconf
       if (isel_conf(iconf) == 1) then
          countP=countP+1
          offsetP(countP)=totalP
          totalP=totalP+offset(iconf+1)-offset(iconf)
       else
          countQ=countQ+1
          offsetQ(countQ)=totalQ
          totalQ=totalQ+offset(iconf+1)-offset(iconf)
       endif
    enddo

    offsetP(nP+1)=totalP
    offsetQ(nQ+1)=totalQ

    return
    
  end subroutine update_partitioning
  
!######################################################################

  function check_conv(nroots,WP,thrsh) result(converged)

    use constants

    implicit none

    ! Function result
    logical                 :: converged

    ! Number of roots
    integer(is), intent(in) :: nroots
    
    ! P space weights
    real(dp), intent(in)    :: WP(nroots)

    ! Convergence threshold
    real(dp), intent(in)    :: thrsh

    ! Everything else
    integer(is)             :: iroot
    integer(is)             :: iconv(nroots)

    iconv=0
    
    do iroot=1,nroots
       if (WP(iroot) >= thrsh) iconv(iroot)=1
    enddo

    if (sum(iconv) == nroots) then
       converged=.true.
    else
       converged=.false.
    endif
    
    return
    
  end function check_conv
  
!######################################################################
  
end module ref_sci
