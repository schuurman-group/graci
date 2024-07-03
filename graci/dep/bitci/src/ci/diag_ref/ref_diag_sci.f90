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
    use ref_guess
    use full_diag
    use utils
    use iomod
    use conftype
  
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

    ! On-diagonal Hamiltonian matrix elements for the selected
    ! and unselected subspaces
    integer(is)                :: csfdim_sel
    real(dp), allocatable      :: hii_sel(:),averageii_sel(:)
    real(dp), allocatable      :: hii_unsel(:)
    
    ! Selected subspace information
    integer(is)                :: nsel,nunsel,csfdim_unsel
    integer(is), allocatable   :: isel(:),iunsel(:)
    integer(is), allocatable   :: offset_sel(:),offset_unsel(:)
    
    ! Eigenpairs of the projected reference space Hamiltonian
    real(dp), allocatable      :: vec_sel(:,:),ener_sel(:)

    ! A-vectors
    real(dp), allocatable      :: Avec(:,:)

    ! 2nd-order energy corrections
    real(dp), allocatable      :: E2(:)
    
    ! Temporary hard-wiring of the maximum number of iterations
    integer(is), parameter     :: maxiter=15

    ! Work arrays
    integer(is)                :: harr2dim
    real(dp), allocatable      :: harr2(:)
    
    ! Everything else
    integer(is)                :: i

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
    hii=0.0d0
    averageii=0.0d0
    E2=0.0d0

!----------------------------------------------------------------------
! Compute the on-diagonal Hamiltonian matrix elements and the
! their spin-coupling averaged values
!----------------------------------------------------------------------
    call hii_double(nconf1,csfdim,offset,conf,sop,n_int_I,m2c,irrep,&
         hii,averageii)

!----------------------------------------------------------------------
! Determine the initial subspace
!----------------------------------------------------------------------
    call init_subspace(nroots,nconf1,offset,csfdim,hii,nsel,isel,&
         nunsel,iunsel,offset_sel,offset_unsel)

!----------------------------------------------------------------------
! Diagonalisation of the reference space Hamiltonian projected onto
! the initial subspace
!----------------------------------------------------------------------
    ! Fill in the selected on-diagonal Hamiltonian matrix elements
    call hii_selected(csfdim,hii,csfdim_sel,csfdim_unsel,hii_sel,&
         hii_unsel,nconf1,offset,nsel,isel,nunsel,iunsel,averageii,&
         averageii_sel)
    
    ! Save to disk the non-zero selected off-diagonal Hamiltonian
    ! matrix elements
    call save_hij_double_selected(nsel,isel,offset_sel,nconf1,&
         csfdim,offset,averageii,conf,sop,n_int_I,m2c,irrep,hscr,&
         nrec,'hij_ref')

    ! Diagonalise the projected reference space Hamiltonian
    call diag_full(hscr,nrec,csfdim_sel,hii_sel,irrep,nroots,vecscr,&
         'refvec')

!----------------------------------------------------------------------
! Retrieve the eigenpairs of the projected reference space Hamiltonian
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(vec_sel(csfdim_sel,nroots))
    allocate(ener_sel(nroots))
    vec_sel=0.0d0
    ener_sel=0.0d0

    ! Read in the eigenpairs
    call read_all_eigenpairs(vecscr,vec_sel,ener_sel,csfdim_sel,nroots)

    ! Subtract off E_SCF from the energies to get the true eigenvalues
    ener_sel=ener_sel-escf
    
!----------------------------------------------------------------------
! Iterative improvement of the subspace
!----------------------------------------------------------------------
    do i=1,maxiter

       ! Allocate arrays
       if (allocated(Avec)) deallocate(Avec)
       allocate(Avec(csfdim_unsel,nroots))
       Avec=0.0d0
       
       ! Compute the ENPT2 wave function and energy corrections
       call pt2_corrections(csfdim_sel,csfdim_unsel,nroots,vec_sel,&
            ener_sel,hii_unsel,Avec,E2,n_int_I,nconf1,conf,sop,nsel,&
            nunsel,isel,iunsel,offset,averageii,harr2dim,harr2,m2c,&
            offset_sel,offset_unsel)

       ! Find the most important unselected configurations
       
       STOP
       
    enddo
    

    STOP

    
    return
  
  end subroutine ref_diag_mrci_sci

!######################################################################

  subroutine init_subspace(nroots,nconf,offset,ncsf,hii,nsel,isel,&
       nunsel,iunsel,offset_sel,offset_unsel)

    use constants
    use bitglobal
    use utils
  
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: nroots,nconf
  
    ! CSF offsets
    integer(is), intent(in)  :: offset(nconf+1)

    ! On-diagonal Hamiltonian matrix elements
    integer(is), intent(in)  :: ncsf
    real(dp), intent(in)     :: hii(ncsf)

    ! Selected configurations
    integer(is), intent(out) :: nsel,nunsel
    integer(is), allocatable :: isel(:),iunsel(:)
    integer(is), allocatable :: offset_sel(:),offset_unsel(:)
    
    ! Sorting arrays
    integer(is), allocatable :: indx(:)
    
    ! Everything else
    integer(is)              :: iconf,icsf,i
    integer(is)              :: count_sel,count_unsel
    integer(is)              :: total_sel,total_unsel
    integer(is)              :: nlow
    integer(is), allocatable :: confmap(:)
    integer(is), allocatable :: isurvive(:)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(confmap(ncsf))
    confmap=0

    allocate(isurvive(nconf))
    isurvive=0

    allocate(indx(ncsf))
    indx=0
    
!----------------------------------------------------------------------
! Sort the on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
    call dsortindxa1('A',ncsf,hii,indx)

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
! Determine the configurations corresponding to the selected CSFs
!----------------------------------------------------------------------
    ! Initialisation
    isurvive=0

    ! For now we shall hardwire the number of initially selected
    ! CSFs to be 2x the no. roots
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
! Number of selected configurations
!----------------------------------------------------------------------
    nsel=sum(isurvive)
    nunsel=nconf-nsel
    
!----------------------------------------------------------------------
! Fill in the array of selected configurations
!----------------------------------------------------------------------
    allocate(isel(nsel))
    allocate(iunsel(nunsel))
    
    count_sel=0
    count_unsel=0
  
    do i=1,nconf
       if (isurvive(i) == 1) then
          count_sel=count_sel+1
          isel(count_sel)=i
       else
          count_unsel=count_unsel+1
          iunsel(count_unsel)=i
       endif
    enddo

!----------------------------------------------------------------------
! Fill in the CSF offsets for the selected configurations
!----------------------------------------------------------------------
    allocate(offset_sel(nsel+1))
    allocate(offset_unsel(nunsel+1))
    offset_sel=0
    offset_unsel=0
    
    total_sel=1
    total_unsel=1
    count_sel=0
    count_unsel=0
        
    do iconf=1,nconf
       if (isurvive(iconf) == 1) then
          count_sel=count_sel+1
          offset_sel(count_sel)=total_sel
          total_sel=total_sel+offset(iconf+1)-offset(iconf)
       else
          count_unsel=count_unsel+1
          offset_unsel(count_unsel)=total_unsel
          total_unsel=total_unsel+offset(iconf+1)-offset(iconf)
       endif
    enddo

    offset_sel(nsel+1)=total_sel
    offset_unsel(nunsel+1)=total_unsel
    
    return
  
  end subroutine init_subspace

!######################################################################

  subroutine hii_selected(csfdim,hii,csfdim_sel,csfdim_unsel,hii_sel,&
       hii_unsel,nconf,offset,nsel,isel,nunsel,iunsel,averageii,&
       averageii_sel)

    use constants
    use bitglobal

    implicit none
    
    ! On-diagonal Hamiltonian matrix elements for the full space
    integer(is), intent(in)  :: csfdim
    real(dp), intent(in)     :: hii(csfdim)
 
    ! On-diagonal Hamiltonian matrix elements for the selected
    ! and unselected subspaces
    integer(is), intent(out) :: csfdim_sel,csfdim_unsel
    real(dp), allocatable    :: hii_sel(:),hii_unsel(:)

    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: averageii(nconf)

    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    ! for the selected subspace
    real(dp), allocatable    :: averageii_sel(:)
    
    ! CSF offsets
    integer(is), intent(in)  :: nconf
    integer(is), intent(in)  :: offset(nconf+1)

    ! Selected and unselected configurations
    integer(is), intent(in)  :: nsel,nunsel
    integer(is), intent(in)  :: isel(nsel),iunsel(nunsel)
    
    ! Everything else
    integer(is)              :: i,iconf,icsf,ncsf,count

!----------------------------------------------------------------------
! Determine the dimension of the Hamiltonian projected onto the
! selected subspace
!----------------------------------------------------------------------
    ! Initialisation
    csfdim_sel=0

    ! Loop over the selected configurations
    do i=1,nsel
       iconf=isel(i)

       ! Number of CSFs generated by this configuration
       ncsf=offset(iconf+1)-offset(iconf)

       ! Running total number of selected CSFs
       csfdim_sel=csfdim_sel+ncsf
       
    enddo

    ! Number of unselected CSFs
    csfdim_unsel=csfdim-csfdim_sel

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(hii_sel(csfdim_sel))
    allocate(hii_unsel(csfdim_unsel))
    allocate(averageii_sel(nsel))
    hii_sel=0.0d0
    hii_unsel=0.0d0
    averageii_sel=0.0d0
    
!----------------------------------------------------------------------    
! Fill in the selected on-diagonal matrix element array
!----------------------------------------------------------------------
    ! Selected CSF counter
    count=0
    
    ! Loop over selected configurations
    do i=1,nsel
       iconf=isel(i)

       ! Loop over the CSFs generated by this configuration
       do icsf=offset(iconf),offset(iconf+1)-1

          ! Update the selected CSF counter
          count=count+1

          ! Save the on-diagonal Hamiltonian matrix element
          ! for this CSF
          hii_sel(count)=hii(icsf)
          
       enddo

       ! Save the spin-coupling averaged on-diagonal Hamiltonian
       ! matrix element value for this configuration
       averageii_sel(i)=averageii(iconf)
              
    enddo

!----------------------------------------------------------------------    
! Fill in the unselected on-diagonal matrix element array
!----------------------------------------------------------------------
    ! Unselected CSF counter
    count=0

    ! Loop over unselected configurations
    do i=1,nunsel
       iconf=iunsel(i)

       ! Loop over the CSFs generated by this configuration
       do icsf=offset(iconf),offset(iconf+1)-1

          ! Update the unselected CSF counter
          count=count+1

          ! Save the on-diagonal Hamiltonian matrix element
          ! for this CSF
          hii_unsel(count)=hii(icsf)
          
       enddo
       
    enddo
    
    return
    
  end subroutine hii_selected
  
!######################################################################
! pt2_corrections: Calculation of the ENPT2 corrections
!######################################################################
  subroutine pt2_corrections(csfdim_sel,csfdim_unsel,nroots,vec_sel,&
       ener_sel,hii_unsel,Avec,E2,n_int_I,nconf,conf,sop,nsel,nunsel,&
       isel,iunsel,offset,averageii,harr2dim,harr2,m2c,offset_sel,&
       offset_unsel)

    use constants
    use bitglobal
    use conftype
    use mrci_integrals
    use mrciutils
    use hbuild_mrci
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: csfdim_sel,csfdim_unsel,nroots

    ! Eigenpairs of the reference space Hamiltonian projected
    ! onto the selected subspace
    real(dp), intent(in)    :: vec_sel(csfdim_sel,nroots)
    real(dp), intent(in)    :: ener_sel(nroots)

    ! On-diagonal Hamiltonian matrix elements for the unselected
    ! CSFs
    real(dp), intent(in)    :: hii_unsel(csfdim_unsel)
    
    ! A-vectors
    real(dp), intent(out)   :: Avec(csfdim_unsel,nroots)

    ! 2nd-order energy corrections
    real(dp), intent(out)   :: E2(nroots)
    
    ! Configurations and SOPs
    integer(is), intent(in) :: nconf,n_int_I
    integer(ib), intent(in) :: conf(n_int_I,2,nconf)
    integer(ib), intent(in) :: sop(n_int_I,2,nconf)

    ! Selected/unselected confs
    integer(is), intent(in) :: nsel,nunsel
    integer(is), intent(in) :: isel(nsel),iunsel(nunsel)

    ! CSF offsets
    integer(is), intent(in) :: offset(nconf+1)

    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: averageii(nconf)

    ! Work arrays
    integer(is), intent(in) :: harr2dim
    real(dp)                :: harr2(harr2dim)

    ! MO mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! CSF offsets for the selected and unselected configurations
    integer(is), intent(in) :: offset_sel(nsel+1)
    integer(is), intent(in) :: offset_unsel(nunsel+1)
    
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
    real(dp)                :: ediff
    
!----------------------------------------------------------------------
! Compute the matrix elements <w omega|H|Psi_I^0>
!----------------------------------------------------------------------
    ! Initialisation
    Avec=0.0d0

    ! Loop over the ket (selected) configurations
    do iket=1,nsel
       kconf=isel(iket)
       
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

       ! Loop over the bra (unselected) configurations
       do ibra=1,nunsel
          bconf=iunsel(ibra)

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

             ! Loop over the ket (selected) CSFs
             do ikcsf=offset_sel(iket),offset_sel(iket+1)-1

                ! Cycle if the ket CSF coefficient is tiny
                if (abs(vec_sel(ikcsf,iroot)) < epshij) cycle
                
                ! Loop over the ket (selected) CSFs
                do ibcsf=offset_unsel(ibra),offset_unsel(ibra+1)-1
                   counter=counter+1

                   Avec(ibcsf,iroot)=Avec(ibcsf,iroot)&
                        +harr2(counter)*vec_sel(ikcsf,iroot)
                   
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

       ! Loop over the unselected CSFs
       do icsf=1,csfdim_unsel

          ! E^(0) - H_ii
          ediff=ener_sel(iroot)-hii_unsel(icsf)

          ! Energy correction
          E2(iroot)=E2(iroot)+Avec(icsf,iroot)**2/ediff

          ! A-vector element
          Avec(icsf,iroot)=Avec(icsf,iroot)/ediff
          
       enddo

    enddo
    
    return
    
  end subroutine pt2_corrections

!######################################################################
  
end module ref_sci
