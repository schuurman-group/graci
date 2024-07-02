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
    use gendav
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
    integer(is)                :: hdim
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
    ! subspace
    integer(is)                :: hdim_sel
    real(dp), allocatable      :: hii_sel(:)
    
    ! Selected subspace information
    integer(is)                :: nsel_conf
    integer(is), allocatable   :: isel_conf(:)
    
    ! Everything else
    integer(is)                :: i
  
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    if (verbose) then
       write(6,'(/,72a)') ('-',i=1,52)
       write(6,'(3(x,a))') &
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

    call basis_dimensions(hdim,offset,sop,n_int_I,nconf1)

!----------------------------------------------------------------------
! Sanity check on the number of roots
!----------------------------------------------------------------------
    if (hdim < nroots) then
       errmsg='Error in ref_diag_mrci: N_roots > N_CSF'
       call error_control
    endif

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(hii(hdim))
    allocate(averageii(nconf1))
    hii=0.0d0
    averageii=0.0d0

!----------------------------------------------------------------------
! Compute the on-diagonal Hamiltonian matrix elements and the
! their spin-coupling averaged values
!----------------------------------------------------------------------
    call hii_double(nconf1,hdim,offset,conf,sop,n_int_I,m2c,irrep,hii,&
         averageii)

!----------------------------------------------------------------------
! Determine the initial subspace
!----------------------------------------------------------------------
    call init_subspace(nroots,nconf1,offset,hdim,hii,nsel_conf,&
         isel_conf)

!----------------------------------------------------------------------
! Diagonalisation of the reference space Hamiltonian projected onto
! the initial subspace
!----------------------------------------------------------------------
    ! Fill in the selected on-diagonal Hamiltonian matrix elements
    call hii_selected(hdim,hii,hdim_sel,hii_sel,nconf1,offset,&
         nsel_conf,isel_conf)
    
    ! Save to disk the non-zero selected off-diagonal Hamiltonian
    ! matrix elements
    call save_hij_double_selected(nsel_conf,isel_conf,nconf1,hdim,&
         offset,averageii,conf,sop,n_int_I,m2c,irrep,hscr,nrec,&
         'hij_ref')

    ! Diagonalise the projected reference space Hamiltonian
    call diag_full(hscr,nrec,hdim,hii,irrep,nroots,vecscr,'refvec')

!----------------------------------------------------------------------
! Iterative improvement of the subspace
!----------------------------------------------------------------------
    
    

    STOP

    
    return
  
  end subroutine ref_diag_mrci_sci

!######################################################################

  subroutine init_subspace(nroots,nconf,offset,ncsf,hii,nsel_conf,&
       isel_conf)

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
    integer(is), intent(out) :: nsel_conf
    integer(is), allocatable :: isel_conf(:)

    ! Sorting arrays
    integer(is), allocatable :: indx(:)
    
    ! Everything else
    integer(is)              :: iconf,icsf,i,n
    integer(is)              :: nsel_csf
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
    nsel_csf=2*nroots

    ! Loop over the lowest energy CSFs
    do i=1,nsel_csf

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
    nsel_conf=sum(isurvive)

!----------------------------------------------------------------------
! Fill in the array of selected configurations
!----------------------------------------------------------------------
    allocate(isel_conf(nsel_conf))

    n=0
  
    do i=1,nconf
       if (isurvive(i) == 1) then
          n=n+1
          isel_conf(n)=i
       endif
    enddo
    
    return
  
  end subroutine init_subspace

!######################################################################

  subroutine hii_selected(hdim,hii,hdim_sel,hii_sel,nconf,offset,&
       nsel_conf,isel_conf)

    use constants
    use bitglobal

    implicit none
    
    ! On-diagonal Hamiltonian matrix elements for the full space
    integer(is), intent(in)  :: hdim
    real(dp), intent(in)     :: hii(hdim)

    ! On-diagonal Hamiltonian matrix elements for the selected
    ! subspace
    integer(is), intent(out) :: hdim_sel
    real(dp), allocatable    :: hii_sel(:)

    ! CSF offsets
    integer(is), intent(in)  :: nconf
    integer(is), intent(in)  :: offset(nconf+1)

    ! Selected configurations
    integer(is), intent(in)  :: nsel_conf
    integer(is), intent(in)  :: isel_conf(nsel_conf)
    
    ! Everything else
    integer(is)              :: i,iconf,icsf,ncsf,count

!----------------------------------------------------------------------
! Determine the dimension of the Hamiltonian projected onto the
! selected subspace
!----------------------------------------------------------------------
    ! Initialisation
    hdim_sel=0

    ! Loop over the selected configurations
    do i=1,nsel_conf
       iconf=isel_conf(i)

       ! Number of CSFs generated by this configuration
       ncsf=offset(iconf+1)-offset(iconf)

       ! Running total number of selected CSFs
       hdim_sel=hdim_sel+ncsf
       
    enddo

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(hii_sel(hdim_sel))
    hii_sel=0.0d0
    
!----------------------------------------------------------------------    
! Fill in the selected on-diagonal matrix element array
!----------------------------------------------------------------------
    ! Selected CSF counter
    count=0
    
    ! Loop over selected configurations
    do i=1,nsel_conf
       iconf=isel_conf(i)

       ! Loop over the CSFs generated by this configuration
       do icsf=offset(iconf),offset(iconf+1)-1

          ! Update the selected CSF counter
          count=count+1

          ! Save the on-diagonal Hamiltonian matrix element
          ! for this CSF
          hii_sel(count)=hii(icsf)
          
       enddo
       
    enddo

    return
    
  end subroutine hii_selected
  
!######################################################################
  
end module ref_sci
