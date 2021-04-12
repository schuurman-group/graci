!**********************************************************************
! Module for the calculation of MRCI guess vectors for use in an
! interative diagonalisation calculation
!**********************************************************************
module mrci_guess

  implicit none
  
contains

!######################################################################
! mrci_guess_unit: Generation of guess vectors corresponding to single
!                  CSFs
!######################################################################
  subroutine mrci_guess_unit(guessscr,nguess,dim,subdim,hdiag)

    use constants
    use bitglobal
    use utils
    use iomod
    
    implicit none

    ! Guess vector scratch file number
    integer(is), intent(out) :: guessscr

    ! Number of guess vectors
    integer(is), intent(in)  :: nguess
    
    ! Full space and subspace dimensions
    integer(is), intent(in)  :: dim,subdim

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: hdiag(dim)

    ! Indices of the diagonals of the full Hamiltonian in
    ! ascending order
    integer(is), allocatable :: indx(:)

    ! Guess vectors
    real(dp), allocatable    :: guessvec(:,:)

    ! I/O variables
    integer(is)              :: iscratch
    character(len=60)        :: guessfile
    
    ! Everything else
    integer(is)              :: i

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(indx(dim))    
    indx=0

    allocate(guessvec(dim,nguess))
    guessvec=0.0d0
    
!----------------------------------------------------------------------
! Sort the on-diagonal matrix elements in order of increasing value
!----------------------------------------------------------------------
    call dsortindxa1('A',dim,hdiag,indx)

!----------------------------------------------------------------------
! Single CSF guess vectors
!----------------------------------------------------------------------
    do i=1,nguess
       guessvec(indx(i),i)=1.0d0
    enddo

!----------------------------------------------------------------------
! Register the guess vector scratch file
!----------------------------------------------------------------------
    call scratch_name('guessvecs',guessfile)
    call register_scratch_file(guessscr,guessfile)
    
!----------------------------------------------------------------------
! Write the guess vectors to disk
!----------------------------------------------------------------------
    ! Open the scratch file
    iscratch=scrunit(guessscr)
    open(iscratch,file=scrname(guessscr),form='unformatted',&
         status='unknown')
    
    ! Dimensions
    write(iscratch) dim
    write(iscratch) nguess
    
    ! Guess vectors
    do i=1,nguess
       write(iscratch) guessvec(:,i)
    enddo
    
    ! Close the scratch file
    close(iscratch)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(indx)
    
    return
    
  end subroutine mrci_guess_unit
  
!######################################################################
! mrci_guess_subspace: Generation of guess vectors via a subspace
!                      diagonalisation
!######################################################################
  subroutine mrci_guess_subspace(guessscr,nguess,cfg,dim,subdim,hdiag,&
       averageii,confdim)
    
    use constants
    use bitglobal
    use conftype
    use utils
    use iomod
        
    implicit none

    ! Guess vector scratch file number
    integer(is), intent(out) :: guessscr

    ! Number of guess vectors
    integer(is), intent(in)  :: nguess
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg
    
    ! Full space and subspace dimensions
    integer(is), intent(in)  :: dim,subdim

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: hdiag(dim)

    ! On-diagonal Hamiltonian matrix elements
    ! averaged over spin couplings
    integer(is), intent(in)  :: confdim
    real(dp), intent(in)     :: averageii(confdim)
    
    ! Indices of the diagonals of the full Hamiltonian in
    ! ascending order
    integer(is), allocatable :: indx(:)

    ! Working arrays
    integer(is), allocatable :: isub(:)

    ! Subspace Hamiltonian matrix
    real(dp), allocatable    :: subhmat(:,:)

    ! Subspace eigenpairs
    real(dp), allocatable    :: subvec(:,:),subeig(:)
    
    ! Subspace configurations and spin-couplings
    integer(ib), allocatable :: conf(:,:,:),sop(:,:,:)
    integer(is), allocatable :: omega(:)

    ! Subspace on-diagonal Hamiltonian matrix elements
    real(dp), allocatable    :: subhii(:)
    
    ! Subspace spin-coupling-averaged on-diagonal
    ! Hamiltonian matrix element values
    real(dp), allocatable    :: subavii(:)
    
    ! Everything else
    integer(is)              :: i
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(indx(dim))    
    indx=0

    allocate(isub(dim))
    isub=0

    allocate(conf(n_int,2,subdim))
    conf=0_ib

    allocate(sop(n_int,2,subdim))
    sop=0_ib

    allocate(omega(subdim))
    omega=0

    allocate(subhii(subdim))
    subhii=0.0d0
    
    allocate(subavii(subdim))
    subavii=0.0d0
    
    allocate(subhmat(subdim,subdim))
    subhmat=0.0d0

    allocate(subvec(subdim,subdim))
    subvec=0.0d0

    allocate(subeig(subdim))
    subeig=0.0d0
    
!----------------------------------------------------------------------
! Sort the on-diagonal matrix elements in order of increasing value
!----------------------------------------------------------------------
    call dsortindxa1('A',dim,hdiag,indx)

!----------------------------------------------------------------------
! Flag the subspace CSFs
!----------------------------------------------------------------------
    do i=1,subdim
       isub(indx(i))=1
    enddo
    
!----------------------------------------------------------------------
! Get the subspace CSF information (configuration bit strings and
! spin coupling indices)
!----------------------------------------------------------------------
    call subspace_info(cfg,dim,subdim,isub,conf,sop,omega,&
         hdiag,confdim,averageii,subhii,subavii)
    
!----------------------------------------------------------------------
! Calculate the elements of the subspace Hamiltonian matrix
!----------------------------------------------------------------------
    call subspace_hamiltonian(cfg,subdim,conf,sop,omega,subhii,&
         subavii,subhmat)
    
!----------------------------------------------------------------------
! Diagonalise the subspace Hamiltonian matrix
!----------------------------------------------------------------------
    call subspace_diagonalise(subdim,subhmat,subvec,subeig)

!----------------------------------------------------------------------
! Save the guess vectors to disk
!----------------------------------------------------------------------
    call write_guess_vectors(guessscr,nguess,dim,subdim,subvec,&
         subeig,isub)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(indx)
    deallocate(isub)
    deallocate(conf)
    deallocate(sop)
    deallocate(omega)
    deallocate(subhmat)
    deallocate(subvec)
    deallocate(subeig)
    deallocate(subhii)
    deallocate(subavii)
    
    return
    
  end subroutine mrci_guess_subspace
    
!######################################################################
! subspace_info: determines the spatial configurations and
!                spin couplings for the subspace CSFs
!######################################################################
  subroutine subspace_info(cfg,dim,subdim,isub,conf,sop,omega,&
       hdiag,confdim,averageii,subhii,subavii)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg
    
    ! Full space and subspace dimensions
    integer(is), intent(in)  :: dim,subdim

    ! CSF subspace flags
    integer(is), intent(in)  :: isub(dim)

    ! Subspace configurations and spin-couplings
    integer(ib), intent(out) :: conf(n_int,2,subdim)
    integer(ib), intent(out) :: sop(n_int,2,subdim)
    integer(is), intent(out) :: omega(subdim)

    ! On-diagonal Hamiltonian matrix elements
    integer(is), intent(in)  :: confdim
    real(dp), intent(in)     :: hdiag(dim)
    real(dp), intent(in)     :: averageii(confdim)
    
    ! Spin-coupling-averaged on-diagonal
    ! Hamiltonian matrix element values
    real(dp), intent(out)    :: subhii(subdim)
    real(dp), intent(out)    :: subavii(subdim)
    
    ! Everything else
    integer(is)              :: n,ioff,csf,counter,isp,k
    integer(is)              :: n_int_I
    integer(is)              :: shift
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    counter=0
    k=0

!----------------------------------------------------------------------
! Reference space configurations
!----------------------------------------------------------------------
    n_int_I=cfg%n_int_I

    ! Loop over ref space configurations
    do n=1,cfg%n0h

       ! Loop over the CSFs generated by this configuration
       isp=0
       do csf=cfg%csfs0h(n),cfg%csfs0h(n+1)-1

          ! Spin-coupling index
          isp=isp+1
          
          ! Increment the CSF counter
          k=k+1

          ! Save the subspace CSF information
          if (isub(k) == 1) then
             counter=counter+1
             omega(counter)=isp
             conf(1:n_int_I,:,counter)=cfg%conf0h(:,:,n)
             subhii(counter)=hdiag(k)
             subavii(counter)=averageii(n)
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
    shift=cfg%n0h

    if (cfg%n1I > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1I configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1I(n),cfg%off1I(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             isp=0
             do csf=cfg%csfs1I(ioff),cfg%csfs1I(ioff+1)-1

                ! Spin-coupling index
                isp=isp+1
                
                ! Increment the CSF counter
                k=k+1
                
                ! Save the subspace CSF information
                if (isub(k) == 1) then
                   counter=counter+1
                   omega(counter)=isp
                   conf(:,:,counter)=cfg%conf1I(:,:,ioff)
                   subhii(counter)=hdiag(k)
                   subavii(counter)=averageii(ioff+shift)
                endif
                
             enddo

          enddo
          
       enddo

    endif

!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
    shift=cfg%n0h+cfg%n1I

    if (cfg%n2I > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2I configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2I(n),cfg%off2I(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             isp=0
             do csf=cfg%csfs2I(ioff),cfg%csfs2I(ioff+1)-1

                ! Spin-coupling index
                isp=isp+1
                
                ! Increment the CSF counter
                k=k+1
                
                ! Save the subspace CSF information
                if (isub(k) == 1) then
                   counter=counter+1
                   omega(counter)=isp
                   conf(:,:,counter)=cfg%conf2I(:,:,ioff)
                   subhii(counter)=hdiag(k)
                   subavii(counter)=averageii(ioff+shift)
                endif
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
    shift=cfg%n0h+cfg%n1I+cfg%n2I

    if (cfg%n1E > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1E configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1E(n),cfg%off1E(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             isp=0
             do csf=cfg%csfs1E(ioff),cfg%csfs1E(ioff+1)-1

                ! Spin-coupling index
                isp=isp+1
                
                ! Increment the CSF counter
                k=k+1
                
                ! Save the subspace CSF information
                if (isub(k) == 1) then
                   counter=counter+1
                   omega(counter)=isp
                   conf(:,:,counter)=cfg%conf1E(:,:,ioff)
                   subhii(counter)=hdiag(k)
                   subavii(counter)=averageii(ioff+shift)
                endif
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
    shift=cfg%n0h+cfg%n1I+cfg%n2I+cfg%n1E

    if (cfg%n2E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2E(n),cfg%off2E(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             isp=0
             do csf=cfg%csfs2E(ioff),cfg%csfs2E(ioff+1)-1

                ! Spin-coupling index
                isp=isp+1
                
                ! Increment the CSF counter
                k=k+1
                
                ! Save the subspace CSF information
                if (isub(k) == 1) then
                   counter=counter+1
                   omega(counter)=isp
                   conf(:,:,counter)=cfg%conf2E(:,:,ioff)
                   subhii(counter)=hdiag(k)
                   subavii(counter)=averageii(ioff+shift)
                endif
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
    shift=cfg%n0h+cfg%n1I+cfg%n2I+cfg%n1E+cfg%n2E

    if (cfg%n1I1E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 1I1E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             isp=0
             do csf=cfg%csfs1I1E(ioff),cfg%csfs1I1E(ioff+1)-1

                ! Spin-coupling index
                isp=isp+1
                
                ! Increment the CSF counter
                k=k+1
                
                ! Save the subspace CSF information
                if (isub(k) == 1) then
                   counter=counter+1
                   omega(counter)=isp
                   conf(:,:,counter)=cfg%conf1I1E(:,:,ioff)
                   subhii(counter)=hdiag(k)
                   subavii(counter)=averageii(ioff+shift)
                endif
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! Subspace SOPs
!----------------------------------------------------------------------
    ! Loop over subspace configurations
    do n=1,subdim

       ! Compute the SOP
       sop(:,:,n)=conf_to_sop(conf(:,:,n),n_int)
       
    enddo

    return
    
  end subroutine subspace_info

!######################################################################
! subspace_hamiltonian: constructs the Hamiltonian projected onto
!                       the subspace
!######################################################################
  subroutine subspace_hamiltonian(cfg,subdim,conf,sop,omega,subhii,&
       subavii,subhmat)
    
    use constants
    use bitglobal
    use conftype
    use hparam
    use hbuild_mrci
    use mrci_integrals
    use dftmrci
    use mrciutils
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg
    
    ! Subspace dimension
    integer(is), intent(in)  :: subdim

    ! Subspace configurations and spin-couplings
    integer(ib), intent(out) :: conf(n_int,2,subdim)
    integer(ib), intent(out) :: sop(n_int,2,subdim)
    integer(is), intent(out) :: omega(subdim)

    ! Diagonal matrix elements of the subspace Hamiltonian
    real(dp), intent(in)     :: subhii(subdim)

    ! Subspace spin-coupling-averaged on-diagonal
    ! Hamiltonian matrix element values
    real(dp), intent(in)     :: subavii(subdim)
    
    ! Subspace Hamiltonian matrix
    real(dp), intent(out)    :: subhmat(subdim,subdim)

    ! Hamiltonian build variables
    integer(is)              :: nexci
    integer(is)              :: socc(nmo),docc(nmo),unocc(nmo)
    integer(is)              :: nopen,nsocc,ndocc,nunocc
    integer(is)              :: knopen,bnopen
    integer(is), parameter   :: maxexci=2
    integer(is)              :: hlist(maxexci),plist(maxexci)

    ! Difference configuration information
    integer(is)              :: ndiff
    integer(is)              :: Dw(nmo,2)

    ! Number of open shells preceding each MO
    integer(is)              :: nbefore(nmo)

    ! Spin-coupling sub-case bitsting encodings
    integer(ib)              :: pairindx(nmo)
    integer(ib)              :: icase
    
    ! Pattern indices
    integer(is)              :: bpattern(nmo+1),kpattern(nmo+1)

    ! Integrals
    real(dp)                 :: Vpqrs(nmo)
    
    ! Everything else
    integer(is)              :: i,bcsf,kcsf
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    subhmat=0.0d0

!----------------------------------------------------------------------
! Diagonal elements
!----------------------------------------------------------------------
    do i=1,subdim
       subhmat(i,i)=subhii(i)
    enddo

!----------------------------------------------------------------------
! Off-diagonal elements
!----------------------------------------------------------------------
    ! Loop over ket CSFs
    do kcsf=1,subdim-1
    
       ! Number of open shells in the ket CSF
       knopen=sop_nopen(sop(:,:,kcsf),n_int)

       ! Package the ket configuration information
       call package_confinfo(sop(:,:,kcsf),conf(:,:,kcsf),&
            unocc,socc,docc,nunocc,nsocc,ndocc,Dw,ndiff,nbefore)
       
       ! Loop over bra CSFs
       do bcsf=kcsf+1,subdim
       
          ! Compute the excitation degree between the two
          ! configurations
          nexci=exc_degree_conf(conf(:,:,kcsf),conf(:,:,bcsf),n_int)

          ! Cycle if the excitation degree is greater than 2
          if (nexci > 2) cycle

          ! Number of open shells in the bra configuration
          bnopen=sop_nopen(sop(:,:,bcsf),n_int)

          ! Get the indices of the MOs involved in the excitation
          hlist=0
          plist=0
          call get_exci_indices(conf(:,:,kcsf),conf(:,:,bcsf),&
               n_int,hlist(1:nexci),plist(1:nexci),nexci)

          !************************************************************
          ! Note that here the lists of holes/particles are the indices
          ! of the annihilation/creation operators operating on the ket
          ! configuration to yield the bra configuration
          !************************************************************

          ! Fill the integrals, pattern index, and pair index arrays
          select case(nexci)
          case(1)
             call package_integrals_nexci1(&
                  sop(:,:,bcsf),sop(:,:,kcsf),&
                  pairindx(1),hlist(1),plist(1),bnopen,knopen,&
                  bpattern,kpattern,Vpqrs,cfg%m2c,socc,nsocc,nbefore,&
                  Dw,ndiff,icase)
          case(2)
             call package_integrals_nexci2(&
                  sop(:,:,bcsf),sop(:,:,kcsf),&
                  pairindx(1:2),hlist(1:2),plist(1:2),bnopen,knopen,&
                  bpattern(1:2),kpattern(1:2),Vpqrs(1:2),cfg%m2c,&
                  nbefore)
          end select

          ! Compute the matrix element
          select case(nexci)
          case(0) ! Same spatial configuration, different
                  ! spin couplings
             subhmat(bcsf,kcsf)=hij_same_mrci_1element(&
                  omega(bcsf),omega(kcsf),sop(:,:,kcsf),socc,&
                  nsocc,nbefore,cfg%m2c)
          case(1) ! Bra and ket configurations linked by a single
                  ! excitation
             subhmat(bcsf,kcsf)=hij_single_mrci(&
                  omega(bcsf),omega(kcsf),bnopen,knopen,pairindx(1),&
                  icase,bpattern,kpattern,Vpqrs,socc,nsocc,&
                  ndiff,hlist(1),plist(1))
          case(2) ! Bra and ket configurations linked by two
                  ! excitations
             subhmat(bcsf,kcsf)=hij_double_mrci(&
                  omega(bcsf),omega(kcsf),bnopen,knopen,pairindx(1:2),&
                  bpattern(1:2),kpattern(1:2),Vpqrs(1:2),&
                  plist(1:2),hlist(1:2))
          end select
          
          ! DFT/MRCI corrections
          ! Note that we will have to take account of the nexci=0
          ! case once we implement the Dusseldorf group Hamiltonians...
          if (nexci > 0) then
             call hij_dftmrci_batch(subhmat(bcsf:bcsf,kcsf:kcsf),1,1,&
                  subavii(bcsf),subavii(kcsf))
          endif
          
          ! Fill in the upper triangle
          subhmat(kcsf,bcsf)=subhmat(bcsf,kcsf)

       enddo
       
    enddo
    
    return
    
  end subroutine subspace_hamiltonian

!######################################################################
! subspace_diagonalise: full diagonalisation of the subspace
!                       Hamiltonian matrix
!######################################################################
  subroutine subspace_diagonalise(subdim,subhmat,subvec,subeig)

    use constants
    use bitglobal
    use utils
    
    implicit none

    ! No. subspace CSFs
    integer(is), intent(in) :: subdim

    ! Subspace Hamiltonian matrix
    real(dp), intent(in)    :: subhmat(subdim,subdim)

    ! Subspace eigenpairs
    real(dp), intent(out)   :: subvec(subdim,subdim)
    real(dp), intent(out)   :: subeig(subdim)

    ! Everything else
    integer(is)             :: i
        
!----------------------------------------------------------------------
! Diagonalise the subspace Hamiltonian matrix
!----------------------------------------------------------------------
    call diag_matrix_real(subhmat,subeig,subvec,subdim)

    ! Add on E_SCF + E_nuc
    subeig=subeig+escf+enuc
    
    return
    
  end subroutine subspace_diagonalise
    
!######################################################################
! write_guess_vectors: writes the guess vectors to disk
!######################################################################
  subroutine write_guess_vectors(guessscr,nguess,dim,subdim,subvec,&
       subeig,isub)

    use constants
    use bitglobal
    use iomod
    
    implicit none

    ! Guess vector scratch file number
    integer(is), intent(out) :: guessscr

    ! Number of guess vectors
    integer(is), intent(in)  :: nguess
    
    ! Full and subspace dimensions
    integer(is), intent(in)  :: dim,subdim

    ! Subspace eigenpairs
    real(dp), intent(in)     :: subvec(subdim,subdim)
    real(dp), intent(in)     :: subeig(subdim)
    
    ! Indices of the subspace CSFs
    integer(is), intent(in)  :: isub(dim)
    integer(is)              :: imap(subdim)
    
    ! Guess vectors
    real(dp), allocatable    :: guess(:,:)

    ! I/O variables
    integer(is)              :: iscratch
    character(len=60)        :: guessfile
    
    ! Everything else
    integer(is)              :: i,j,j1,counter

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(guess(dim,nguess))
    guess=0.0d0

!----------------------------------------------------------------------
! Get the subspace-to-full space mapping
!----------------------------------------------------------------------
    imap=0
    counter=0
    do i=1,dim
       if (isub(i) == 1) then
          counter=counter+1
          imap(counter)=i
       endif
    enddo
    
!----------------------------------------------------------------------
! Construct the guess vectors
!----------------------------------------------------------------------
    ! Loop over guess vectors
    do i=1,nguess

       ! Loop over subspace CSFs
       do j=1,subdim

          ! Index in the full space of the current subspace CSF
          j1=imap(j)

          ! Fill in the guess vector
          guess(j1,i)=subvec(j,i)
          
       enddo
       
    enddo
    
!----------------------------------------------------------------------
! Register the guess vector scratch file
!----------------------------------------------------------------------
    call scratch_name('guessvecs',guessfile)
    call register_scratch_file(guessscr,guessfile)

!----------------------------------------------------------------------
! Write the guess vectors to disk
!----------------------------------------------------------------------
    ! Open the scratch file
    iscratch=scrunit(guessscr)
    open(iscratch,file=scrname(guessscr),form='unformatted',&
         status='unknown')

    ! Dimensions
    write(iscratch) dim
    write(iscratch) nguess

    ! Guess vectors
    do i=1,nguess
       write(iscratch) guess(:,i)
    enddo

    ! Subspace eigenpairs (needed if the generalised Davidson
    ! eigensolver is being used)
    write(iscratch) subdim
    write(iscratch) imap
    write(iscratch) subeig
    write(iscratch) subvec
    
    ! Close the scratch file
    close(iscratch)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(guess)

    return
    
  end subroutine write_guess_vectors
  
!######################################################################
  
end module mrci_guess
