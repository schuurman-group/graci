!***********************************************************************
! Routines for the application of the DFT/MRCI corrections to the
! Hamiltonian matrix
!***********************************************************************
module dftmrci

  implicit none

contains
  
!######################################################################
! hii_dftmrci: applies DFT/MRCI corrections to a batch of on-diagonal
!              Hamiltonian matrix elements
!######################################################################
  subroutine hii_dftmrci(harr,nsp,Dw,ndiff,nopen,m2c,&
       sop,socc,nsocc,nbefore)
    
    use constants
    use bitglobal
    use hparam
    use iomod
    
    implicit none

    ! Array of on-diagonal Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: harr(nsp)

    ! Difference configuration information
    integer(is), intent(in) :: ndiff
    integer(is), intent(in) :: Dw(nmo,2)

    ! Number of open shells
    integer(is), intent(in) :: nopen
    
    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! SOP characterising the spatial occupation
    integer(ib), intent(in) :: sop(n_int,2)

    ! Indices of the singly-occupied MOs
    integer(is), intent(in) :: nsocc
    integer(is), intent(in) :: socc(nmo)

    ! Numbers of open shells preceding each MO
    integer(is), intent(in) :: nbefore(nmo)
    
    select case(ihamiltonian)
       
    case(2:3)
       ! Grimme's parameterisation
       call hii_dftmrci_grimme(harr,nsp,Dw,ndiff,nopen,m2c)

    case(4:5)
       ! Lyskov's parameterisation
       call hii_dftmrci_lyskov(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
            nsocc,nbefore)

    case(6:9)
       ! Heil's parameterisations
       call hii_dftmrci_heil(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
            nsocc,nbefore)

    case(10:11)
       ! Ottawa CVS parameterisation
       call hii_dftmrci_cvs(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
            nsocc,nbefore)
       
    case default
       errmsg='Your Hamiltonian choice has not been implemented yet'
       call error_control
       
    end select
    
    return
    
  end subroutine hii_dftmrci

!######################################################################
! hij_dftmrci_batch: applies DFT/MRCI corrections to a batch of
!                    off-diagonal Hamiltonian matrix elements
!######################################################################
  subroutine hij_dftmrci_batch(hij,bdim,kdim,bav,kav)

    use constants
    use bitglobal
    use hparam
    use iomod
    
    implicit none

    ! Hamiltonian matrix elements
    integer(is), intent(in) :: bdim,kdim
    real(dp), intent(inout) :: hij(:)
    real(dp), intent(in)    :: bav,kav
    real(dp)                :: damp
    
!----------------------------------------------------------------------
! Compute the damping factor
!----------------------------------------------------------------------
    select case(ihamiltonian)
       
    case(2:3)
       ! Grimme's parameterisation
       damp=damping_grimme(bav,kav)

    case(4:7, 10:11)
       ! Lyskov's parameterisation
       ! Note that this is also used for Heil's 2017 Hamiltonian
       ! as well as the Ottawa CVS Hamiltonian
       damp=damping_lyskov(bav,kav)

    case(8:9)
       ! Heil's 2018 parameterisation
       damp=damping_heil18(bav,kav)
       
    case default
       errmsg='Your Hamiltonian choice has not been implemented yet'
       call error_control
       
    end select

!----------------------------------------------------------------------
! Apply the damping factor
!----------------------------------------------------------------------
    hij(1:bdim*kdim)=hij(1:bdim*kdim)*damp
    
    return
    
  end subroutine hij_dftmrci_batch

!######################################################################
! hij_same_dftmrci: applies DFT/MRCI corrections to a batch of
!                   off-diagonal Hamiltonian matrix elements with the
!                   same spatial part but different spin couplings
!######################################################################
  subroutine hij_same_dftmrci(hij,nsp,sop,socc,nsocc,nbefore,m2c)

    use constants
    use bitglobal
    use hparam
    use iomod
    
    implicit none

    ! Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: hij(:)

    ! SOP
    integer(ib), intent(in) :: sop(n_int,2)

    ! Singly-occupied MOs
    integer(is), intent(in) :: socc(nmo)
    integer(is), intent(in) :: nsocc

    ! Numbers of open shells preceding each MO
    integer(is), intent(in) :: nbefore(nmo)

    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)
    
    ! Everything else
    integer(is)             :: nij
        
    select case(ihamiltonian)
       
    case(2:3)
       ! Grimme's parameterisation: do nothing
       return

    case(4:9)
       ! Lyskov's parameterisation
       ! Note that this is also used for Heil's Hamiltonians
       nij=nsp*(nsp-1)/2
       hij(1:nij)=(1.0d0-hpar(2))*hij(1:nij)
       return

    case(10:11)
       ! Ottawa CVS parameterisation
       call hij_same_dftmrci_cvs(hij,nsp,sop,socc,nsocc,nbefore,m2c)
       
    case default
       errmsg='Your Hamiltonian choice has not been implemented yet'
       call error_control
       
    end select
    
    return
    
  end subroutine hij_same_dftmrci

!######################################################################
! hij_same_dftmrci_cvs: applies the Ottawa CVS-DFT/MRCI corrections to
!                       a batch of off-diagonal Hamiltonian matrix
!                       elements with the same spatial part but
!                       different spin couplings
!######################################################################
  subroutine hij_same_dftmrci_cvs(hij,nsp,sop,socc,nsocc,nbefore,m2c)

    use constants
    use bitglobal
    use pattern_indices
    use hparam

    implicit none

    ! Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: hij(:)

    ! SOP
    integer(ib), intent(in) :: sop(n_int,2)

    ! Singly-occupied MOs
    integer(is), intent(in) :: socc(nmo)
    integer(is), intent(in) :: nsocc

    ! Numbers of open shells preceding each MO
    integer(is), intent(in) :: nbefore(nmo)

    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! Temporary K-edge MO energy threshold
    real(dp), parameter     :: ethrsh=-10.0d0
    
    ! Everything else
    integer(is)             :: insp,nopen
    integer(is)             :: i,i1,j,j1,ic,ja
    integer(is)             :: bomega,komega
    integer(is)             :: pattern,kstart,bstart
    integer(is)             :: count,n
    real(dp)                :: Vijji,product
    real(dp)                :: pF,pFvv,pFcv

!----------------------------------------------------------------------
! Parameter values
!----------------------------------------------------------------------
    ! Valence-valence interactions
    pFvv=hpar(2)

    ! Core-valence interactions
    pFcv=hpar(6)
    
!----------------------------------------------------------------------
! Numbers of 'intermediate' CSFs entering into the contractions of the
! fibers of the spin-coupling coefficient tensor
!----------------------------------------------------------------------
    nopen=nsocc

    if (nopen > 1) then
       insp=ncsfs(nopen-2)
    else
       insp=0
    endif

!----------------------------------------------------------------------
! Case 2b spin-coupling coefficients
!----------------------------------------------------------------------
! Note that all other <w omega' | E_i^j E_j^i | w omega> terms are zero
!----------------------------------------------------------------------
    ! Loop over singly-occupied MOs (creation operator)
    do i=1,nsocc-1

       ! Creation operator index
       ic=socc(i)
       
       ! DFT/HF MO index
       i1=m2c(ic)

       ! Loop over singly-occupied MOs (annihilation operator)
       do j=i+1,nsocc

          ! Annihilation operator index
          ja=socc(j)
          
          ! DFT/HF MO index
          j1=m2c(ja)

          ! Get the spin coupling coefficient pattern indices
          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)
          
          ! V_ijji
          Vijji=Vx(i1,j1)

          ! Are we at a valence-valence or core-valence interaction?
          if (moen(i1) < ethrsh .or. moen(j1) < ethrsh) then
             ! core-valence
             pF=pFcv
          else
             ! valence-valence
             pF=pFvv
          endif
          
          ! Contributions to hij
          count=0
          kstart=pattern
          do komega=1,nsp
             bstart=pattern
             do bomega=1,nsp
                if (bomega > komega) then
                   count=count+1
                   product=dot_product(&
                        spincp(bstart:bstart+insp-1),&
                        spincp(kstart:kstart+insp-1))
                   hij(count)=hij(count)-pF*Vijji*product
                endif
                bstart=bstart+insp
             enddo
             kstart=kstart+insp
          enddo
          
       enddo

    enddo
    
    return

  end subroutine hij_same_dftmrci_cvs
  
!######################################################################
! hii_dftmrci_grimme: applies Grimme's DFT/MRCI correction to a batch
!                     of on-diagonal Hamiltonian matrix elements
!######################################################################
  subroutine hii_dftmrci_grimme(harr,nsp,Dw,ndiff,nopen,m2c)
    
    use constants
    use bitglobal
    use hparam

    implicit none

    ! Array of on-diagonal Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: harr(nsp)

    ! Difference configuration information
    integer(is), intent(in) :: ndiff
    integer(is), intent(in) :: Dw(nmo,2)

    ! Number of open shells
    integer(is), intent(in) :: nopen
    
    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)
    
    ! Everything else
    integer(is)             :: i,j,i1,j1,Dwi,Dwj,ipos
    integer(is)             :: nexci
    real(dp)                :: Viijj,Vijji
    real(dp)                :: pJ,pN0
    real(dp)                :: contrib
    
!----------------------------------------------------------------------
! Return if we are at the base configuration
!----------------------------------------------------------------------
    if (ndiff == 0) return

!----------------------------------------------------------------------
! Parameter values
!----------------------------------------------------------------------
    pJ=hpar(3)
    pN0=hpar(4)+nopen*hpar(5)
    
!----------------------------------------------------------------------
! Sum_i F_ii^KS - F_ii^HF Delta w_i
!----------------------------------------------------------------------
    ! Loop over non-zero Delta w_i values
    do i=1,ndiff
    
       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
    
       ! Sum the contribution
       harr=harr+(moen(i1)-Fii(i1))*Dwi
       
    enddo

!----------------------------------------------------------------------
!  1/nexc Sum_i Sum_j (pJ V_iijj - p[N0] V_ijji) |Delta w_i| Delta w_j
!         w_i<0 w_j>0
!----------------------------------------------------------------------
    !
    ! Find the start of the postive Delta w_i values
    !
    do i=1,ndiff
       if (Dw(i,2) > 0) then
          ipos=i
          exit
       endif
    enddo

    !
    ! Excitation degree
    !
    nexci=sum(Dw(ipos:ndiff,2))

    !
    ! Sum_i Sum_j (pJ V_iijj - p[N0] V_ijji) |Delta w_i| Delta w_j
    ! w_i<0 w_j<0
    !
    contrib=0.0d0
    ! Loop over negative Delta w_i values
    do i=1,ipos-1

       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       ! Loop over positive Delta w_j values
       do j=ipos,ndiff

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_j value
          Dwj=Dw(j,2)

          ! V_iijj
          Viijj=Vc(i1,j1)
          
          ! V_ijji
          Vijji=Vx(i1,j1)

          ! Sum the contibution
          contrib=contrib-(pJ*Viijj-pN0*Vijji)*Dwi*Dwj
          
       enddo
    enddo

    !
    ! Make the correction
    !
    harr=harr+contrib/nexci
    
    return
    
  end subroutine hii_dftmrci_grimme

!######################################################################
! hii_dftmrci_lyskov: applies Lyskov's DFT/MRCI correction to a batch
!                     of on-diagonal Hamiltonian matrix elements
!######################################################################
  subroutine hii_dftmrci_lyskov(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
       nsocc,nbefore)
    
    use constants
    use bitglobal
    use pattern_indices
    use hparam

    implicit none

    ! Array of on-diagonal Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: harr(nsp)

    ! Difference configuration information
    integer(is), intent(in) :: ndiff
    integer(is), intent(in) :: Dw(nmo,2)

    ! Number of open shells
    integer(is), intent(in) :: nopen
    
    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! SOP characterising the spatial occupation
    integer(ib), intent(in) :: sop(n_int,2)

    ! Indices of the singly-occupied MOs
    integer(is), intent(in) :: nsocc
    integer(is), intent(in) :: socc(nmo)

    ! Numbers of open shells preceding each MO
    integer(is), intent(in) :: nbefore(nmo)
    
    ! Everything else
    integer(is)             :: i,j,i1,j1,Dwi,Dwj,ipos,insp
    integer(is)             :: ic,ja,omega,pattern,start
    real(dp)                :: Viijj,Vijji
    real(dp)                :: contrib(nsp)
    real(dp)                :: product
    real(dp)                :: pJ,pF
    
!----------------------------------------------------------------------
! Return if we are at the base configuration
!----------------------------------------------------------------------
    if (ndiff == 0) return

!----------------------------------------------------------------------
! Parameter values
!----------------------------------------------------------------------
    pJ=hpar(1)
    pF=hpar(2)

!----------------------------------------------------------------------
! Sum_i F_ii^KS - F_ii^HF Delta w_i
!----------------------------------------------------------------------
    ! Loop over non-zero Delta w_i values
    do i=1,ndiff
    
       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
    
       ! Sum the contribution
       harr=harr+(moen(i1)-Fii(i1))*Dwi
       
    enddo

!----------------------------------------------------------------------
! Find the start of the postive Delta w_i values
!----------------------------------------------------------------------
    do i=1,ndiff
       if (Dw(i,2) > 0) then
          ipos=i
          exit
       endif
    enddo
    
!----------------------------------------------------------------------
! Coulomb corrections
!----------------------------------------------------------------------
! -pJ Sum_i Sum_j Viijj, Delta w_i < 0, Delta w_j < 0
! -pJ Sum_i Sum_j Viijj, Delta w_i > 0, Delta w_j > 0
! +pJ Sum_i Sum_j Viijj, Delta w_i < 0, Delta w_j > 0
!----------------------------------------------------------------------
    contrib=0.0d0

    ! Loop over pairs of created/annihilated MOs (relative to the base
    ! configuration)
    do i=1,ndiff

       ! MO index
       i1=m2c(Dw(i,1))
       
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       do j=i,ndiff

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_i value
          Dwj=Dw(j,2)

          if (i == j .and. abs(Dwi) == 2) then
             ! Same MO index and double excitation

             ! V_iijj
             Viijj=Vc(i1,i1)

             ! Sum the contribution
             contrib=contrib-pJ*Viijj

          else if (i /= j) then
             ! Different MO indices

             ! V_iijj
             Viijj=Vc(i1,j1)

             ! Sum the contribution
             contrib=contrib-pJ*Viijj*Dwi*Dwj             

          endif
          
       enddo

    enddo

!----------------------------------------------------------------------
! Exchange correction 1
!----------------------------------------------------------------------
! -pF/2 Sum_i Sum_j V_ijji , Delta w_i < 0, Delta w_j > 0
!----------------------------------------------------------------------
    ! Loop over negative Delta w_i values
    do i=1,ipos-1

       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       ! Loop over positive Delta w_j values
       do j=ipos,ndiff

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_j value
          Dwj=Dw(j,2)

          ! V_ijij
          Vijji=Vx(i1,j1)

          ! Sum the contribution (note that the plus sign arises
          ! because Dwi * Dwj < 0)
          contrib=contrib+0.5d0*pF*Vijji*Dwi*Dwj
                    
       enddo

    enddo

!----------------------------------------------------------------------
! Exchange correction 2
!----------------------------------------------------------------------
! -pF Sum_i Sum_j V_ijji <w omega|E_i^j E_j^i|w omega>,
! j > i, i and j singly-occupied
!----------------------------------------------------------------------
    ! Numbers of 'intermediate' CSFs entering into the contractions of
    ! the fibers of the spin-coupling coefficient tensor
    if (nopen > 1) then
       insp=ncsfs(nopen-2)
    else
       insp=0
    endif

    ! Loop over singly-occupied MOs (creation operator)
    do i=1,nsocc-1
       
       ! Creation operator index
       ic=socc(i)
       
       ! DFT/HF MO index
       i1=m2c(ic)

       ! Loop over singly-occupied MOs (annihilation operator)
       do j=i+1,nsocc
          
          ! Annihilation operator index
          ja=socc(j)
          
          ! DFT/HF MO index
          j1=m2c(ja)

          ! Get the spin coupling coefficient pattern index
          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)

          ! V_ijji
          Vijji=Vx(i1,j1)

          ! Sum the contributions
          start=pattern
          do omega=1,nsp
             product=dot_product(&
                  spincp(start:start+insp-1),&
                  spincp(start:start+insp-1))
             contrib(omega)=contrib(omega)-pF*Vijji*product
             start=start+insp
          enddo
          
       enddo

    enddo

!----------------------------------------------------------------------
! Add the Coulomb and exchange corrections
!----------------------------------------------------------------------
    harr=harr+contrib
    
    return
    
  end subroutine hii_dftmrci_lyskov

!######################################################################
! hii_dftmrci_heil: applies Heil's DFT/MRCI correction to a batch of
!                   on-diagonal Hamiltonian matrix elements
!######################################################################
  subroutine hii_dftmrci_heil(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
       nsocc,nbefore)
    
    use constants
    use bitglobal
    use pattern_indices
    use hparam

    implicit none

    ! Array of on-diagonal Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: harr(nsp)

    ! Difference configuration information
    integer(is), intent(in) :: ndiff
    integer(is), intent(in) :: Dw(nmo,2)

    ! Number of open shells
    integer(is), intent(in) :: nopen
    
    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! SOP characterising the spatial occupation
    integer(ib), intent(in) :: sop(n_int,2)

    ! Indices of the singly-occupied MOs
    integer(is), intent(in) :: nsocc
    integer(is), intent(in) :: socc(nmo)

    ! Numbers of open shells preceding each MO
    integer(is), intent(in) :: nbefore(nmo)
    
    ! Everything else
    integer(is)             :: i,j,i1,j1,Dwi,Dwj,ipos,insp
    integer(is)             :: ic,ja,omega,pattern,start
    real(dp)                :: Viijj,Vijji,Viiii
    real(dp)                :: contrib(nsp)
    real(dp)                :: product
    real(dp)                :: pJ,pF

!----------------------------------------------------------------------
! Diagonal shift: 1/4 Sum_i V_iiii, i singly occupied in the base
! configuration    
!----------------------------------------------------------------------
    do i=1,nmo
       if (iopen0(i) == 1) harr=harr+0.25d0*Vc(i,i)
    enddo
    
!----------------------------------------------------------------------
! Return if we are at the base configuration
!----------------------------------------------------------------------
    if (ndiff == 0) return

!----------------------------------------------------------------------
! Parameter values
!----------------------------------------------------------------------
    pJ=hpar(1)
    pF=hpar(2)

!----------------------------------------------------------------------
! Sum_i F_ii^KS - F_ii^HF Delta w_i
!----------------------------------------------------------------------
    ! Loop over non-zero Delta w_i values
    do i=1,ndiff
    
       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
    
       ! Sum the contribution
       harr=harr+(moen(i1)-Fii(i1))*Dwi
       
    enddo

!----------------------------------------------------------------------
! Find the start of the postive Delta w_i values
!----------------------------------------------------------------------
    do i=1,ndiff
       if (Dw(i,2) > 0) then
          ipos=i
          exit
       endif
    enddo

!----------------------------------------------------------------------
! Coulomb correction 1
!----------------------------------------------------------------------
! -pJ Sum_i Sum_j Viijj, Delta w_i < 0, Delta w_j < 0
! -pJ Sum_i Sum_j Viijj, Delta w_i > 0, Delta w_j > 0
! +pJ Sum_i Sum_j Viijj, Delta w_i < 0, Delta w_j > 0
!----------------------------------------------------------------------
    contrib=0.0d0

    ! Loop over pairs of created/annihilated MOs (relative to the base
    ! configuration)
    do i=1,ndiff

       ! MO index
       i1=m2c(Dw(i,1))
       
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       do j=i,ndiff

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_i value
          Dwj=Dw(j,2)

          if (i == j .and. abs(Dwi) == 2) then
             ! Same MO index and double excitation

             ! V_iijj
             Viijj=Vc(i1,i1)

             ! Sum the contribution
             contrib=contrib-pJ*Viijj

          else if (i /= j) then
             ! Different MO indices

             ! V_iijj
             Viijj=Vc(i1,j1)

             ! Sum the contribution
             contrib=contrib-pJ*Viijj*Dwi*Dwj             

          endif
          
       enddo

    enddo

!----------------------------------------------------------------------
! Coulomb correction 2
!----------------------------------------------------------------------
! -pJ sum_i Viiii, |Delta w_i| = 1 and i indexes an open-shell in the
! base configuration
!----------------------------------------------------------------------
    ! Loop over created/annihilated MOs (relative to the base
    ! configuration)
    do i=1,ndiff

       ! MO index
       i1=m2c(Dw(i,1))
       
       ! Cycle if this MO is not singly-occupied in the base
       ! configuration
       if (iopen0(i1) == 0) cycle

       ! V_iiii
       Viiii=Vc(i1,i1)
       
       ! Sum the contribution
       ! Note that if we are here, then |Delta w_i| = 1
       contrib=contrib-0.5d0*pJ*Viiii
       
    enddo
    
!----------------------------------------------------------------------
! Exchange correction 1
!----------------------------------------------------------------------
! -pF/2 Sum_i Sum_j V_ijji , Delta w_i < 0, Delta w_j > 0
!----------------------------------------------------------------------
    ! Loop over negative Delta w_i values
    do i=1,ipos-1

       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       ! Loop over positive Delta w_j values
       do j=ipos,ndiff

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_j value
          Dwj=Dw(j,2)

          ! V_ijij
          Vijji=Vx(i1,j1)

          ! Sum the contribution (note that the plus sign arises
          ! because Dwi * Dwj < 0)
          contrib=contrib+0.5d0*pF*Vijji*Dwi*Dwj
                    
       enddo

    enddo

!----------------------------------------------------------------------
! Exchange correction 2
!----------------------------------------------------------------------
! -pF Sum_i Sum_j V_ijji <w omega|E_i^j E_j^i|w omega>,
! j > i, i and j singly-occupied
!----------------------------------------------------------------------
    ! Numbers of 'intermediate' CSFs entering into the contractions of
    ! the fibers of the spin-coupling coefficient tensor
    if (nopen > 1) then
       insp=ncsfs(nopen-2)
    else
       insp=0
    endif

    ! Loop over singly-occupied MOs (creation operator)
    do i=1,nsocc-1
       
       ! Creation operator index
       ic=socc(i)
       
       ! DFT/HF MO index
       i1=m2c(ic)

       ! Loop over singly-occupied MOs (annihilation operator)
       do j=i+1,nsocc
          
          ! Annihilation operator index
          ja=socc(j)
          
          ! DFT/HF MO index
          j1=m2c(ja)

          ! Get the spin coupling coefficient pattern index
          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)

          ! V_ijji
          Vijji=Vx(i1,j1)

          ! Sum the contributions
          start=pattern
          do omega=1,nsp
             product=dot_product(&
                  spincp(start:start+insp-1),&
                  spincp(start:start+insp-1))
             contrib(omega)=contrib(omega)-pF*Vijji*product
             start=start+insp
          enddo
          
       enddo

    enddo

!----------------------------------------------------------------------
! Exchange correction 3
!----------------------------------------------------------------------
! +pF/2 Sum_i Sum_j Vijji, Delta w_i > 0, Delta w_j > 0, and i indexes
! an open shell in the base configuration
!----------------------------------------------------------------------
    ! Loop over pairs of created MOs
    do i=ipos,ndiff

       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       ! Cycle if i does not index an open shell in the base
       ! configuration
       if (iopen0(i1) == 0) cycle
       
       do j=ipos,ndiff

          ! Cycle if j indexes an open shell in the base
          ! configuration
          if (i == j) cycle

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_j value
          Dwj=Dw(j,2)

          ! V_ijij
          Vijji=Vx(i1,j1)

          ! Sum the contribution
          contrib=contrib+0.5d0*pF*abs(Dwj)*Vijji
          
       enddo
    enddo

!----------------------------------------------------------------------
! Exchange correction 4
!----------------------------------------------------------------------
! +pF/2 Sum_i Sum_j Vijji, Delta w_i < 0, Delta w_j < 0, and i indexes
! an open shell in the base configuration
!----------------------------------------------------------------------
    ! Loop over pairs of annihilated MOs
    do i=1,ipos-1

       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       ! Cycle if i does not index an open shell in the base
       ! configuration
       if (iopen0(i1) == 0) cycle
       
       do j=1,ipos-1

          ! Cycle if j indexes an open shell in the base
          ! configuration
          if (i == j) cycle

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_j value
          Dwj=Dw(j,2)

          ! V_ijij
          Vijji=Vx(i1,j1)

          ! Sum the contribution
          contrib=contrib+0.5d0*pF*abs(Dwj)*Vijji
          
       enddo
    enddo

!----------------------------------------------------------------------
! Add the Coulomb and exchange corrections
!----------------------------------------------------------------------
    harr=harr+contrib
    
    return
    
  end subroutine hii_dftmrci_heil

!######################################################################
! hii_dftmrci_cvs: applies the Ottawa CVS-DFT/MRCI correction to a
!                  batch of on-diagonal Hamiltonian matrix elements
!######################################################################
  subroutine hii_dftmrci_cvs(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
       nsocc,nbefore)
    
    use constants
    use bitglobal
    use pattern_indices
    use hparam

    implicit none

    ! Array of on-diagonal Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: harr(nsp)

    ! Difference configuration information
    integer(is), intent(in) :: ndiff
    integer(is), intent(in) :: Dw(nmo,2)

    ! Number of open shells
    integer(is), intent(in) :: nopen
    
    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! SOP characterising the spatial occupation
    integer(ib), intent(in) :: sop(n_int,2)

    ! Indices of the singly-occupied MOs
    integer(is), intent(in) :: nsocc
    integer(is), intent(in) :: socc(nmo)

    ! Numbers of open shells preceding each MO
    integer(is), intent(in) :: nbefore(nmo)

    ! Temporary K-edge MO energy threshold
    real(dp), parameter     :: ethrsh=-10.0d0
    
    ! Everything else
    integer(is)             :: i,j,i1,j1,Dwi,Dwj,ipos,insp
    integer(is)             :: ic,ja,omega,pattern,start
    real(dp)                :: Viijj,Vijji,Viiii
    real(dp)                :: contrib(nsp)
    real(dp)                :: product
    real(dp)                :: pJ,pF,pJvv,pFvv,pJcv,pFcv

!----------------------------------------------------------------------
! Diagonal shift: 1/4 Sum_i V_iiii, i singly occupied in the base
! configuration    
!----------------------------------------------------------------------
    do i=1,nmo
       if (iopen0(i) == 1) harr=harr+0.25d0*Vc(i,i)
    enddo

!----------------------------------------------------------------------
! Return if we are at the base configuration
!----------------------------------------------------------------------
    if (ndiff == 0) return

!----------------------------------------------------------------------
! Parameter values
!----------------------------------------------------------------------
    ! Valence-valence interactions
    pJvv=hpar(1)
    pFvv=hpar(2)

    ! Core-valence interactions
    pJcv=hpar(5)
    pFcv=hpar(6)

!----------------------------------------------------------------------
! Sum_i F_ii^KS - F_ii^HF Delta w_i
!----------------------------------------------------------------------
    ! Loop over non-zero Delta w_i values
    do i=1,ndiff
    
       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
    
       ! Sum the contribution
       harr=harr+(moen(i1)-Fii(i1))*Dwi
       
    enddo

!----------------------------------------------------------------------
! Find the start of the postive Delta w_i values
!----------------------------------------------------------------------
    do i=1,ndiff
       if (Dw(i,2) > 0) then
          ipos=i
          exit
       endif
    enddo

!----------------------------------------------------------------------
! Coulomb correction 1
!----------------------------------------------------------------------
! -pJ Sum_i Sum_j Viijj, Delta w_i < 0, Delta w_j < 0
! -pJ Sum_i Sum_j Viijj, Delta w_i > 0, Delta w_j > 0
! +pJ Sum_i Sum_j Viijj, Delta w_i < 0, Delta w_j > 0
!----------------------------------------------------------------------
    contrib=0.0d0

    ! Loop over pairs of created/annihilated MOs (relative to the base
    ! configuration)
    do i=1,ndiff

       ! MO index
       i1=m2c(Dw(i,1))
       
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       do j=i,ndiff

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_i value
          Dwj=Dw(j,2)

          ! Are we at a valence-valence or core-valence interaction?
          if (moen(i1) < ethrsh .or. moen(j1) < ethrsh) then
             ! core-valence
             pJ=pJcv
          else
             ! valence-valence
             pJ=pJvv
          endif
          
          if (i == j .and. abs(Dwi) == 2) then
             ! Same MO index and double excitation

             ! V_iijj
             Viijj=Vc(i1,i1)

             ! Sum the contribution
             contrib=contrib-pJ*Viijj

          else if (i /= j) then
             ! Different MO indices

             ! V_iijj
             Viijj=Vc(i1,j1)

             ! Sum the contribution
             contrib=contrib-pJ*Viijj*Dwi*Dwj             

          endif
          
       enddo

    enddo

!----------------------------------------------------------------------
! Coulomb correction 2
!----------------------------------------------------------------------
! -pJ sum_i Viiii, |Delta w_i| = 1 and i indexes an open-shell in the
! base configuration
!----------------------------------------------------------------------
    ! Loop over created/annihilated MOs (relative to the base
    ! configuration)
    do i=1,ndiff

       ! MO index
       i1=m2c(Dw(i,1))
       
       ! Cycle if this MO is not singly-occupied in the base
       ! configuration
       if (iopen0(i1) == 0) cycle

       ! V_iiii
       Viiii=Vc(i1,i1)
       
       ! Sum the contribution
       ! Note that if we are here, then |Delta w_i| = 1
       contrib=contrib-0.5d0*pJ*Viiii
       
    enddo

!----------------------------------------------------------------------
! Exchange correction 1
!----------------------------------------------------------------------
! -pF/2 Sum_i Sum_j V_ijji , Delta w_i < 0, Delta w_j > 0
!----------------------------------------------------------------------
    ! Loop over negative Delta w_i values
    do i=1,ipos-1

       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       ! Loop over positive Delta w_j values
       do j=ipos,ndiff

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_j value
          Dwj=Dw(j,2)

          ! V_ijij
          Vijji=Vx(i1,j1)

          ! Are we at a valence-valence or core-valence interaction?
          if (moen(i1) < ethrsh .or. moen(j1) < ethrsh) then
             ! core-valence
             pF=pFcv
          else
             ! valence-valence
             pF=pFvv
          endif
          
          ! Sum the contribution (note that the plus sign arises
          ! because Dwi * Dwj < 0)
          contrib=contrib+0.5d0*pF*Vijji*Dwi*Dwj
                    
       enddo

    enddo

!----------------------------------------------------------------------
! Exchange correction 2
!----------------------------------------------------------------------
! -pF Sum_i Sum_j V_ijji <w omega|E_i^j E_j^i|w omega>,
! j > i, i and j singly-occupied
!----------------------------------------------------------------------
    ! Numbers of 'intermediate' CSFs entering into the contractions of
    ! the fibers of the spin-coupling coefficient tensor
    if (nopen > 1) then
       insp=ncsfs(nopen-2)
    else
       insp=0
    endif

    ! Loop over singly-occupied MOs (creation operator)
    do i=1,nsocc-1
       
       ! Creation operator index
       ic=socc(i)
       
       ! DFT/HF MO index
       i1=m2c(ic)

       ! Loop over singly-occupied MOs (annihilation operator)
       do j=i+1,nsocc
          
          ! Annihilation operator index
          ja=socc(j)
          
          ! DFT/HF MO index
          j1=m2c(ja)

          ! Get the spin coupling coefficient pattern index
          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)

          ! V_ijji
          Vijji=Vx(i1,j1)

          ! Are we at a valence-valence or core-valence interaction?
          if (moen(i1) < ethrsh .or. moen(j1) < ethrsh) then
             ! core-valence
             pF=pFcv
          else
             ! valence-valence
             pF=pFvv
          endif
          
          ! Sum the contributions
          start=pattern
          do omega=1,nsp
             product=dot_product(&
                  spincp(start:start+insp-1),&
                  spincp(start:start+insp-1))
             contrib(omega)=contrib(omega)-pF*Vijji*product
             start=start+insp
          enddo
          
       enddo

    enddo

!----------------------------------------------------------------------
! Exchange correction 3
!----------------------------------------------------------------------
! +pF/2 Sum_i Sum_j Vijji, Delta w_i > 0, Delta w_j > 0, and i indexes
! an open shell in the base configuration
!----------------------------------------------------------------------
    ! Loop over pairs of created MOs
    do i=ipos,ndiff

       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       ! Cycle if i does not index an open shell in the base
       ! configuration
       if (iopen0(i1) == 0) cycle
       
       do j=ipos,ndiff

          ! Cycle if j indexes an open shell in the base
          ! configuration
          if (i == j) cycle

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_j value
          Dwj=Dw(j,2)

          ! V_ijij
          Vijji=Vx(i1,j1)

          ! Are we at a valence-valence or core-valence interaction?
          if (moen(j1) < ethrsh) then
             ! core-valence
             pF=pFcv
          else
             ! valence-valence
             pF=pFvv
          endif
          
          ! Sum the contribution
          contrib=contrib+0.5d0*pF*abs(Dwj)*Vijji
          
       enddo
    enddo

!----------------------------------------------------------------------
! Exchange correction 4
!----------------------------------------------------------------------
! +pF/2 Sum_i Sum_j Vijji, Delta w_i < 0, Delta w_j < 0, and i indexes
! an open shell in the base configuration
!----------------------------------------------------------------------
    ! Loop over pairs of annihilated MOs
    do i=1,ipos-1

       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       ! Cycle if i does not index an open shell in the base
       ! configuration
       if (iopen0(i1) == 0) cycle
       
       do j=1,ipos-1

          ! Cycle if j indexes an open shell in the base
          ! configuration
          if (i == j) cycle

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_j value
          Dwj=Dw(j,2)

          ! V_ijij
          Vijji=Vx(i1,j1)

          ! Are we at a valence-valence or core-valence interaction?
          if (moen(j1) < ethrsh) then
             ! core-valence
             pF=pFcv
          else
             ! valence-valence
             pF=pFvv
          endif
          
          ! Sum the contribution
          contrib=contrib+0.5d0*pF*abs(Dwj)*Vijji
          
       enddo
    enddo

!----------------------------------------------------------------------
! Add the Coulomb and exchange corrections
!----------------------------------------------------------------------
    harr=harr+contrib
    
    return
    
  end subroutine hii_dftmrci_cvs
    
!######################################################################
! damping_grimme: for two CSF-averaged on-diagonal matrix element
!                 values, returns the value of Grimme's original
!                 DFT/MRCI damping function
!######################################################################
  function damping_grimme(av1,av2) result(func)

    use constants
    use bitglobal
    use hparam
    
    implicit none

    ! Function result
    real(dp)             :: func

    ! CSF-averaged on-diagonal matrix elements
    real(dp), intent(in) :: av1,av2

    ! Everything else
    real(dp)             :: DE4
    
    !
    !  p1 exp(-p2 DeltaE^4)
    !
    DE4=(av1-av2)**4
    func=hpar(1)*exp(-hpar(2)*DE4)
    
    return
    
  end function damping_grimme

!######################################################################
! damping_lyskov: for two CSF-averaged on-diagonal matrix element
!                 values, returns the value of Lyskov's redesigned
!                 DFT/MRCI damping function
!######################################################################
  function damping_lyskov(av1,av2) result(func)

    use constants
    use bitglobal
    use hparam
    
    implicit none

    ! Function result
    real(dp)             :: func

    ! CSF-averaged on-diagonal matrix elements
    real(dp), intent(in) :: av1,av2

    ! Everything else
    real(dp)             :: p2DE5

    !
    ! p1 / {1 + (p2 DeltaE)^5 arctan([p2 DeltaE]^5)}
    !
    p2DE5=hpar(4)*abs(av1-av2)
    p2DE5=p2DE5**5
    func=hpar(3)/(1.0d0+p2DE5*atan(p2DE5))

    return
    
  end function damping_lyskov

!######################################################################
! damping_heil18: for two CSF-averaged on-diagonal matrix element
!                 values, returns the value of Grimme's original
!                 DFT/MRCI damping function
!######################################################################
  function damping_heil18(av1,av2) result(func)

    use constants
    use bitglobal
    use hparam
    
    implicit none

    ! Function result
    real(dp)             :: func

    ! CSF-averaged on-diagonal matrix elements
    real(dp), intent(in) :: av1,av2

    ! Everything else
    real(dp)             :: DE6
    
    !
    !  p1 exp(-p2 DeltaE^6)
    !
    DE6=(av1-av2)**6
    func=hpar(3)*exp(-hpar(4)*DE6)
    
    return
    
  end function damping_heil18
  
!######################################################################
  
end module dftmrci
