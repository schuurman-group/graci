!**********************************************************************
! Routines for the construction of the MRCI and DFT/MRCI Hamiltonian
! matrix elements
!**********************************************************************
module hbuild_mrci

  implicit none

contains

!######################################################################
! hii_mrci: Computes a batch of on-diagonal Hamiltonian matrix
!           elements < w omega | H - E_SCF | w omega > for a single
!           spatial occupation w and all spin-couplings omega.
!######################################################################
  subroutine hii_mrci(harr,nsp,sop,nopen,socc,nsocc,docc,ndocc,&
       unocc,nunocc,nbefore,Dw,ndiff,m2c)

    use constants
    use bitglobal
    use pattern_indices
    use int_pyscf
    
    implicit none

    ! Number of spin-couplings (i.e., the no. CSFs)
    integer(is), intent(in) :: nsp

    ! SOP characterising the spatial occupation
    integer(ib), intent(in) :: sop(n_int,2)

    ! Number of open shells in the spatial configuration
    integer(is), intent(in) :: nopen

    ! Indices of the unoccupied, singly- and doubly-occupied
    ! MOs in the spatial configuration
    integer(is), intent(in) :: nsocc,ndocc,nunocc
    integer(is), intent(in) :: socc(nmo),docc(nmo),unocc(nmo)

    ! Numbers of open shells preceding each MO
    integer(is), intent(in) :: nbefore(nmo)

    ! Difference configuration information
    integer(is)             :: ndiff
    integer(is)             :: Dw(nmo,2)
    
    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)
    
    ! Hamiltonian matrix elements
    real(dp), intent(out)   :: harr(nsp)

    ! Everything else
    integer(is)             :: i,j,i1,j1,ic,ja,k
    integer(is)             :: Dwi,Dwj,wi,wj
    integer(is)             :: omega,indx
    integer(is)             :: nsp2b
    real(dp)                :: focksum,xsum1,xsum2
    real(dp)                :: product,Vijji

!----------------------------------------------------------------------
! Numbers of 'intermediate' CSFs entering into the contractions of the
! fibers of the spin-coupling coefficient tensor
!----------------------------------------------------------------------
    if (nopen > 1) then
       nsp2b=ncsfs(nopen-2)
    else
       nsp2b=0
    endif
    
!----------------------------------------------------------------------
! Common, spin-independent contribution
!----------------------------------------------------------------------
    !
    ! Initialisation
    !
    ! Fock-matrix contribution:
    ! Sum_i F_ii Delta w_i
    focksum=0.0d0
    ! Exchange integral contribution 1:
    ! Sum_ij (V_iijj - 1/2 V_ijji) Delta w_i Delta w_j
    xsum1=0.0d0
    ! Exchange integral contribution 2:
    ! Sum_ij V_ijji (1/2 w_i w_j - w_i)
    xsum2=0.0d0

    !
    ! Sum_i F_ii Delta w_i
    !
    ! Loop over non-zero Delta w_i values
    do i=1,ndiff

       ! MO index
       i1=m2c(Dw(i,1))

       ! Delta w_i value
       Dwi=Dw(i,2)
       
       ! Sum the contribution
       focksum=focksum+Fii(i1)*Dwi
       
    enddo

    !
    ! Sum_ij (V_iijj - 1/2 V_ijji) Delta w_i Delta w_j
    !
    do j=1,ndiff

       ! jth MO index
       j1=m2c(Dw(j,1))

       ! Delta w_j value
       Dwj=Dw(j,2)

       ! Sum the i=j contribution
       xsum1=xsum1+0.5d0*Vc(j1,j1)*Dwj**2
       
       do i=1,j-1
          
          ! ith MO index
          i1=m2c(Dw(i,1))

          ! Delta w_i value
          Dwi=Dw(i,2)

          ! Sum the i/=j contribution
          xsum1=xsum1+2.0d0*(Vc(i1,j1)-0.5d0*Vx(i1,j1))*Dwi*Dwj

       enddo

    enddo

    !
    ! Sum_ij V_ijji (1/2 w_i w_j - w_i)
    !
    ! Note that, due to a partial cancellation with terms
    ! V_ijji <w omega| E_i^j E_j^i | w omega>, only the
    ! contributions with i and j singly-occupied survive
    do i=1,nsocc
       i1=m2c(socc(i))
       xsum2=xsum2-0.5d0*Vc(i1,i1)
       do j=i+1,nsocc
          j1=m2c(socc(j))
          xsum2=xsum2-Vx(i1,j1)
       enddo
    enddo
    
    !
    ! Initialise all Hamiltonian matrix elements to the common,
    ! spin-independent contribution value
    !
    harr=focksum+0.5d0*(xsum1+xsum2)
    
!----------------------------------------------------------------------
! Case 2b spin-coupling coefficients
!----------------------------------------------------------------------
! Note that all other V_ijji <w omega| E_i^j E_j^i | w omega>
! terms cancel out with V_ijji (1/2 wi wj - wi) terms
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

          ! Get the spin coupling coefficient pattern index
          indx=pattern_index_case2b(sop,ic,ja,nbefore(ic),nbefore(ja),&
               nopen)
          
          ! V_ijji
          Vijji=Vx(i1,j1)
          
          ! Contributions to hii
          do omega=1,nsp
             product=dot_product(spincp2(1:nsp2b,omega,indx),&
                  spincp2(1:nsp2b,omega,indx))
             harr(omega)=harr(omega)+Vijji*product
          enddo
          
       enddo
       
    enddo

    return
    
  end subroutine hii_mrci

!######################################################################
! hij_mrci: Computes a batch of off-diagonal matrix elements
!           corresponding to a given pair of bra and ket configurations  
!######################################################################
  subroutine hij_mrci(harr,harrdim,nexci,bconf,kconf,bsop,ksop,&
       bnsp,knsp,bnopen,knopen,hlist,plist,m2c,&
       socc,nsocc,nbefore,Dw,ndiff,&
       bcsfs,kcsfs,bdim,kdim,bavii,kavii)

    use constants
    use bitglobal
    use mrci_integrals
    use dftmrci
    
    implicit none

    ! Array of off-diagonal matrix elements
    integer(is), intent(in) :: harrdim
    real(dp), intent(out)   :: harr(harrdim)

    ! Bra and ket configuration indices
    integer(is), intent(in) :: bconf,kconf
    
    ! Excitation degree
    integer(is), intent(in) :: nexci

    ! Bra and ket SOPs
    integer(ib), intent(in) :: bsop(n_int,2),ksop(n_int,2)

    ! Number of bra and ket CSFs
    integer(is), intent(in) :: bnsp,knsp
    
    ! Number of open shells in the bra and ket configurations
    integer(is), intent(in) :: bnopen,knopen

    ! Creation and annihilation operator indices
    integer(is), intent(in) :: hlist(2),plist(2)

    ! MO mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! Indices of the singly-occupied ket MOs
    integer(is), intent(in) :: socc(nmo)
    integer(is), intent(in) :: nsocc

    ! Number of open shells preceding each MO
    integer(is), intent(in) :: nbefore(nmo)
    
    ! Difference configuration information
    integer(is),intent(in)  :: ndiff
    integer(is),intent(in)  :: Dw(nmo,2)
    
    ! CSF offsets
    integer(is), intent(in) :: kdim,bdim
    integer(is), intent(in) :: bcsfs(bdim),kcsfs(kdim)

    ! Bra and ket spin-coupling averaged Hii values
    real(dp), intent(in)    :: bavii,kavii
    
    ! Spin-coupling sub-case bitsting encodings
    integer(ib)             :: pairindx(nmo)
    integer(ib)             :: icase
    
    ! Pattern indices
    integer(is)             :: bpattern(nmo+1),kpattern(nmo+1)

    ! No. CSFs for the intermediate configuration in the case
    ! of double excitations
    integer(is)             :: insp(2)
    
    ! Integrals
    real(dp)                :: Vpqrs(nmo)
    
!*********************************************************************
! Note that here the lists of holes/particles are the indices of the
! annihilation/creation operators operating on the ket configuration
! to yield the bra configuration
!**********************************************************************

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    harr(1:bnsp*knsp)=0.0d0
    
!----------------------------------------------------------------------
! Fill the integrals, pattern index, and pair index arrays
!----------------------------------------------------------------------
    select case(nexci)
    case(1)
       call package_integrals_nexci1(bsop,ksop,pairindx(1),hlist(1),&
            plist(1),bnopen,knopen,bpattern,kpattern,Vpqrs,m2c,&
            socc,nsocc,nbefore,Dw,ndiff,icase)
    case(2)
       !call package_integrals_nexci2(bsop,ksop,&
       !     pairindx(1:2),hlist(1:2),plist(1:2),bnopen,knopen,&
       !     bpattern(1:2),kpattern(1:2),Vpqrs(1:2),m2c,nbefore)

       call package_integrals_nexci2_new(bsop,ksop,hlist(1:2),&
            plist(1:2),bnopen,knopen,bpattern(1:2),kpattern(1:2),&
            Vpqrs(1:2),m2c,nbefore,insp)
       
    end select

!----------------------------------------------------------------------
! Compute the matrix elements between the CSFs generated by the two
! configurations
!----------------------------------------------------------------------
    select case(nexci)
    case(1)
       call hij_single_mrci_batch(bnopen,knopen,pairindx(1),icase,&
            bpattern,kpattern,Vpqrs,socc,nsocc,ndiff,hlist(1),&
            plist(1),harr,harrdim,bcsfs,kcsfs,bdim,kdim,bconf,kconf)
    case(2)
       !call hij_double_mrci_batch(bnopen,knopen,pairindx(1:2),&
       !     bpattern(1:2),kpattern(1:2),Vpqrs(1:2),plist(1:2),&
       !     hlist(1:2),harr,harrdim,bcsfs,kcsfs,bdim,kdim,bconf,kconf)

       call hij_double_mrci_batch_new(bnopen,knopen,bpattern(1:2),&
            kpattern(1:2),Vpqrs(1:2),plist(1:2),hlist(1:2),harr,&
            harrdim,bcsfs,kcsfs,bdim,kdim,bconf,kconf,insp)
       
    end select
    
!----------------------------------------------------------------------
! DFT/MRCI corrections
!----------------------------------------------------------------------
    if (ldftmrci) call hij_dftmrci_batch(harr(1:bnsp*knsp),bnsp,&
         knsp,bavii,kavii)

    return
    
  end subroutine hij_mrci
  
!######################################################################
! hij_single_mrci: Computes an off-diagonal Hamiltonian matrix element
!                  for a pair of CSFs differing by one pair of spatial
!                  orbital occupation. The spin couplings of the bra
!                  and ket CSFs are given by the indices bomega and
!                  komega, respectively.
!######################################################################
  function hij_single_mrci(bomega,komega,bnopen,knopen,pairindx,&
       icase,bpattern,kpattern,Vpqrs,socc,nsocc,ndiff,ia,ac) result(hij)

    use constants
    use bitglobal
    use bitstrings
    use mrci_integrals
    use iomod
    
    implicit none

    ! Function result
    real(dp)                :: hij

    ! Bra and ket spin couplings
    integer(is), intent(in) :: bomega,komega

    ! No. open shells in the bra and ket configurations
    integer(is), intent(in) :: bnopen,knopen

    ! Spin-coupling sub-case bitstring encodings
    integer(ib), intent(in) :: pairindx,icase

    ! Pattern indices
    integer(is), intent(in) :: bpattern(nmo+1),kpattern(nmo+1)

    ! Integrals and functions of integrals
    real(dp), intent(in)    :: Vpqrs(nmo)

    ! Singly-occupied MOs in the ket configuration
    integer(is), intent(in) :: nsocc
    integer(is), intent(in) :: socc(nmo)
    
    ! Number of excitations relative to the base configuration
    integer(is), intent(in) :: ndiff

    ! Indices of the annihilation and creation operators
    integer(is), intent(in) :: ia,ac
    
    ! Everything else
    integer(is)             :: indx,k,k1
    real(dp)                :: scp,halfscp,product

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    hij=0.0d0
    
!----------------------------------------------------------------------
! Get the spin-coupling coefficient <w' omega'|E_a^i|w omega>
!----------------------------------------------------------------------
    indx=kpattern(nsocc+1)
    select case(icase)
    case(i1a)
       scp=spincp1(komega,bomega,indx)
    case(i1b)
       scp=-spincp1(komega,bomega,indx)
    case(i2a)
       scp=spincp2(komega,bomega,indx)
    case(i2b)
       scp=spincp2(bomega,komega,indx)
    case default
       errmsg='Unrecognised icase value in hij_single_mrci'
       call error_control
    end select
       
!----------------------------------------------------------------------
! Sum_k V_ikka <w' omega'|E_a^k E_k^i - 1/2E_a^i|w omega>, k singly-
! occupied in the ket
!----------------------------------------------------------------------
    halfscp=0.5d0*scp
    
    ! Loop over singly-occupied MOs
    do k=1,nsocc

       ! MO index
       k1=socc(k)

       ! Cycle if the current MO corresponds to either the creation
       ! or annihilation operator
       if (k1 == ia) cycle
       if (k1 == ac) cycle

       ! Contraction of the fibers of the spin-coupling coefficient
       ! tensor
       product=contract_spincp(bomega,komega,bpattern(k),kpattern(k),&
            pairindx,bnopen,knopen)

       ! Sum the contribution
       hij=hij+Vpqrs(k)*(product-halfscp)
       
    enddo

!----------------------------------------------------------------------
! [F_ia + Sum_k (V_iakk - 1/2 V_ikka) Delta w_k]
! x <w' omega'|E_a^i|w omega>
!----------------------------------------------------------------------
    hij=hij+Vpqrs(nsocc+1)*scp
    
!----------------------------------------------------------------------
! 1/2 [V_aaai w_a + Vaiii (w_i -2)] <w' omega'|E_a^i|w omega>
!----------------------------------------------------------------------
    hij=hij+Vpqrs(nsocc+2)*scp
    
    return
    
  end function hij_single_mrci

!######################################################################
! hij_single_mrci_batch: Computes a batch of off-diagonal Hamiltonian
!                        matrix element for a pair of configurations
!                        differing by one pair of spatial orbital
!                        occupation.
!######################################################################
  subroutine hij_single_mrci_batch(bnopen,knopen,pairindx,icase,&
       bpattern,kpattern,Vpqrs,socc,nsocc,ndiff,ia,ac,harr,harrdim,&
       bcsfs,kcsfs,bdim,kdim,bconf,kconf)

    use constants
    use bitglobal
    use bitstrings
    use mrci_integrals
    use iomod
    
    implicit none

    ! No. open shells in the bra and ket configurations
    integer(is), intent(in) :: bnopen,knopen

    ! Spin-coupling sub-case bitstring encodings
    integer(ib), intent(in) :: pairindx,icase

    ! Pattern indices
    integer(is), intent(in) :: bpattern(nmo+1),kpattern(nmo+1)

    ! Integrals and functions of integrals
    real(dp), intent(in)    :: Vpqrs(nmo)

    ! Singly-occupied MOs in the ket configuration
    integer(is), intent(in) :: nsocc
    integer(is), intent(in) :: socc(nmo)
    
    ! Number of excitations relative to the base configuration
    integer(is), intent(in) :: ndiff

    ! Indices of the annihilation and creation operators
    integer(is), intent(in) :: ia,ac

    ! Array of off-diagonal matrix elements
    integer(is), intent(in) :: harrdim
    real(dp), intent(out)   :: harr(harrdim)

    ! CSF offsets
    integer(is), intent(in) :: kdim,bdim
    integer(is), intent(in) :: bcsfs(bdim),kcsfs(kdim)

    ! Bra and ket configuration indices
    integer(is), intent(in) :: bconf,kconf
    
    ! Everything else
    integer(is)             :: bomega,komega,bcsf,kcsf
    integer(is)             :: indx,k,k1,counter
    real(dp)                :: scp,halfscp,product

    ! Initialise the harr counter
    counter=0
    
    ! Loop over ket CSFs
    komega=0
    do kcsf=kcsfs(kconf),kcsfs(kconf+1)-1
       
       ! Ket CSF spin coupling index
       komega=komega+1
       
       ! Loop over bra CSFs
       bomega=0
       do bcsf=bcsfs(bconf),bcsfs(bconf+1)-1

          ! Increment the harr counter
          counter=counter+1
          
          ! Bra CSF spin coupling index
          bomega=bomega+1

          ! Get the spin-coupling coefficient <w' omega'|E_a^i|w omega>
          indx=kpattern(nsocc+1)
          select case(icase)
          case(i1a)
             scp=spincp1(komega,bomega,indx)
          case(i1b)
             scp=-spincp1(komega,bomega,indx)
          case(i2a)
             scp=spincp2(komega,bomega,indx)
          case(i2b)
             scp=spincp2(bomega,komega,indx)
          case default
             errmsg='Unrecognised icase value in hij_single_mrci'
             call error_control
          end select

          !
          ! Sum_k V_ikka <w' omega'|E_a^k E_k^i - 1/2E_a^i|w omega>,
          ! k singly-occupied in the ket
          !          
          halfscp=0.5d0*scp
          
          ! Loop over singly-occupied MOs
          do k=1,nsocc
             
             ! MO index
             k1=socc(k)
             
             ! Cycle if the current MO corresponds to either the creation
             ! or annihilation operator
             if (k1 == ia) cycle
             if (k1 == ac) cycle
             
             ! Contraction of the fibers of the spin-coupling coefficient
             ! tensor
             product=contract_spincp(bomega,komega,bpattern(k),&
                  kpattern(k),pairindx,bnopen,knopen)
             
             ! Sum the contribution
             harr(counter)=harr(counter)+Vpqrs(k)*(product-halfscp)
       
          enddo

          !
          ! [F_ia + Sum_k (V_iakk - 1/2 V_ikka) Delta w_k]
          ! x <w' omega'|E_a^i|w omega>
          !
          harr(counter)=harr(counter)+Vpqrs(nsocc+1)*scp
          
          !
          ! 1/2 [V_aaai w_a + Vaiii (w_i -2)] <w' omega'|E_a^i|w omega>
          !
          harr(counter)=harr(counter)+Vpqrs(nsocc+2)*scp
          
       enddo

    enddo
       
    return
    
  end subroutine hij_single_mrci_batch

!######################################################################
! hij_double_mrci: Computes an off-diagonal Hamiltonian matrix element
!                  for a pair of CSFs differing by two pairs of spatial
!                  orbital occupations. The spin couplings of the bra
!                  and ket CSFs are given by the indices bomega and
!                  komega, respectively.
!######################################################################
  function hij_double_mrci(bomega,komega,bnopen,knopen,pairindx,&
       bpattern,kpattern,Vpqrs,plist,hlist) result(hij)

    use constants
    use bitglobal
    use mrci_integrals
    
    implicit none

    ! Function result
    real(dp)                :: hij

    ! Bra and ket spin couplings
    integer(is), intent(in) :: bomega,komega

    ! No. open shells in the bra and ket configurations
    integer(is), intent(in) :: bnopen,knopen

    ! Spin-coupling sub-case pair bitstrings
    integer(ib), intent(in) :: pairindx(2)

    ! Pattern indices
    integer(is), intent(in) :: bpattern(2),kpattern(2)

    ! Integrals (pre-scaled by the 1/[(1+delta_ab)*(1+delta_ij)]
    ! prefactor)
    real(dp), intent(in)    :: Vpqrs(2)

    ! Indices of the creation and annihilation operators
    integer(is), intent(in) :: plist(2),hlist(2)
    integer(is)             :: ic1,ic2,ja1,ja2
    
    ! Everything else
    real(dp)                :: product

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    hij=0.0d0

!----------------------------------------------------------------------
! Creation and annihilation operator indices
!----------------------------------------------------------------------
    ic1=plist(1)
    ja1=hlist(1)

    ic2=plist(2)
    ja2=hlist(2)
    
!----------------------------------------------------------------------
! V_aibj contribution    
!----------------------------------------------------------------------
    ! Contraction of the fibers of the spin-coupling coefficient
    ! tensor
    product=contract_spincp(bomega,komega,bpattern(1),kpattern(1),&
         pairindx(1),bnopen,knopen)

    ! Contribution to hij
    hij=hij+Vpqrs(1)*product

!----------------------------------------------------------------------
! V_ajbi contribution    
!----------------------------------------------------------------------
    ! Contraction of the fibers of the spin-coupling coefficient
    ! tensor
    product=contract_spincp(bomega,komega,bpattern(2),kpattern(2),&
         pairindx(2),bnopen,knopen)

    ! Contribution to hij
    hij=hij+Vpqrs(2)*product
    
    return
    
  end function hij_double_mrci
  
!######################################################################
! hij_double_mrci_batch: Computes a batch of off-diagonal Hamiltonian
!                        matrix element for a pair of configurations
!                        differing by two pairs of spatial orbital
!                        occupations.
!######################################################################
  subroutine hij_double_mrci_batch(bnopen,knopen,pairindx,&
       bpattern,kpattern,Vpqrs,plist,hlist,harr,harrdim,&
       bcsfs,kcsfs,bdim,kdim,bconf,kconf)

    use constants
    use bitglobal
    use mrci_integrals
    
    implicit none

    ! No. open shells in the bra and ket configurations
    integer(is), intent(in) :: bnopen,knopen

    ! Spin-coupling sub-case pair bitstrings
    integer(ib), intent(in) :: pairindx(2)

    ! Pattern indices
    integer(is), intent(in) :: bpattern(2),kpattern(2)

    ! Integrals (pre-scaled by the 1/[(1+delta_ab)*(1+delta_ij)]
    ! prefactor)
    real(dp), intent(in)    :: Vpqrs(2)

    ! Indices of the creation and annihilation operators
    integer(is), intent(in) :: plist(2),hlist(2)
    integer(is)             :: ic1,ic2,ja1,ja2

    ! Array of off-diagonal matrix elements
    integer(is), intent(in) :: harrdim
    real(dp), intent(out)   :: harr(harrdim)

    ! CSF offsets
    integer(is), intent(in) :: kdim,bdim
    integer(is), intent(in) :: bcsfs(bdim),kcsfs(kdim)

    ! Bra and ket configuration indices
    integer(is), intent(in) :: bconf,kconf
    
    ! CSFs
    integer(is)             :: bomega,komega,bcsf,kcsf
    
    ! Everything else
    integer(is)             :: counter
    real(dp)                :: product

!----------------------------------------------------------------------
! Creation and annihilation operator indices
!----------------------------------------------------------------------
    ic1=plist(1)
    ja1=hlist(1)

    ic2=plist(2)
    ja2=hlist(2)
    
!----------------------------------------------------------------------
! Compute the matrix elements
!----------------------------------------------------------------------
    ! Initialise the harr counter
    counter=0

    ! Loop over ket CSFs
    komega=0
    do kcsf=kcsfs(kconf),kcsfs(kconf+1)-1
       
       ! Ket CSF spin coupling index
       komega=komega+1
    
       ! Loop over bra CSFs
       bomega=0
       do bcsf=bcsfs(bconf),bcsfs(bconf+1)-1

          ! Increment the harr counter
          counter=counter+1
          
          ! Bra CSF spin coupling index
          bomega=bomega+1

          ! V_aibj contribution    
          product=contract_spincp(bomega,komega,bpattern(1),&
               kpattern(1),pairindx(1),bnopen,knopen)
          harr(counter)=harr(counter)+Vpqrs(1)*product
          
          ! V_ajbi contribution    
          product=contract_spincp(bomega,komega,bpattern(2),&
               kpattern(2),pairindx(2),bnopen,knopen)
          harr(counter)=harr(counter)+Vpqrs(2)*product
          
       enddo

    enddo
       
    return
    
  end subroutine hij_double_mrci_batch

!######################################################################
! hij_double_mrci_batch: Computes a batch of off-diagonal Hamiltonian
!                        matrix element for a pair of configurations
!                        differing by two pairs of spatial orbital
!                        occupations.
!######################################################################
  subroutine hij_double_mrci_batch_new(bnopen,knopen,bpattern,&
       kpattern,Vpqrs,plist,hlist,harr,harrdim,bcsfs,kcsfs,bdim,kdim,&
       bconf,kconf,insp)

    use constants
    use bitglobal
    use mrci_integrals
    
    implicit none

    ! No. open shells in the bra and ket configurations
    integer(is), intent(in) :: bnopen,knopen

    ! Pattern indices
    integer(is), intent(in) :: bpattern(2),kpattern(2)

    ! Integrals (pre-scaled by the 1/[(1+delta_ab)*(1+delta_ij)]
    ! prefactor)
    real(dp), intent(in)    :: Vpqrs(2)

    ! Indices of the creation and annihilation operators
    integer(is), intent(in) :: plist(2),hlist(2)
    integer(is)             :: ic1,ic2,ja1,ja2

    ! Array of off-diagonal matrix elements
    integer(is), intent(in) :: harrdim
    real(dp), intent(out)   :: harr(harrdim)

    ! CSF offsets
    integer(is), intent(in) :: kdim,bdim
    integer(is), intent(in) :: bcsfs(bdim),kcsfs(kdim)

    ! Bra and ket configuration indices
    integer(is), intent(in) :: bconf,kconf

    ! Number of CSFs for the intermediate configuration obtained
    ! by acting on the ket CSF with the first singlet excitation
    ! operator
    integer(is), intent(in) :: insp(2)
    
    ! CSFs
    integer(is)             :: bomega,komega,bcsf,kcsf,bnsp,knsp
    
    ! Everything else
    integer(is)             :: counter,bstart1,bstart2,kstart1,kstart2
    real(dp)                :: product

!----------------------------------------------------------------------
! Creation and annihilation operator indices
!----------------------------------------------------------------------
    ic1=plist(1)
    ja1=hlist(1)

    ic2=plist(2)
    ja2=hlist(2)

!----------------------------------------------------------------------
! Number of bra and ket CSFs
!----------------------------------------------------------------------
    bnsp=ncsfs(bnopen)
    knsp=ncsfs(knopen)

!----------------------------------------------------------------------
! Compute the matrix elements
!----------------------------------------------------------------------
    ! Initialise counters
    kstart1=kpattern(1)
    kstart2=kpattern(2)
    komega=0
    counter=0

    ! Loop over ket CSFs
    do kcsf=kcsfs(kconf),kcsfs(kconf+1)-1
       
       ! Ket CSF spin coupling index
       komega=komega+1

       ! Loop over bra CSFs
       bstart1=bpattern(1)
       bstart2=bpattern(2)
       bomega=0
       do bcsf=bcsfs(bconf),bcsfs(bconf+1)-1

          ! Increment the harr counter
          counter=counter+1
          
          ! Bra CSF spin coupling index
          bomega=bomega+1

          ! V_aibj contribution
          harr(counter)=harr(counter)+Vpqrs(1) &
               *dot_product(&
               spincp(bstart1:bstart1+insp(1)-1),&
               spincp(kstart1:kstart1+insp(1)-1))
          
          ! V_ajbi contribution    
          harr(counter)=harr(counter)+Vpqrs(2) &
               *dot_product(&
               spincp(bstart2:bstart2+insp(2)-1),&
               spincp(kstart2:kstart2+insp(2)-1))
          
          ! Update the bra starting point in the spincp array
          bstart1=bstart1+insp(1)
          bstart2=bstart2+insp(2)
          
       enddo
          
       ! Update the ket starting point in the spincp array
       kstart1=kstart1+insp(1)
       kstart2=kstart2+insp(2)
       
    enddo
    
    return
    
  end subroutine hij_double_mrci_batch_new
  
!######################################################################
! hij_same_mrci: Computes a batch of off-diagonal Hamiltonian matrix
!                elements < w omega' | H - E_SCF | w omega > for CSFs
!                with the spatial occupation w but different
!                spin couplings omega' and omega
!######################################################################
  subroutine hij_same_mrci(harr,arrdim,sop,socc,nsocc,nbefore,m2c)

    use constants
    use bitglobal
    use pattern_indices
    use dftmrci
    
    implicit none

    ! SOP
    integer(ib), intent(in) :: sop(n_int,2)

    ! Singly-occupied MOs
    integer(is), intent(in) :: socc(nmo)
    integer(is), intent(in) :: nsocc

    ! Numbers of open shells preceding each MO
    integer(is), intent(in) :: nbefore(nmo)

    ! Hamiltonian matrix elements
    integer(is), intent(in) :: arrdim
    real(dp), intent(out)   :: harr(arrdim)

    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)
    
    ! Everything else
    integer(is)             :: nsp,nsp2b,nopen
    integer(is)             :: i,i1,j,j1,ic,ja
    integer(is)             :: bomega,komega,indx
    integer(is)             :: count,n
    real(dp)                :: Vijji,product

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    harr=0.0d0
    
!----------------------------------------------------------------------
! Number of CSFs
!----------------------------------------------------------------------
    nopen=nsocc
    nsp=ncsfs(nopen)
    
!----------------------------------------------------------------------
! Numbers of 'intermediate' CSFs entering into the contractions of the
! fibers of the spin-coupling coefficient tensor
!----------------------------------------------------------------------
    if (nopen > 1) then
       nsp2b=ncsfs(nopen-2)
    else
       nsp2b=0
    endif

!----------------------------------------------------------------------
! Case 2b spin-coupling coefficients
!----------------------------------------------------------------------
! Note that all other <w omega| E_i^j E_j^i | w omega> terms are zero
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

          ! Get the spin coupling coefficient pattern index
          indx=pattern_index_case2b(sop,ic,ja,nbefore(ic),nbefore(ja),&
               nopen)
          
          ! V_ijji
          Vijji=Vx(i1,j1)
          
          ! Contributions to hij
          count=0
          do bomega=1,nsp-1
             do komega=bomega+1,nsp
                count=count+1
                product=dot_product(spincp2(1:nsp2b,bomega,indx),&
                     spincp2(1:nsp2b,komega,indx))
                harr(count)=harr(count)+Vijji*product
             enddo
          enddo

       enddo

    enddo

!----------------------------------------------------------------------
! DFT/MRCI corrections
!----------------------------------------------------------------------
    if (ldftmrci) call hij_same_dftmrci(harr,nsp)
       
    return
    
  end subroutine hij_same_mrci

!######################################################################
! hij_same_mrci_1element: Computes a single off-diagonal Hamiltonian
!                         matrix element
!                         < w omega' | H - E_SCF | w omega >
!                         between bra and ket CSFs with the spatial
!                         occupation w but different spin couplings
!                         omega' and omega
!######################################################################
  function hij_same_mrci_1element(bomega,komega,sop,socc,nsocc,&
       nbefore,m2c) result(hij)

    use constants
    use bitglobal
    use pattern_indices
    
    implicit none

    ! Function result
    real(dp)                :: hij
    
    ! Bra and ket spin-coupling indices
    integer(is), intent(in) :: bomega,komega
    
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
    integer(is)             :: nsp,nsp2b,nopen
    integer(is)             :: i,i1,j,j1,ic,ja
    integer(is)             :: indx,count
    real(dp)                :: Vijji,product

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    hij=0.0d0
    
!----------------------------------------------------------------------
! Number of CSFs
!----------------------------------------------------------------------
    nopen=nsocc
    nsp=ncsfs(nopen)
    
!----------------------------------------------------------------------
! Numbers of 'intermediate' CSFs entering into the contractions of the
! fibers of the spin-coupling coefficient tensor
!----------------------------------------------------------------------
    if (nopen > 1) then
       nsp2b=ncsfs(nopen-2)
    else
       nsp2b=0
    endif

!----------------------------------------------------------------------
! Case 2b spin-coupling coefficients
!----------------------------------------------------------------------
! Note that all other <w omega| E_i^j E_j^i | w omega> terms are zero
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

          ! Get the spin coupling coefficient pattern index
          indx=pattern_index_case2b(sop,ic,ja,nbefore(ic),nbefore(ja),&
               nopen)
          
          ! V_ijji
          Vijji=Vx(i1,j1)
          
          ! Contributions to hij
          product=dot_product(spincp2(1:nsp2b,bomega,indx),&
               spincp2(1:nsp2b,komega,indx))
          hij=hij+Vijji*product

       enddo

    enddo

    return
    
  end function hij_same_mrci_1element
  
!######################################################################
  
end module hbuild_mrci
