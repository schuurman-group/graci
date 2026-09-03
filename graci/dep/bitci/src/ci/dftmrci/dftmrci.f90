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
    
    ! MO index mapping arrays
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

    case(6:9,11:12)
       ! Heil's parameterisations
       ! Also used by the the QE8 parameterisation
       call hii_dftmrci_heil(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
            nsocc,nbefore)

    case(13,14)
       ! CVS-QE8 parameterisation
       call hii_dftmrci_cvsqe8(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
            nsocc,nbefore)

    case(10)
       ! R2022 parameterisation
       call hii_dftmrci_r2022(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
            nsocc,nbefore)

    case(15)
       ! 2026 parameterisation
       call hii_dftmrci_2026(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
            nsocc,nbefore)

    case(16)
       ! CVS-2026 parameterisation
       call hii_dftmrci_cvs2026(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
            nsocc,nbefore)

    case(17)
       ! QE8 with symmetry-averaged exchange
       call hii_dftmrci_heil_sym(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
            nsocc,nbefore)

    case(18,19,20,21)
       ! RC DFT/MRCI: range-corrected exchange. ihamiltonian=19 and 20
       ! are restricted forms; hii_dftmrci_rc remaps their parameter
       ! indices.
       call hii_dftmrci_rc(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
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
  subroutine hij_dftmrci_batch(hij,bdim,kdim,bav,kav,hij_lr)

    use constants
    use bitglobal
    use hparam
    use iomod
    use lrsr_diag

    implicit none

    ! Hamiltonian matrix elements
    integer(is), intent(in)           :: bdim,kdim
    real(dp), intent(inout)           :: hij(:)
    real(dp), intent(in)              :: bav,kav
    real(dp), optional, intent(in)    :: hij_lr(:)
    real(dp)                          :: damp,damp_lr

!----------------------------------------------------------------------
! Compute the damping factor
!----------------------------------------------------------------------
    select case(ihamiltonian)

    case(2:3)
       ! Grimme's parameterisation
       damp=damping_grimme(bav,kav)

    case(4:7)
       ! Lyskov's parameterisation
       ! Note that this is also used for Heil's 2017 Hamiltonian
       damp=damping_lyskov(bav,kav)

    case(8:9)
       ! Heil's 2018 parameterisation
       damp=damping_heil18(bav,kav)

    case(10)
       ! R2022 parameterisation
       damp=damping_r2022(bav,kav)

    case(11:14,17)
       ! QE8 parameterisations
       damp=damping_qe8(bav,kav)

    case(18)
       ! RC DFT/MRCI parameterisation
       damp=damping_rc(bav,kav)
       damp_lr=damping_rc_lr(bav,kav)

    case(19)
       ! Reduced form: the SR contribution is damped as usual and the LR
       ! contribution is taken at its ab initio value
       damp=damping_rc5p(bav,kav)
       damp_lr=1.0d0

    case(20)
       ! Physically-motivated form: two prefactors, one shared damping
       ! function. p1_SR=hpar(4), p1_LR=hpar(5), p2=hpar(6), n=hpar(7).
       damp=damping_rc7p(bav,kav,hpar(4))
       damp_lr=damping_rc7p(bav,kav,hpar(5))

    case(21)
       ! Split Coulomb and exchange, shared damping shape.
       ! p1_SR=hpar(5), p1_LR=hpar(6), p2=hpar(7), n=hpar(8).
       damp=damping_rc8p(bav,kav,hpar(5))
       damp_lr=damping_rc8p(bav,kav,hpar(6))

    case(15:16)
       ! 2026 and CVS-2026 parameterisations
       damp=damping_2026(bav,kav)

    case default
       errmsg='Your Hamiltonian choice has not been implemented yet'
       call error_control

    end select

!----------------------------------------------------------------------
! Apply the damping factor
!----------------------------------------------------------------------
    if (ihamiltonian >= 18 .and. ihamiltonian <= 21 &
         .and. present(hij_lr)) then
       ! LR/SR split, each with its own damping function:
       !
       !   H_IJ = damp * H_SR + damp_lr * H_LR
       !
       ! hij holds H_full and hij_lr holds H_LR, so H_SR = hij - hij_lr and
       ! the combination rearranges to the form below. Writing it this way
       ! avoids forming H_SR as a temporary.
       !
       ! The previous version used (1 - damp) as the LR coefficient, which
       ! is algebraically damp*H_SR + H_LR: the LR term was undamped at
       ! every gap, and hpar(5) was computed but never applied.
       !
       ! Bin the undamped SR and LR contributions first: once the line
       ! below runs, H_SR is no longer recoverable from hij.
       if (lrsr_active) &
            call lrsr_accumulate(hij(1:bdim*kdim),hij_lr(1:bdim*kdim),&
            bdim*kdim,bav,kav,damp,damp_lr)
       hij(1:bdim*kdim)=damp*hij(1:bdim*kdim) &
            +(damp_lr-damp)*hij_lr(1:bdim*kdim)
    else
       hij(1:bdim*kdim)=hij(1:bdim*kdim)*damp
    endif

    return

  end subroutine hij_dftmrci_batch

!######################################################################
! hij_same_dftmrci: applies DFT/MRCI corrections to a batch of
!                   off-diagonal Hamiltonian matrix elements with the
!                   same spatial part but different spin couplings
!######################################################################
  subroutine hij_same_dftmrci(hij,nsp,Dw,ndiff,sop,socc,nsocc,nbefore,&
       m2c)

    use constants
    use bitglobal
    use hparam
    use iomod
    
    implicit none

    ! Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: hij(:)

    ! Difference configuration information
    integer(is), intent(in) :: ndiff
    integer(is), intent(in) :: Dw(nmo,2)
    
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

    case(4:9,11:12,14,17)
       ! Lyskov's parameterisation
       ! Note that this is also used for Heil's Hamiltonian,
       ! the QE8 Hamiltonian, and qe8_sym.
       ! ihamiltonian == 14 is cvs-test and uses same exchange damping
       !                    for valence-valence and core-valence
       nij=nsp*(nsp-1)/2
       hij(1:nij)=(1.0d0-hpar(2))*hij(1:nij)
       return

    case(18,21)
       ! RC DFT/MRCI: same-config off-diagonal uses (1-pF_SR).
       ! ihamiltonian=21 carries pK_SR in the same slot.
       nij=nsp*(nsp-1)/2
       hij(1:nij)=(1.0d0-hpar(3))*hij(1:nij)
       return

    case(19,20)
       ! Restricted forms: pK_SR is hpar(2) here, not hpar(3)
       nij=nsp*(nsp-1)/2
       hij(1:nij)=(1.0d0-hpar(2))*hij(1:nij)
       return

    case(15)
       ! 2026 parameterisation: distinguish he-type vs hhee-type pairs
       call hij_same_dftmrci_2026(hij,nsp,Dw,ndiff,sop,socc,nsocc,&
            nbefore,m2c)

    case(16)
       ! CVS-2026 parameterisation
       call hij_same_dftmrci_cvs2026(hij,nsp,sop,socc,nsocc,nbefore,m2c)

    case(10)
       ! R2022 parameterisation
       call hij_same_dftmrci_r2022(hij,nsp,Dw,ndiff,sop,socc,nsocc,&
            nbefore,m2c)
       
    case(13)
       ! CVS-QE8 parameterisation
       call hij_same_dftmrci_cvsqe8(hij,nsp,sop,socc,nsocc,nbefore,m2c)
       
    case default
       errmsg='Your Hamiltonian choice has not been implemented yet'
       call error_control
       
    end select
    
    return
    
  end subroutine hij_same_dftmrci

!######################################################################
! hij_same_dftmrci_cvs: applies the CVS-QE8 corrections to a batch of
!                       off-diagonal Hamiltonian matrix elements with
!                       the same spatial part but different spin
!                       couplings
!######################################################################
  subroutine hij_same_dftmrci_cvsqe8(hij,nsp,sop,socc,nsocc,nbefore,m2c)

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
    if (nhpar == 6) then
       pFcv=hpar(6)
    else if (nhpar == 7) then
       pFcv=hpar(7)
    endif

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
          if (icvs(i1) == 1 .or. icvs(j1) == 1) then
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

  end subroutine hij_same_dftmrci_cvsqe8

!######################################################################
! hij_same_dftmrci_cvs2026: applies the CVS-2026 DFT/MRCI corrections
!                           to off-diagonal Hamiltonian matrix elements
!                           with the same spatial part but different
!                           spin couplings.
!                           Uses pF_he_vv for valence-valence pairs and
!                           pFcv for core-valence pairs.
!######################################################################
  subroutine hij_same_dftmrci_cvs2026(hij,nsp,sop,socc,nsocc,nbefore,m2c)

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

    ! Everything else
    integer(is)             :: insp,nopen
    integer(is)             :: i,i1,j,j1,ic,ja
    integer(is)             :: bomega,komega
    integer(is)             :: pattern,kstart,bstart
    integer(is)             :: count,n
    real(dp)                :: Vijji,product
    real(dp)                :: pF,pF_he_vv,pFcv

!----------------------------------------------------------------------
! Parameter values
!----------------------------------------------------------------------
    pF_he_vv=hpar(2)
    pFcv=hpar(7)

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
    do i=1,nsocc-1

       ic=socc(i)
       i1=m2c(ic)

       do j=i+1,nsocc

          ja=socc(j)
          j1=m2c(ja)

          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)

          Vijji=Vx(i1,j1)

          if (icvs(i1) == 1 .or. icvs(j1) == 1) then
             pF=pFcv
          else
             pF=pF_he_vv
          endif

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

  end subroutine hij_same_dftmrci_cvs2026

!######################################################################
! hij_same_dftmrci_r2022: applies the R2022 DFT/MRCI corrections to
!                         a batch of off-diagonal Hamiltonian matrix
!                         elements with the same spatial part but
!                         different spin couplings
!######################################################################
  subroutine hij_same_dftmrci_r2022(hij,nsp,Dw,ndiff,sop,socc,nsocc,&
       nbefore,m2c)

    use constants
    use bitglobal
    use pattern_indices
    use hparam

    implicit none

    ! Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: hij(:)

    ! Difference configuration information
    integer(is), intent(in) :: ndiff
    integer(is), intent(in) :: Dw(nmo,2)
    
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
    integer(is)             :: insp,nopen
    integer(is)             :: i,i1,j,j1,ic,ja
    integer(is)             :: bomega,komega
    integer(is)             :: pattern,kstart,bstart
    integer(is)             :: count,n
    integer(is)             :: Dwi_open(nomax)
    integer(is)             :: Dwi,Dwj
    real(dp)                :: Vijji,product
    real(dp)                :: px_he,px_hhee,px
    
!----------------------------------------------------------------------
! Parameter values
!----------------------------------------------------------------------
    pX_he=hpar(4)
    pX_hhee=hpar(5)

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
! Determine the Delta w_i values for the open shells
!----------------------------------------------------------------------
    Dwi_open=0
    
    ! Loop over created/annihilated MOs
    do i=1,ndiff
       i1=Dw(i,1)

       ! Loop over open shells
       do j=1,nsocc
          j1=socc(j)

          ! Does the created/annihilated MO match an open shell
          ! index?
          if (i1 == j1) Dwi_open(j)=Dw(i,2)

       enddo
    enddo
    
!----------------------------------------------------------------------
! Case 2b spin-coupling coefficients
!----------------------------------------------------------------------
! Note that all other <w omega' | E_i^j E_j^i | w omega> terms are zero
!----------------------------------------------------------------------
    ! Loop over singly-occupied MOs (creation operator)
    do i=1,nsocc-1

       ! Cycle if Delta w_i = 0
       Dwi=Dwi_open(i)
       
       ! Creation operator index
       ic=socc(i)
       
       ! DFT/HF MO index
       i1=m2c(ic)

       ! Loop over singly-occupied MOs (annihilation operator)
       do j=i+1,nsocc

          ! Cycle if Delta w_j = 0
          Dwj=Dwi_open(j)
          
          ! Annihilation operator index
          ja=socc(j)
          
          ! DFT/HF MO index
          j1=m2c(ja)

          ! Get the spin coupling coefficient pattern indices
          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)
          
          ! V_ijji
          Vijji=Vx(i1,j1)

          ! Exchange scaling parameter
          if (Dwi*Dwj <= 0) then
             px=px_he
          else
             px=px_hhee
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
                   hij(count)=hij(count)-px*Vijji*product
                endif
                bstart=bstart+insp
             enddo
             kstart=kstart+insp
          enddo
          
       enddo

    enddo
    
    return
    
  end subroutine hij_same_dftmrci_r2022

!######################################################################
! hij_same_dftmrci_2026: applies the R2026 DFT/MRCI corrections to
!                        a batch of off-diagonal Hamiltonian matrix
!                        elements with the same spatial part but
!                        different spin couplings.
!
!                        he-type pairs  (Dwi*Dwj <= 0): scale by pF_he   = hpar(2)
!                        hhee-type pairs (Dwi*Dwj > 0): scale by pF_hhee = hpar(3)
!######################################################################
  subroutine hij_same_dftmrci_2026(hij,nsp,Dw,ndiff,sop,socc,nsocc,&
       nbefore,m2c)

    use constants
    use bitglobal
    use pattern_indices
    use hparam

    implicit none

    ! Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: hij(:)

    ! Difference configuration information
    integer(is), intent(in) :: ndiff
    integer(is), intent(in) :: Dw(nmo,2)

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
    integer(is)             :: insp,nopen
    integer(is)             :: i,i1,j,j1,ic,ja
    integer(is)             :: bomega,komega
    integer(is)             :: pattern,kstart,bstart
    integer(is)             :: count,n
    integer(is)             :: Dwi_open(nomax)
    integer(is)             :: Dwi,Dwj
    real(dp)                :: Vijji,product
    real(dp)                :: pF_he,pF_hhee,px

!----------------------------------------------------------------------
! Parameter values: hpar(2)=pF_he, hpar(3)=pF_hhee
!----------------------------------------------------------------------
    pF_he  =hpar(2)
    pF_hhee=hpar(3)

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
! Determine the Delta w_i values for the open shells
!----------------------------------------------------------------------
    Dwi_open=0

    ! Loop over created/annihilated MOs
    do i=1,ndiff
       i1=Dw(i,1)

       ! Loop over open shells
       do j=1,nsocc
          j1=socc(j)

          ! Does the created/annihilated MO match an open shell index?
          if (i1 == j1) Dwi_open(j)=Dw(i,2)

       enddo
    enddo

!----------------------------------------------------------------------
! Case 2b spin-coupling coefficients
!----------------------------------------------------------------------
! Note that all other <w omega' | E_i^j E_j^i | w omega> terms are zero
!----------------------------------------------------------------------
    ! Loop over singly-occupied MOs (creation operator)
    do i=1,nsocc-1

       ! Cycle if Delta w_i = 0
       Dwi=Dwi_open(i)

       ! Creation operator index
       ic=socc(i)

       ! DFT/HF MO index
       i1=m2c(ic)

       ! Loop over singly-occupied MOs (annihilation operator)
       do j=i+1,nsocc

          ! Cycle if Delta w_j = 0
          Dwj=Dwi_open(j)

          ! Annihilation operator index
          ja=socc(j)

          ! DFT/HF MO index
          j1=m2c(ja)

          ! Get the spin coupling coefficient pattern indices
          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)

          ! V_ijji
          Vijji=Vx(i1,j1)

          ! Exchange scaling parameter: he-type vs hhee-type
          if (Dwi*Dwj <= 0) then
             px=pF_he
          else
             px=pF_hhee
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
                   hij(count)=hij(count)-px*Vijji*product
                endif
                bstart=bstart+insp
             enddo
             kstart=kstart+insp
          enddo

       enddo

    enddo

    return

  end subroutine hij_same_dftmrci_2026

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
! hii_dftmrci_heil_sym: QE8/Heil diagonal with symmetry-averaged
!                       exchange (symvx) to enforce degeneracy of
!                       Abelian-subgroup components of degenerate states
!######################################################################
  subroutine hii_dftmrci_heil_sym(harr,nsp,Dw,ndiff,nopen,m2c,sop,&
       socc,nsocc,nbefore)

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
    do i=1,ndiff
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       harr=harr+(moen(i1)-Fii(i1))*Dwi
    enddo

!----------------------------------------------------------------------
! Find the start of the positive Delta w_i values
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

    do i=1,ndiff
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       do j=i,ndiff
          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)
          if (i == j .and. abs(Dwi) == 2) then
             Viijj=Vc(i1,i1)
             contrib=contrib-pJ*Viijj
          else if (i /= j) then
             Viijj=Vc(i1,j1)
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
    do i=1,ndiff
       i1=m2c(Dw(i,1))
       if (iopen0(i1) == 0) cycle
       Viiii=Vc(i1,i1)
       contrib=contrib-0.5d0*pJ*Viiii
    enddo

!----------------------------------------------------------------------
! Exchange correction 1
!----------------------------------------------------------------------
! -pF/2 Sum_i Sum_j V_ijji , Delta w_i < 0, Delta w_j > 0
!----------------------------------------------------------------------
    do i=1,ipos-1
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       do j=ipos,ndiff
          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)
          Vijji=symvx(i1,j1)
          contrib=contrib+0.5d0*pF*Vijji*Dwi*Dwj
       enddo
    enddo

!----------------------------------------------------------------------
! Exchange correction 2
!----------------------------------------------------------------------
! -pF Sum_i Sum_j V_ijji <w omega|E_i^j E_j^i|w omega>,
! j > i, i and j singly-occupied
!----------------------------------------------------------------------
    if (nopen > 1) then
       insp=ncsfs(nopen-2)
    else
       insp=0
    endif

    do i=1,nsocc-1
       ic=socc(i)
       i1=m2c(ic)
       do j=i+1,nsocc
          ja=socc(j)
          j1=m2c(ja)
          if (degen_orbs(i1) == j1) cycle
          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)
          Vijji=symvx(i1,j1)
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
    do i=ipos,ndiff
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       if (iopen0(i1) == 0) cycle
       do j=ipos,ndiff
          if (i == j) cycle
          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)
          Vijji=symvx(i1,j1)
          contrib=contrib+0.5d0*pF*abs(Dwj)*Vijji
       enddo
    enddo

!----------------------------------------------------------------------
! Exchange correction 4
!----------------------------------------------------------------------
! +pF/2 Sum_i Sum_j Vijji, Delta w_i < 0, Delta w_j < 0, and i indexes
! an open shell in the base configuration
!----------------------------------------------------------------------
    do i=1,ipos-1
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       if (iopen0(i1) == 0) cycle
       do j=1,ipos-1
          if (i == j) cycle
          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)
          Vijji=symvx(i1,j1)
          contrib=contrib+0.5d0*pF*abs(Dwj)*Vijji
       enddo
    enddo

!----------------------------------------------------------------------
! Add the Coulomb and exchange corrections
!----------------------------------------------------------------------
    harr=harr+contrib

    return

  end subroutine hii_dftmrci_heil_sym

!######################################################################
! hii_dftmrci_rc: RC DFT/MRCI diagonal corrections.
!                 Based on hii_dftmrci_heil_sym with range-separated
!                 exchange: EC1-EC4 use
!                   pF_SR * symvx(i,j) + (pF_LR - pF_SR) * v_lr(i,j)
!                 EC2 constraint (skip degenerate pairs) is applied.
!######################################################################
  subroutine hii_dftmrci_rc(harr,nsp,Dw,ndiff,nopen,m2c,sop,&
       socc,nsocc,nbefore)

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
    real(dp)                :: Viijj,Vijji_eff,Viiii
    real(dp)                :: contrib(nsp)
    real(dp)                :: product
    real(dp)                :: pJSR,pJLR,dJLR,pFSR,pFLR,dFLR
    real(dp)                :: Viijj_eff

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
! hpar(1)=pJ_SR, hpar(2)=pJ_LR, hpar(3)=pF_SR,  hpar(4)=pF_LR,
! hpar(5)=p1_LR, hpar(6)=p2_LR, hpar(7)=n_LR,
! hpar(8)=p1,    hpar(9)=p2,    hpar(10)=n
!----------------------------------------------------------------------
    if (ihamiltonian == 19 .or. ihamiltonian == 20) then
       ! Restricted forms: a single Coulomb scaling, so dJLR = 0 and the LR
       ! Coulomb term drops out of the diagonal entirely. Exchange keeps
       ! its SR/LR split, which the range-separated functional justifies.
       pJSR = hpar(1)
       pJLR = hpar(1)
       pFSR = hpar(2)
       pFLR = hpar(3)
    else
       pJSR = hpar(1)
       pJLR = hpar(2)
       pFSR = hpar(3)
       pFLR = hpar(4)
    endif
    dJLR = pJLR - pJSR   ! increment applied to LR Coulomb
    dFLR = pFLR - pFSR   ! increment applied to LR exchange

!----------------------------------------------------------------------
! Sum_i F_ii^KS - F_ii^HF Delta w_i
!----------------------------------------------------------------------
    do i=1,ndiff
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       harr=harr+(moen(i1)-Fii(i1))*Dwi
    enddo

!----------------------------------------------------------------------
! Find the start of the positive Delta w_i values
!----------------------------------------------------------------------
    do i=1,ndiff
       if (Dw(i,2) > 0) then
          ipos=i
          exit
       endif
    enddo

!----------------------------------------------------------------------
! Coulomb correction 1
! V_iijj_eff = pJ_SR * Vc(i,j) + dJLR * j_lr(i,j)
!----------------------------------------------------------------------
    contrib=0.0d0

    do i=1,ndiff
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       do j=i,ndiff
          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)
          if (i == j .and. abs(Dwi) == 2) then
             Viijj_eff=pJSR*Vc(i1,i1)+dJLR*Vc_lr(i1,i1)
             contrib=contrib-Viijj_eff
          else if (i /= j) then
             Viijj_eff=pJSR*Vc(i1,j1)+dJLR*Vc_lr(i1,j1)
             contrib=contrib-Viijj_eff*Dwi*Dwj
          endif
       enddo
    enddo

!----------------------------------------------------------------------
! Coulomb correction 2
!----------------------------------------------------------------------
    do i=1,ndiff
       i1=m2c(Dw(i,1))
       if (iopen0(i1) == 0) cycle
       Viijj_eff=pJSR*Vc(i1,i1)+dJLR*Vc_lr(i1,i1)
       contrib=contrib-0.5d0*Viijj_eff
    enddo

!----------------------------------------------------------------------
! Exchange correction 1
! -pF/2 Sum_i Sum_j V_ijji, Delta w_i < 0, Delta w_j > 0
! V_ijji_eff = pF_SR * symvx(i,j) + dFLR * v_lr(i,j)
!----------------------------------------------------------------------
    do i=1,ipos-1
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       do j=ipos,ndiff
          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)
          Vijji_eff = pFSR*symvx(i1,j1) + dFLR*symvx_lr(i1,j1)
          contrib=contrib+0.5d0*Vijji_eff*Dwi*Dwj
       enddo
    enddo

!----------------------------------------------------------------------
! Exchange correction 2
! -pF Sum_i Sum_j V_ijji <w omega|E_i^j E_j^i|w omega>
! EC2 degeneracy constraint: skip pairs where i and j are degenerate
!----------------------------------------------------------------------
    if (nopen > 1) then
       insp=ncsfs(nopen-2)
    else
       insp=0
    endif

    do i=1,nsocc-1
       ic=socc(i)
       i1=m2c(ic)
       do j=i+1,nsocc
          ja=socc(j)
          j1=m2c(ja)
          if (degen_orbs(i1) == j1) cycle
          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)
          Vijji_eff = pFSR*symvx(i1,j1) + dFLR*symvx_lr(i1,j1)
          start=pattern
          do omega=1,nsp
             product=dot_product(&
                  spincp(start:start+insp-1),&
                  spincp(start:start+insp-1))
             contrib(omega)=contrib(omega)-Vijji_eff*product
             start=start+insp
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Exchange correction 3
! +pF/2 Sum_i Sum_j Vijji, Delta w_i > 0, Delta w_j > 0, i open in base
!----------------------------------------------------------------------
    do i=ipos,ndiff
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       if (iopen0(i1) == 0) cycle
       do j=ipos,ndiff
          if (i == j) cycle
          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)
          Vijji_eff = pFSR*symvx(i1,j1) + dFLR*symvx_lr(i1,j1)
          contrib=contrib+0.5d0*Vijji_eff*abs(Dwj)
       enddo
    enddo

!----------------------------------------------------------------------
! Exchange correction 4
! +pF/2 Sum_i Sum_j Vijji, Delta w_i < 0, Delta w_j < 0, i open in base
!----------------------------------------------------------------------
    do i=1,ipos-1
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       if (iopen0(i1) == 0) cycle
       do j=1,ipos-1
          if (i == j) cycle
          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)
          Vijji_eff = pFSR*symvx(i1,j1) + dFLR*symvx_lr(i1,j1)
          contrib=contrib+0.5d0*Vijji_eff*abs(Dwj)
       enddo
    enddo

!----------------------------------------------------------------------
! Add the Coulomb and exchange corrections
!----------------------------------------------------------------------
    harr=harr+contrib

    return

  end subroutine hii_dftmrci_rc

!######################################################################
! hii_dftmrci_2026: on-diagonal DFT/MRCI corrections for the 2026
!                   Hamiltonian. Based on the Heil functional form but
!                   with separate exchange scaling for hole-electron
!                   (he, pF_he=hpar(2)) and hole-hole/electron-electron
!                   (hhee, pF_hhee=hpar(3)) orbital pairs.
!######################################################################
  subroutine hii_dftmrci_2026(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
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
    real(dp)                :: pJ,pF_he,pF_hhee

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
    pF_he=hpar(2)
    pF_hhee=hpar(3)

!----------------------------------------------------------------------
! Sum_i (F_ii^KS - F_ii^HF) Delta w_i
!----------------------------------------------------------------------
    do i=1,ndiff
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       harr=harr+(moen(i1)-Fii(i1))*Dwi
    enddo

!----------------------------------------------------------------------
! Find the start of the positive Delta w_i values
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
! -pJ Sum_i<j Viijj Delta w_i Delta w_j  (all pairs)
!----------------------------------------------------------------------
    contrib=0.0d0

    do i=1,ndiff

       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)

       do j=i,ndiff

          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)

          if (i == j .and. abs(Dwi) == 2) then
             Viijj=Vc(i1,i1)
             contrib=contrib-pJ*Viijj

          else if (i /= j) then
             Viijj=Vc(i1,j1)
             contrib=contrib-pJ*Viijj*Dwi*Dwj

          endif

       enddo

    enddo

!----------------------------------------------------------------------
! Coulomb correction 2
!----------------------------------------------------------------------
! -pJ Sum_i Viiii, |Delta w_i| = 1 and i indexes an open-shell in the
! base configuration
!----------------------------------------------------------------------
    do i=1,ndiff

       i1=m2c(Dw(i,1))
       if (iopen0(i1) == 0) cycle

       Viiii=Vc(i1,i1)
       contrib=contrib-0.5d0*pJ*Viiii

    enddo

!----------------------------------------------------------------------
! Exchange correction 1
!----------------------------------------------------------------------
! -pF_he/2 Sum_i<j V_ijji Delta w_i Delta w_j,
! Delta w_i < 0, Delta w_j > 0  (hole-electron pairs)
!----------------------------------------------------------------------
    do i=1,ipos-1

       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)

       do j=ipos,ndiff

          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)

          Vijji=symvx(i1,j1)
          contrib=contrib+0.5d0*pF_he*Vijji*Dwi*Dwj

       enddo

    enddo

!----------------------------------------------------------------------
! Exchange correction 2
!----------------------------------------------------------------------
! -pF_he Sum_i<j V_ijji <w omega|E_i^j E_j^i|w omega>,
! i and j both singly-occupied in the current configuration
! (spin-coupling-dependent exchange for open-shell pairs)
!----------------------------------------------------------------------
    if (nopen > 1) then
       insp=ncsfs(nopen-2)
    else
       insp=0
    endif

    do i=1,nsocc-1

       ic=socc(i)
       i1=m2c(ic)

       do j=i+1,nsocc

          ja=socc(j)
          j1=m2c(ja)

          if (degen_orbs(i1) == j1) cycle

          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)

          Vijji=symvx(i1,j1)

          start=pattern
          do omega=1,nsp
             product=dot_product(&
                  spincp(start:start+insp-1),&
                  spincp(start:start+insp-1))
             contrib(omega)=contrib(omega)-pF_he*Vijji*product
             start=start+insp
          enddo

       enddo

    enddo

!----------------------------------------------------------------------
! Exchange correction 3
!----------------------------------------------------------------------
! +pF_hhee/2 Sum_i Sum_j V_ijji |Delta w_j|,
! Delta w_i > 0, Delta w_j > 0, i indexes a base open shell
! (electron-electron pairs)
!----------------------------------------------------------------------
    do i=ipos,ndiff

       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)

       if (iopen0(i1) == 0) cycle

       do j=ipos,ndiff

          if (i == j) cycle

          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)

          Vijji=symvx(i1,j1)
          contrib=contrib+0.5d0*pF_hhee*abs(Dwj)*Vijji

       enddo
    enddo

!----------------------------------------------------------------------
! Exchange correction 4
!----------------------------------------------------------------------
! +pF_hhee/2 Sum_i Sum_j V_ijji |Delta w_j|,
! Delta w_i < 0, Delta w_j < 0, i indexes a base open shell
! (hole-hole pairs)
!----------------------------------------------------------------------
    do i=1,ipos-1

       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)

       if (iopen0(i1) == 0) cycle

       do j=1,ipos-1

          if (i == j) cycle

          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)

          Vijji=symvx(i1,j1)
          contrib=contrib+0.5d0*pF_hhee*abs(Dwj)*Vijji

       enddo
    enddo

!----------------------------------------------------------------------
! Add the Coulomb and exchange corrections
!----------------------------------------------------------------------
    harr=harr+contrib

    return

  end subroutine hii_dftmrci_2026

!######################################################################
! hii_dftmrci_cvs: applies the CVS-QE8 correction to a batch of
! on-diagonal Hamiltonian matrix elements
!######################################################################
  subroutine hii_dftmrci_cvsqe8(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
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
    ! ihamiltonian == 14 is for cvs-test -- eventually should be standardized
    if (ihamiltonian == 14) then
       pJcv=hpar(6)
       pFcv=hpar(2)
    else if (nhpar == 6) then
       pJcv=hpar(5)
       pFcv=hpar(6)
    else if (nhpar == 7) then
       pJcv=hpar(6)
       pFcv=hpar(7)
    endif

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
          if (icvs(i1) == 1 .or. icvs(j1) == 1) then
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
    pJ = pJvv
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
          if (icvs(i1) == 1 .or. icvs(j1) == 1) then
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
          if (icvs(i1) == 1 .or. icvs(j1) == 1) then
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
          if (icvs(j1) == 1) then
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
          if (icvs(j1) == 1) then
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

  end subroutine hii_dftmrci_cvsqe8

!######################################################################
! hii_dftmrci_cvs2026: on-diagonal DFT/MRCI corrections for the
!                      CVS-2026 Hamiltonian. As hii_dftmrci_2026 but
!                      with separate core-valence Coulomb (pJ_cv) and
!                      exchange (pF_cv) scaling.
!######################################################################
  subroutine hii_dftmrci_cvs2026(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
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
    real(dp)                :: pJvv,pF_he_vv,pF_hhee_vv,pJcv,pFcv

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
    pJvv=hpar(1)
    pF_he_vv=hpar(2)
    pF_hhee_vv=hpar(3)
    ! hpar(4) = p2, hpar(5) = n  (used in damping only)
    pJcv=hpar(6)
    pFcv=hpar(7)

!----------------------------------------------------------------------
! Sum_i (F_ii^KS - F_ii^HF) Delta w_i
!----------------------------------------------------------------------
    do i=1,ndiff
       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)
       harr=harr+(moen(i1)-Fii(i1))*Dwi
    enddo

!----------------------------------------------------------------------
! Find the start of the positive Delta w_i values
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
! -pJ Sum_i<j Viijj Delta w_i Delta w_j  (all pairs)
! pJ = pJcv if either MO is a core orbital, else pJvv
!----------------------------------------------------------------------
    contrib=0.0d0

    do i=1,ndiff

       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)

       do j=i,ndiff

          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)

          if (icvs(i1) == 1 .or. icvs(j1) == 1) then
             pJ=pJcv
          else
             pJ=pJvv
          endif

          if (i == j .and. abs(Dwi) == 2) then
             Viijj=Vc(i1,i1)
             contrib=contrib-pJ*Viijj

          else if (i /= j) then
             Viijj=Vc(i1,j1)
             contrib=contrib-pJ*Viijj*Dwi*Dwj

          endif

       enddo

    enddo

!----------------------------------------------------------------------
! Coulomb correction 2
!----------------------------------------------------------------------
! -pJvv Sum_i Viiii, |Delta w_i| = 1 and i indexes an open-shell in
! the base configuration (base open shells are always valence)
!----------------------------------------------------------------------
    do i=1,ndiff

       i1=m2c(Dw(i,1))
       if (iopen0(i1) == 0) cycle

       Viiii=Vc(i1,i1)
       contrib=contrib-0.5d0*pJvv*Viiii

    enddo

!----------------------------------------------------------------------
! Exchange correction 1
!----------------------------------------------------------------------
! -pF_he/2 Sum_i<j V_ijji Delta w_i Delta w_j,
! Delta w_i < 0, Delta w_j > 0  (hole-electron pairs)
! pF = pFcv if either MO is a core orbital, else pF_he_vv
!----------------------------------------------------------------------
    do i=1,ipos-1

       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)

       do j=ipos,ndiff

          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)

          if (icvs(i1) == 1 .or. icvs(j1) == 1) then
             pF=pFcv
          else
             pF=pF_he_vv
          endif

          Vijji=Vx(i1,j1)
          contrib=contrib+0.5d0*pF*Vijji*Dwi*Dwj

       enddo

    enddo

!----------------------------------------------------------------------
! Exchange correction 2
!----------------------------------------------------------------------
! -pF Sum_i<j V_ijji <w omega|E_i^j E_j^i|w omega>,
! i and j both singly-occupied in the current configuration
! pF = pFcv if either MO is a core orbital, else pF_he_vv
!----------------------------------------------------------------------
    if (nopen > 1) then
       insp=ncsfs(nopen-2)
    else
       insp=0
    endif

    do i=1,nsocc-1

       ic=socc(i)
       i1=m2c(ic)

       do j=i+1,nsocc

          ja=socc(j)
          j1=m2c(ja)

          if (degen_orbs(i1) == j1) cycle

          if (icvs(i1) == 1 .or. icvs(j1) == 1) then
             pF=pFcv
          else
             pF=pF_he_vv
          endif

          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)

          Vijji=Vx(i1,j1)

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
! +pF/2 Sum_i Sum_j V_ijji |Delta w_j|,
! Delta w_i > 0, Delta w_j > 0, i indexes a base open shell
! (electron-electron pairs; i is always valence; pF = pFcv if j is core)
!----------------------------------------------------------------------
    do i=ipos,ndiff

       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)

       if (iopen0(i1) == 0) cycle

       do j=ipos,ndiff

          if (i == j) cycle

          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)

          if (icvs(j1) == 1) then
             pF=pFcv
          else
             pF=pF_hhee_vv
          endif

          Vijji=Vx(i1,j1)
          contrib=contrib+0.5d0*pF*abs(Dwj)*Vijji

       enddo
    enddo

!----------------------------------------------------------------------
! Exchange correction 4
!----------------------------------------------------------------------
! +pF/2 Sum_i Sum_j V_ijji |Delta w_j|,
! Delta w_i < 0, Delta w_j < 0, i indexes a base open shell
! (hole-hole pairs; i is always valence; pF = pFcv if j is core)
!----------------------------------------------------------------------
    do i=1,ipos-1

       i1=m2c(Dw(i,1))
       Dwi=Dw(i,2)

       if (iopen0(i1) == 0) cycle

       do j=1,ipos-1

          if (i == j) cycle

          j1=m2c(Dw(j,1))
          Dwj=Dw(j,2)

          if (icvs(j1) == 1) then
             pF=pFcv
          else
             pF=pF_hhee_vv
          endif

          Vijji=Vx(i1,j1)
          contrib=contrib+0.5d0*pF*abs(Dwj)*Vijji

       enddo
    enddo

!----------------------------------------------------------------------
! Add the Coulomb and exchange corrections
!----------------------------------------------------------------------
    harr=harr+contrib

    return

  end subroutine hii_dftmrci_cvs2026

!######################################################################
! hii_dftmrci_r2022: applies the R2022 DFT/MRCI correction to a batch
!                    of on-diagonal Hamiltonian matrix elements
!######################################################################
  subroutine hii_dftmrci_r2022(harr,nsp,Dw,ndiff,nopen,m2c,sop,socc,&
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
    integer(is)             :: Dwi_open(nomax)
    real(dp)                :: Viijj,Vijji,Viiii
    real(dp)                :: contrib(nsp)
    real(dp)                :: product
    real(dp)                :: pJ_he,pJ_hhee,pJ_eeee,pX_hhee,pX_he
    real(dp)                :: px
    
!----------------------------------------------------------------------
! Parameter values
!----------------------------------------------------------------------
    pJ_he=hpar(2)
    pJ_hhee=hpar(3)
    pJ_eeee=hpar(3)
    pX_he=hpar(4)
    pX_hhee=hpar(5)

!----------------------------------------------------------------------
! Return here if we are at the base configuration
!----------------------------------------------------------------------
    if (ndiff == 0) return
    
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
! Determine the Delta w_i values for the open shells
!----------------------------------------------------------------------
    Dwi_open=0

    ! Loop over created/annihilated MOs
    do i=1,ndiff
       i1=Dw(i,1)

       ! Loop over open shells
       do j=1,nsocc
          j1=socc(j)

          ! Does the created/annihilated MO match an open shell
          ! index?
          if (i1 == j1) Dwi_open(j)=Dw(i,2)

       enddo
    enddo
    
!----------------------------------------------------------------------
! Term (1)
!----------------------------------------------------------------------
! Sum_i (F_ii^KS - F_ii^HF) Delta w_i
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
! Term (2a)
!----------------------------------------------------------------------
! - Sum_i<j (pJ_he Viijj -1/2 pX_he Vijji) Delta w_i Delta w_j,
! Delta w_i * Delta w_j < 0
!----------------------------------------------------------------------
    contrib=0.0d0
    
    ! Loop over hole MOs
    do i=1,ipos-1

       ! Hole MO index
       i1=m2c(Dw(i,1))

       ! Hole Delta w_i value
       Dwi=Dw(i,2)
       
       ! Loop over particle MOs
       do j=ipos,ndiff

          ! Particle MO index
          j1=m2c(Dw(j,1))

          ! Particle Delta w_j value
          Dwj=Dw(j,2)

          ! Integral values
          Viijj=Vc(i1,j1)
          Vijji=Vx(i1,j1)

          ! Sum the contribution
          contrib=contrib-(pJ_he*Viijj-0.5d0*pX_he*Vijji)*Dwi*Dwj
          
       enddo
          
    enddo
    
!----------------------------------------------------------------------
! Term (2b)
!----------------------------------------------------------------------
! - Sum_i<j (pJ_he Viijj -1/2 pX_he Vijji) Delta w_i Delta w_j,
! Delta w_i * Delta w_j > 0
!----------------------------------------------------------------------
! Hole-hole terms
!----------------------------------------------------------------------
    ! Loop over unique pairs of hole MOs
    do i=1,ipos-2
       do j=i+1,ipos-1

          ! Hole MO indices
          i1=m2c(Dw(i,1))
          j1=m2c(Dw(j,1))

          ! Delta w values
          Dwi=Dw(i,2)
          Dwj=Dw(j,2)

          ! Integral values
          Viijj=Vc(i1,j1)
          Vijji=Vx(i1,j1)

          ! Sum the contribution
          contrib=contrib-(pJ_hhee*Viijj-0.5d0*pX_hhee*Vijji)*Dwi*Dwj
          
       enddo

    enddo

!----------------------------------------------------------------------
! Term (2c)
!----------------------------------------------------------------------
! - Sum_i<j (pJ_he Viijj -1/2 pX_he Vijji) Delta w_i Delta w_j,
! Delta w_i * Delta w_j > 0
!----------------------------------------------------------------------
! Particle-particle terms
!----------------------------------------------------------------------
    ! Loop over unique pairs of particle MOs
    do i=ipos,ndiff-1
       do j=i+1,ndiff

          ! Particle MO indices
          i1=m2c(Dw(i,1))
          j1=m2c(Dw(j,1))

          ! Delta w values
          Dwi=Dw(i,2)
          Dwj=Dw(j,2)

          ! Integral values
          Viijj=Vc(i1,j1)
          Vijji=Vx(i1,j1)

          ! Sum the contribution
          contrib=contrib-(pJ_hhee*Viijj-0.5d0*pX_hhee*Vijji)*Dwi*Dwj
          
       enddo

    enddo

!----------------------------------------------------------------------
! Terms (3a) and (3b)
!----------------------------------------------------------------------
! - Sum_i<j pX Vijji (<w omega| E_i^j E_j^i |w omega> - 1/2),
! where both i and j index an open shell in the configuration w, and
! pX = pX_he or pX_hhee depending on the sign of
! Delta w_i * Delta w_j 
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

       ! Cycle if Delta w_i = 0
       Dwi=Dwi_open(i)
       
       ! Creation operator index
       ic=socc(i)
       
       ! DFT/HF MO index
       i1=m2c(ic)
       
       ! Loop over singly-occupied MOs (annihilation operator)
       do j=i+1,nsocc

          ! Cycle if Delta w_j = 0
          Dwj=Dwi_open(j)
          
          ! Annihilation operator index
          ja=socc(j)
          
          ! DFT/HF MO index
          j1=m2c(ja)

          ! Get the spin coupling coefficient pattern index
          pattern=pattern_index_case2b(sop,ic,ja,nbefore(ic),&
               nbefore(ja),nopen)

          ! V_ijji
          Vijji=Vx(i1,j1)

          ! Determine the value of the exchange scaling parameter
          ! for this open shell pair
          if (Dwi*Dwj <= 0) then
             pX=pX_he
          else
             pX=pX_hhee
          endif

          ! Sum the contributions
          start=pattern
          do omega=1,nsp
             product=dot_product(&
                  spincp(start:start+insp-1),&
                  spincp(start:start+insp-1))
             contrib(omega)=contrib(omega)-pX*Vijji*(product-0.5d0)
             start=start+insp
          enddo
          
       enddo

    enddo

!----------------------------------------------------------------------
! Term (4)
!----------------------------------------------------------------------
! -1/4 pJ_eeee Sum_i V_iiii (Delta w_i)^2, for MOs i that are *not*
! open shells in the base configuration
!----------------------------------------------------------------------
    ! Loop over created/annihilated MOs
    do i=1,ndiff

       ! MO index
       i1=m2c(Dw(i,1))

       ! Delta w_i value
       Dwi=Dw(i,2)

       ! Cycle if the MO is an open shell in the base conf
       if (iopen0(i1) == 1) cycle

       ! V_iiii
       Viiii=Vc(i1,i1)

       ! Sum the contribution
       contrib=contrib-0.25d0*pJ_eeee*Viiii*Dwi**2
       
    enddo

!!----------------------------------------------------------------------
!! Term (5)
!----------------------------------------------------------------------
! -1/4 pJ_he Sum_i V_iiii |Delta w_i|, for MOs i that *are*
! open shells in the base configuration
!----------------------------------------------------------------------
    ! Loop over created/annihilated MOs
    do i=1,ndiff

       ! MO index
       i1=m2c(Dw(i,1))

       ! Delta w_i value
       Dwi=Dw(i,2)

       ! Cycle if the MO is not an open shell in the base conf
       if (iopen0(i1) == 0) cycle

       ! V_iiii
       Viiii=Vc(i1,i1)

       ! Sum the contribution
       contrib=contrib-0.25d0*pJ_he*Viiii*abs(Dwi)
       
    enddo
    
!----------------------------------------------------------------------
! Term (6)
!----------------------------------------------------------------------
! 1/4 pJ_eeee Sum_i V_iiii, for MOs i that are *not*
! open shells in the base configuration and are singly-occupied in the
! configuration w
!----------------------------------------------------------------------
    ! Loop over singly-occupied MOs
    do i=1,nsocc

       ! MO index
       i1=m2c(socc(i))

       ! Cycle if this MO is an open shell in the base conf
       if (iopen0(i1) == 1) cycle

       ! V_iiii
       Viiii=Vc(i1,i1)

       ! Sum the contribution
       contrib=contrib+0.25d0*pJ_eeee*Viiii

    enddo

!----------------------------------------------------------------------
! Term (7)
!----------------------------------------------------------------------
! 1/2 pX_he Sum_i<j V_ijji Delta w_i Delta w_j,
! Delta w_i * Delta w_j < 0 and i,j not singly-occupied in the base
! conf
!----------------------------------------------------------------------
    ! Loop over hole MOs
    do i=1,ipos-1

       ! Hole MO index
       i1=m2c(Dw(i,1))

       ! Hole Delta w_i value
       Dwi=Dw(i,2)

       ! Cycle if the hole MO is an open shell in the base conf
       if (iopen0(i1) == 1) cycle
       
       ! Loop over particle MOs
       do j=ipos,ndiff

          ! Particle MO index
          j1=m2c(Dw(j,1))

          ! Particle Delta w_j value
          Dwj=Dw(j,2)

          ! Cycle if the particle MO is an open shell in the base conf
          if (iopen0(j1) == 1) cycle

          ! V_ijji
          Vijji=Vx(i1,j1)

          ! Sum the contribution
          contrib=contrib+0.5d0*pX_he*Vijji*Dwi*Dwj
          
       enddo

    enddo

!----------------------------------------------------------------------
! Add the Coulomb and exchange corrections
!----------------------------------------------------------------------
    harr=harr+contrib
    
    return
    
  end subroutine hii_dftmrci_r2022
  
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
!                 values, returns the value of Heil's 2018
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
! damping_qe8: for two CSF-averaged on-diagonal matrix element
!              values, returns the value of the QE8 exponential
!              damping function
!######################################################################
  function damping_qe8(av1,av2) result(func)

    use constants
    use bitglobal
    use hparam
    
    implicit none

    ! Function result
    real(dp)             :: func

    ! CSF-averaged on-diagonal matrix elements
    real(dp), intent(in) :: av1,av2

    ! Everything else
    real(dp)             :: DEp3
    
    !
    !  p1 exp(-p2 DeltaE^p3)
    !
    DEp3=abs(av1-av2)**hpar(5)
    func=hpar(3)*exp(-hpar(4)*DEp3)
    
    return
    
  end function damping_qe8

!######################################################################
! damping_rc: RC DFT/MRCI SR off-diagonal damping.
!             p1=hpar(8), p2=hpar(9), n=hpar(10).
!######################################################################
  function damping_rc(av1,av2) result(func)

    use constants
    use bitglobal
    use hparam

    implicit none

    real(dp)             :: func
    real(dp), intent(in) :: av1,av2
    real(dp)             :: DEp3

    DEp3=abs(av1-av2)**hpar(10)
    func=hpar(8)*exp(-hpar(9)*DEp3)

    return

  end function damping_rc

!######################################################################
! damping_rc_lr: RC DFT/MRCI LR off-diagonal damping.
!                p1_LR=hpar(5), p2_LR=hpar(6), n_LR=hpar(7).
!
!                Same functional form as the SR damping but with its own
!                parameters. p1_LR is held at 1, so at dE = 0 the LR
!                contribution is the undamped ab initio matrix element and
!                p2_LR alone sets how quickly it falls away with the gap.
!######################################################################
  function damping_rc_lr(av1,av2) result(func)

    use constants
    use bitglobal
    use hparam

    implicit none

    real(dp)             :: func
    real(dp), intent(in) :: av1,av2
    real(dp)             :: DEp3

    DEp3=abs(av1-av2)**hpar(7)
    func=hpar(5)*exp(-hpar(6)*DEp3)

    return

  end function damping_rc_lr

!######################################################################
! damping_rc5p: SR off-diagonal damping for the reduced RC Hamiltonian.
!               p1=hpar(4), p2=hpar(5), n=hpar(6).
!
!               There is no LR counterpart: the long-range contribution
!               is undamped by construction.
!######################################################################
  function damping_rc5p(av1,av2) result(func)

    use constants
    use bitglobal
    use hparam

    implicit none

    real(dp)             :: func
    real(dp), intent(in) :: av1,av2
    real(dp)             :: DEp3

    DEp3=abs(av1-av2)**hpar(6)
    func=hpar(4)*exp(-hpar(5)*DEp3)

    return

  end function damping_rc5p

!######################################################################
! damping_rc7p: off-diagonal damping for the physically-motivated RC
!               Hamiltonian. One damping function, shared by both
!               channels; only the pre-factor differs.
!
!               p1 is passed in (hpar(4) for SR, hpar(5) for LR);
!               p2=hpar(6) and n=hpar(7) are common to both.
!######################################################################
  function damping_rc7p(av1,av2,p1) result(func)

    use constants
    use bitglobal
    use hparam

    implicit none

    real(dp)             :: func
    real(dp), intent(in) :: av1,av2,p1
    real(dp)             :: DEp3

    DEp3=abs(av1-av2)**hpar(7)
    func=p1*exp(-hpar(6)*DEp3)

    return

  end function damping_rc7p

!######################################################################
! damping_rc8p: off-diagonal damping for rc_dftmrci_8p. As rc_dftmrci_7p
!               -- one shared damping function, differing only in the
!               pre-factor -- but the parameter block is two longer
!               because the Coulomb scaling is split.
!
!               p1 is passed in (hpar(5) for SR, hpar(6) for LR);
!               p2=hpar(7) and n=hpar(8) are common to both.
!######################################################################
  function damping_rc8p(av1,av2,p1) result(func)

    use constants
    use bitglobal
    use hparam

    implicit none

    real(dp)             :: func
    real(dp), intent(in) :: av1,av2,p1
    real(dp)             :: DEp3

    DEp3=abs(av1-av2)**hpar(8)
    func=p1*exp(-hpar(7)*DEp3)

    return

  end function damping_rc8p

!######################################################################
! damping_r2022: for two CSF-averaged on-diagonal matrix element
!                values, returns the value of the R2022 DFT/MRCI
!                damping function
!######################################################################
  function damping_r2022(av1,av2) result(func)

    use constants
    use bitglobal
    use hparam
    
    implicit none

    ! Function result
    real(dp)             :: func

    ! CSF-averaged on-diagonal matrix elements
    real(dp), intent(in) :: av1,av2

    ! Everything else
    real(dp)             :: p1,p2,DE4
    
    !
    !  p1 exp(-p2 DeltaE^4)
    !
    p1=1.0d0-2.0d0*hpar(3)+hpar(5)
    p2=hpar(1)
    DE4=(av1-av2)**4
    func=p1*exp(-p2*DE4)
    
    return
    
  end function damping_r2022

!######################################################################
! damping_2026: for two CSF-averaged on-diagonal matrix element
!               values, returns the value of the 2026 DFT/MRCI
!               damping function.
!
!               p1 = 1 - 2*pJ + pF_hhee  (derived from the
!               degeneracy-preserving constraint; not a free parameter)
!######################################################################
  function damping_2026(av1,av2) result(func)

    use constants
    use bitglobal
    use hparam

    implicit none

    ! Function result
    real(dp)             :: func

    ! CSF-averaged on-diagonal matrix elements
    real(dp), intent(in) :: av1,av2

    ! Everything else
    real(dp)             :: p1,DEn

    !
    !  p1 exp(-p2 DeltaE^n),  p1 = 1 - 2*pJ + pF_hhee
    !
    p1=1.0d0-2.0d0*hpar(1)+hpar(3)
    DEn=abs(av1-av2)**hpar(5)
    func=p1*exp(-hpar(4)*DEn)

    return

  end function damping_2026

!######################################################################
! symvx: Exchange integral V_pqqp averaged over degenerate partners
!        of p and/or q (via degen_orbs). Enforces D∞h rotational
!        symmetry by returning the rotationally-averaged exchange:
!          - one partner: average over that partner
!          - both partners (not each other): average over p's partner
!            → gives K_avg = 0.5*(K_par + K_perp), equal for Au and B1u
!              components of Δu states
!          - partners of each other (e.g. πu_x/πu_y): bare Vx
!          - no partners: bare Vx
!######################################################################
  function symvx(p, q) result(val)

    use constants
    use bitglobal

    implicit none

    integer(is), intent(in) :: p, q
    real(dp)                :: val
    integer(is)             :: ip, jq

    ip = degen_orbs(p)
    jq = degen_orbs(q)

    if (ip == 0 .and. jq == 0) then
       val = Vx(p, q)
    else if (ip == q) then
       ! p and q are degenerate partners of each other: skip averaging
       val = Vx(p, q)
    else if (ip /= 0 .and. jq == 0) then
       val = 0.5d0 * (Vx(p, q) + Vx(ip, q))
    else if (ip == 0 .and. jq /= 0) then
       val = 0.5d0 * (Vx(p, q) + Vx(p, jq))
    else
       ! Both p and q have degenerate partners (not each other).
       ! Average over p's partner: gives K_avg = 0.5*(K_par+K_perp)
       ! for both Au and B1u components of Δu states → exact degeneracy.
       val = 0.5d0 * (Vx(p, q) + Vx(ip, q))
    end if

    return

  end function symvx

!######################################################################
! symvx_lr: LR exchange integral V_pqqp^LR averaged over degenerate
!           partners of p and/or q. Mirrors symvx but uses Vx_lr.
!######################################################################
  function symvx_lr(p, q) result(val)

    use constants
    use bitglobal

    implicit none

    integer(is), intent(in) :: p, q
    real(dp)                :: val
    integer(is)             :: ip, jq

    ip = degen_orbs(p)
    jq = degen_orbs(q)

    if (ip == 0 .and. jq == 0) then
       val = Vx_lr(p, q)
    else if (ip == q) then
       val = Vx_lr(p, q)
    else if (ip /= 0 .and. jq == 0) then
       val = 0.5d0 * (Vx_lr(p, q) + Vx_lr(ip, q))
    else if (ip == 0 .and. jq /= 0) then
       val = 0.5d0 * (Vx_lr(p, q) + Vx_lr(p, jq))
    else
       val = 0.5d0 * (Vx_lr(p, q) + Vx_lr(ip, q))
    end if

    return

  end function symvx_lr

!######################################################################

end module dftmrci
