!**********************************************************************
! Diagnostic accumulation of the SR and LR contributions to the
! off-diagonal Hamiltonian matrix elements, resolved by the
! configuration-averaged energy gap dE = |bav - kav|.
!
! This exists to answer three questions about the rc_dftmrci
! (ihamiltonian=18) LR/SR split:
!
!   1. do the SR and LR contributions occupy different ranges of dE?
!   2. do those distributions overlap?
!   3. is damping the LR contribution as a function of dE justified?
!
! The Epstein-Nesbet weight |H_IJ|^2/dE is the quantity that matters for
! (3): it is what the contribution is worth to the correlation energy at
! second order, and so what damping it away actually discards. It is
! accumulated separately for the SR and LR channels, together with the
! cross term, since H = H_SR + H_LR means
!
!   H^2 = H_SR^2 + 2 H_SR H_LR + H_LR^2
!
! and the cross term is not small in general.
!
! Inactive unless the environment variable GRACI_LRSR_DIAG is set, in
! which case its value is the file the histogram is written to.
!**********************************************************************
module lrsr_diag

  use constants

  implicit none

  ! Has the environment been queried yet?
  logical, private            :: lrsr_checked=.false.

  ! Is accumulation switched on?
  logical                     :: lrsr_active=.false.

  ! Output file
  character(len=255), private :: lrsr_file=''

  ! Histogram: linear bins in dE (Eh), plus one overflow bin
  integer(is), parameter      :: lrsr_nbin=500
  real(dp), parameter         :: lrsr_demax=5.0d0
  real(dp), parameter         :: lrsr_dbin=lrsr_demax/lrsr_nbin

  ! dE below which the Epstein-Nesbet weight is not accumulated: the
  ! 1/dE is genuinely singular there, and a handful of near-degenerate
  ! configuration pairs would otherwise swamp the histogram
  real(dp), parameter         :: lrsr_defloor=1.0d-4

  ! Accumulators, indexed (quantity,bin). Bin 0 is dE < defloor and
  ! bin lrsr_nbin+1 is dE > demax.
  integer(is), parameter      :: lrsr_nq=11
  real(dp), private           :: acc(lrsr_nq,0:lrsr_nbin+1)=0.0d0

  ! Quantity indices
  integer(is), parameter      :: iq_cnt   = 1  ! number of elements
  integer(is), parameter      :: iq_sr2   = 2  ! sum H_SR^2
  integer(is), parameter      :: iq_lr2   = 3  ! sum H_LR^2
  integer(is), parameter      :: iq_srlr  = 4  ! sum H_SR*H_LR
  integer(is), parameter      :: iq_ensr  = 5  ! sum H_SR^2/dE
  integer(is), parameter      :: iq_enlr  = 6  ! sum H_LR^2/dE
  integer(is), parameter      :: iq_encr  = 7  ! sum 2*H_SR*H_LR/dE
  integer(is), parameter      :: iq_endmp = 8  ! sum (damped H)^2/dE
  integer(is), parameter      :: iq_abssr = 9  ! sum |H_SR|
  integer(is), parameter      :: iq_abslr =10  ! sum |H_LR|
  integer(is), parameter      :: iq_ndmp  =11  ! sum of (damp,damp_lr) pairs

  ! Damping factors are functions of dE alone, so one representative
  ! value per bin suffices
  real(dp), private           :: dmp(2,0:lrsr_nbin+1)=0.0d0

  ! Per-element "LR-ness": r = |H_LR| / (|H_SR| + |H_LR|), in [0,1].
  ! This is the question the dE histogram cannot answer -- whether an
  ! individual matrix element is exclusively LR or exclusively SR (r
  ! bimodal at 0 and 1) or whether every element is a mixture (r
  ! unimodal). Only in the bimodal case can a per-element criterion
  ! separate the two channels.
  integer(is), parameter      :: lrsr_nfbin=100
  integer(is), parameter      :: lrsr_nfq=5
  real(dp), private           :: facc(lrsr_nfq,0:lrsr_nfbin)=0.0d0

  integer(is), parameter      :: if_cnt = 1  ! number of elements
  integer(is), parameter      :: if_en  = 2  ! sum H_full^2/dE  (EN weight)
  integer(is), parameter      :: if_h2  = 3  ! sum H_full^2
  integer(is), parameter      :: if_de  = 4  ! sum dE (for the mean gap)
  integer(is), parameter      :: if_sgn = 5  ! count with H_SR*H_LR > 0

contains

!######################################################################
! lrsr_init: reads GRACI_LRSR_DIAG on first use
!######################################################################
  subroutine lrsr_init

    implicit none

    integer(is) :: length,status

    if (lrsr_checked) return

    call get_environment_variable('GRACI_LRSR_DIAG',lrsr_file,length,status)

    lrsr_active=(status == 0 .and. length > 0)
    lrsr_checked=.true.

    return

  end subroutine lrsr_init

!######################################################################
! lrsr_accumulate: bins one batch of off-diagonal elements
!
!                  hij_full holds H_SR + H_LR and hij_lr holds H_LR, as
!                  they are held in hij_dftmrci_batch on entry, i.e.
!                  before any damping has been applied.
!######################################################################
  subroutine lrsr_accumulate(hij_full,hij_lr,n,bav,kav,damp,damp_lr)

    implicit none

    integer(is), intent(in) :: n
    real(dp), intent(in)    :: hij_full(n),hij_lr(n)
    real(dp), intent(in)    :: bav,kav,damp,damp_lr

    integer(is)             :: i,ibin,ifbin
    real(dp)                :: de,rde,hsr,hlr,hdmp,denom
    real(dp)                :: c(lrsr_nq)
    real(dp)                :: fc(lrsr_nfq,0:lrsr_nfbin)

    de=abs(bav-kav)

    !
    ! Locate the bin
    !
    if (de < lrsr_defloor) then
       ibin=0
       rde=0.0d0
    else if (de > lrsr_demax) then
       ibin=lrsr_nbin+1
       rde=1.0d0/de
    else
       ibin=min(int(de/lrsr_dbin)+1,lrsr_nbin)
       rde=1.0d0/de
    endif

    !
    ! Accumulate this batch into a local array first, so that the
    ! critical section is entered once per batch rather than once per
    ! matrix element
    !
    c=0.0d0
    fc=0.0d0
    do i=1,n
       hlr=hij_lr(i)
       hsr=hij_full(i)-hlr
       hdmp=damp*hsr+damp_lr*hlr

       ! per-element LR fraction
       denom=abs(hsr)+abs(hlr)
       if (denom > 0.0d0) then
          ifbin=min(int(abs(hlr)/denom*lrsr_nfbin),lrsr_nfbin)
          fc(if_cnt,ifbin)=fc(if_cnt,ifbin)+1.0d0
          fc(if_en ,ifbin)=fc(if_en ,ifbin)+hij_full(i)**2*rde
          fc(if_h2 ,ifbin)=fc(if_h2 ,ifbin)+hij_full(i)**2
          fc(if_de ,ifbin)=fc(if_de ,ifbin)+de
          if (hsr*hlr > 0.0d0) fc(if_sgn,ifbin)=fc(if_sgn,ifbin)+1.0d0
       endif

       c(iq_cnt)  =c(iq_cnt)  +1.0d0
       c(iq_sr2)  =c(iq_sr2)  +hsr*hsr
       c(iq_lr2)  =c(iq_lr2)  +hlr*hlr
       c(iq_srlr) =c(iq_srlr) +hsr*hlr
       c(iq_ensr) =c(iq_ensr) +hsr*hsr*rde
       c(iq_enlr) =c(iq_enlr) +hlr*hlr*rde
       c(iq_encr) =c(iq_encr) +2.0d0*hsr*hlr*rde
       c(iq_endmp)=c(iq_endmp)+hdmp*hdmp*rde
       c(iq_abssr)=c(iq_abssr)+abs(hsr)
       c(iq_abslr)=c(iq_abslr)+abs(hlr)
    enddo
    c(iq_ndmp)=1.0d0

    !$omp critical
    acc(:,ibin)=acc(:,ibin)+c
    dmp(1,ibin)=damp
    dmp(2,ibin)=damp_lr
    facc=facc+fc
    !$omp end critical

    return

  end subroutine lrsr_accumulate

!######################################################################
! lrsr_write: writes the histogram, then zeroes it
!######################################################################
  subroutine lrsr_write

    implicit none

    integer(is) :: unit,ibin
    real(dp)    :: de

    if (.not. lrsr_active) return
    if (sum(acc(iq_cnt,:)) == 0.0d0) return

    open(newunit=unit,file=trim(lrsr_file),form='formatted',status='unknown')

    write(unit,'(a)') '# rc_dftmrci SR/LR off-diagonal diagnostic'
    write(unit,'(a)') '# dE, H_SR and H_LR in atomic units'
    write(unit,'(a)') '# bin 0: dE below the 1/dE floor, EN weights not accumulated'
    write(unit,'(a,es12.5)') '# defloor = ',lrsr_defloor
    write(unit,'(a,es12.5)') '# demax   = ',lrsr_demax
    write(unit,'(a)') '#'
    write(unit,'(a4,1x,a14,12(1x,a18))') &
         '#bin','dE_lo','count','sum_HSR2','sum_HLR2','sum_HSRHLR',&
         'EN_SR','EN_LR','EN_cross','EN_damped','sum_absHSR',&
         'sum_absHLR','damp','damp_lr'

    do ibin=0,lrsr_nbin+1
       if (acc(iq_cnt,ibin) == 0.0d0) cycle
       if (ibin == 0) then
          de=0.0d0
       else if (ibin == lrsr_nbin+1) then
          de=lrsr_demax
       else
          de=(ibin-1)*lrsr_dbin
       endif
       write(unit,'(i4,1x,es14.5e3,12(1x,es18.8e3))') &
            ibin,de,&
            acc(iq_cnt,ibin),acc(iq_sr2,ibin),acc(iq_lr2,ibin),&
            acc(iq_srlr,ibin),acc(iq_ensr,ibin),acc(iq_enlr,ibin),&
            acc(iq_encr,ibin),acc(iq_endmp,ibin),acc(iq_abssr,ibin),&
            acc(iq_abslr,ibin),dmp(1,ibin),dmp(2,ibin)
    enddo

    close(unit)

    !
    ! Per-element LR-fraction distribution, written alongside as .frac
    !
    open(newunit=unit,file=trim(lrsr_file)//'.frac',form='formatted',&
         status='unknown')
    write(unit,'(a)') '# per-element LR fraction r = |H_LR|/(|H_SR|+|H_LR|)'
    write(unit,'(a)') '# bimodal at r=0,1 => elements are exclusively SR or LR'
    write(unit,'(a)') '# unimodal        => every element is a mixture'
    write(unit,'(a5,1x,a12,5(1x,a18))') &
         '#bin','r_lo','count','EN_weight','sum_H2','sum_dE','n_samesign'
    do ibin=0,lrsr_nfbin
       if (facc(if_cnt,ibin) == 0.0d0) cycle
       write(unit,'(i5,1x,es12.5,5(1x,es18.8e3))') &
            ibin,dble(ibin)/dble(lrsr_nfbin),&
            facc(if_cnt,ibin),facc(if_en,ibin),facc(if_h2,ibin),&
            facc(if_de,ibin),facc(if_sgn,ibin)
    enddo
    close(unit)

    ! The accumulators are deliberately not zeroed: bitci_finalise runs
    ! once per CI section, so leaving them to accumulate means the file
    ! left behind at the end of the run covers the whole calculation
    ! rather than only its last section.

    return

  end subroutine lrsr_write

end module lrsr_diag
