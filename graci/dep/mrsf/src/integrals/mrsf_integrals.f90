!**********************************************************************
! mrsf_integrals: ingestion of the density-fitted MO integrals written
! by graci.core.ao2mo (packed (naux, nmo(nmo+1)/2) tensor, Fortran
! unformatted file) into the run-time layouts:
!   Bvv(nvirb,nvirb,nplane): paired symmetric planes, plane k holds
!       B^{2k-1}_{ab} (a>=b, lower) and B^{2k}_{ab} (a<b, upper)
!       [or full planes Bvv(nvirb,nvirb,naux) when vv_full]
!   Dall(naux,nmo): diagonals B^Q_pp
!   Boo(nQ,nocca,nocca,nblk): hole-hole block, Q-blocked
!   Bco(naux,nC,2), Bvo(naux,nV,2): SOMO slices
!   Bcv(nvirb,nC,naux): particle-core block (extended method only)
! plus the pairing-strength blocks G, M, N and the diagonal.
!**********************************************************************
module mrsf_integrals

  use mrsf_constants
  use mrsf_global
  use mrsf_etensor
  use mrsf_io
  use mrsf_space

  implicit none

  integer(is), allocatable, private :: pair_p(:), pair_q(:)

contains

!######################################################################
! load_df_file
!######################################################################
  subroutine load_df_file(fname, prec, vvstore, budget_bytes)

    character(len=*), intent(in) :: fname, prec, vvstore
    real(dp), intent(in)         :: budget_bytes

    integer(is)           :: unit, dims(2), nrec, cpr, irec, c0, c1, ncol1
    integer(is)           :: n_ij, nv_est, p, q, ij
    logical               :: exists
    real(dp), allocatable :: buf(:,:)
    real(sp), allocatable :: buf_sp(:,:)
    real(dp)              :: t0, gb, rq

    t0 = wall_time()

    if (.not. init_done) call mrsf_error('mrsf_initialise must be called first')

    ! exchange term: explicit tensor or DF sweep
    select case (exchange_mode)
    case (1)
       use_etensor = .true.
    case (2)
       use_etensor = .false.
    case default
       use_etensor = (etensor_bytes() <= budget_bytes)
    end select
    if (verbose) then
       if (use_etensor) then
          write(6,'(/,2x,a,f8.3,a)') 'exchange term: explicit (ij|ab) tensor (', &
               etensor_bytes()/1.0e9_dp, ' GB)'
       else
          write(6,'(/,2x,a)') 'exchange term: DF plane sweep'
       endif
    endif

    ! reuse the loaded integrals if the key (file, precision, storage,
    ! orbital classes) is unchanged: only the Fock-dependent diagonal
    ! has to be rebuilt (and the CV block loaded if newly needed)
    if (ints_loaded .and. ints_match(fname, prec, vvstore)) then
       if (verbose) write(6,'(/,2x,a)') 'MRSF integrals already loaded: reusing them'
       call build_diag0()
       if (extended .and. .not. allocated(Bcv)) call load_bcv_pass(fname)
       if (use_etensor) then
          if (.not. etensor_ready .or. etensor_nfc /= nfc .or. &
               (extended .and. .not. allocated(Eo2))) call etensor_build()
       else if (etensor_ready) then
          call etensor_free()
       endif
       return
    endif

    inquire(file=trim(fname), exist=exists)
    if (.not. exists) call mrsf_error('integral file not found: '//trim(fname))

    store_sp = (trim(prec) == 'single')
    vv_full  = (trim(vvstore) == 'full')
    if (extended .and. store_sp) &
         call mrsf_error('the extended method requires double-precision integrals')

    ! header: dimensions, number of records, columns per record
    call freeunit(unit)
    open(unit, file=trim(fname), form='unformatted', status='old')
    read(unit) dims(1)
    read(unit) dims(2)
    naux = dims(1)
    n_ij = dims(2)
    if (n_ij /= nmo*(nmo+1)/2) &
         call mrsf_error('integral file pair dimension does not match nmo')
    read(unit) nrec
    read(unit) cpr

    ! block structure (the T buffer holds ncol columns per vector)
    nplane = (naux + 1) / 2
    nv_est = 16
    rq = budget_bytes / (8.0_dp * real(nvirb,dp) * real(ncol,dp) * real(nv_est,dp))
    if (rq > real(naux + 2, dp)) then
       nQ = naux + mod(naux, 2_is)
    else
       nQ = max(2_is, int(rq, is))
    endif
    nQ = (nQ / 2) * 2
    npblk = nQ / 2
    nblk  = (naux + nQ - 1) / nQ

    ! pair tables (1-based, p >= q, ij = p(p-1)/2 + q)
    if (allocated(pair_p)) deallocate(pair_p, pair_q)
    allocate(pair_p(n_ij), pair_q(n_ij))
    ij = 0
    do p = 1, nmo
       do q = 1, p
          ij = ij + 1
          pair_p(ij) = p
          pair_q(ij) = q
       enddo
    enddo

    ! allocate the integral arrays
    call free_ints()
    if (store_sp) then
       if (vv_full) then
          allocate(Bvv_sp(nvirb,nvirb,naux), source=0.0_sp)
       else
          allocate(Bvv_sp(nvirb,nvirb,nplane), source=0.0_sp)
       endif
       allocate(Boo_sp(nQ,nocca,nocca,nblk), source=0.0_sp)
    else
       if (vv_full) then
          allocate(Bvv(nvirb,nvirb,naux), source=0.0_dp)
       else
          allocate(Bvv(nvirb,nvirb,nplane), source=0.0_dp)
       endif
       allocate(Boo(nQ,nocca,nocca,nblk), source=0.0_dp)
    endif
    allocate(Dall(naux,nmo), source=0.0_dp)
    allocate(Bco(naux,max(nC,1_is),2), source=0.0_dp)
    allocate(Bvo(naux,max(nV,1_is),2), source=0.0_dp)
    if (extended) allocate(Bcv(nvirb,max(nC,1_is),naux), source=0.0_dp)

    if (verbose) then
       gb = real(nvirb,dp)**2 * real(merge(naux, nplane, vv_full),dp) &
            * merge(4.0_dp, 8.0_dp, store_sp) / 1.0e9_dp
       write(6,'(/,2x,a)') 'MRSF integral ingestion'
       write(6,'(2x,a,i0,a,i0,a,i0)') 'naux = ', naux, ', nplane = ', nplane, &
            ', Q-block size = ', nQ
       write(6,'(2x,a,f8.3,a)') 'vir-vir planes: ', gb, ' GB'
       if (extended) write(6,'(2x,a,f8.3,a)') 'particle-core block (extended): ', &
            8.0_dp*real(nvirb,dp)*real(nC,dp)*real(naux,dp)/1.0e9_dp, ' GB'
    endif

    ! stream the records
    allocate(buf(naux,cpr))
    if (store_sp .or. trim(prec) == 'single') allocate(buf_sp(naux,cpr))
    do irec = 1, nrec
       c0    = (irec-1)*cpr + 1
       c1    = min(irec*cpr, n_ij)
       ncol1 = c1 - c0 + 1
       if (trim(prec) == 'single') then
          read(unit) buf_sp(1:naux, 1:ncol1)
          buf(1:naux,1:ncol1) = real(buf_sp(1:naux,1:ncol1), dp)
       else
          read(unit) buf(1:naux, 1:ncol1)
       endif
       call scatter_record(c0, ncol1, buf)
    enddo
    close(unit)
    deallocate(buf)
    if (allocated(buf_sp)) deallocate(buf_sp)

    call build_pairing_blocks()
    call build_diag0()

    ints_loaded = .true.
    loaded_file = trim(fname)
    loaded_prec = trim(prec)
    loaded_vv   = trim(vvstore)
    loaded_nmo  = nmo
    if (allocated(loaded_occ)) deallocate(loaded_occ)
    allocate(loaded_occ(nmo))
    loaded_occ  = occ
    if (use_etensor) call etensor_build()
    time_load = wall_time() - t0
    if (verbose) write(6,'(2x,a,f10.2,a)') 'integral ingestion time: ', time_load, ' s'

  end subroutine load_df_file

!######################################################################
! ints_match: .true. if the loaded integrals correspond to the same
! file, precision, storage and orbital classification
!######################################################################
  function ints_match(fname, prec, vvstore) result(same)

    character(len=*), intent(in) :: fname, prec, vvstore
    logical                      :: same

    same = .false.
    if (.not. ints_loaded) return
    if (trim(fname) /= trim(loaded_file)) return
    if (trim(prec) /= trim(loaded_prec)) return
    if (trim(vvstore) /= trim(loaded_vv)) return
    if (nmo /= loaded_nmo) return
    if (.not. allocated(loaded_occ)) return
    if (maxval(abs(occ - loaded_occ)) > 1.0e-6_dp) return
    ! the T buffer / Q-blocking was sized for the loaded column count
    if (extended .and. nQ < naux .and. ncol > nocca) then
       if (8.0_dp*real(nvirb,dp)*real(nQ,dp)*real(ncol,dp)*16.0_dp > 1.5_dp*mem_budget) return
    endif
    same = .true.

  end function ints_match

!######################################################################
! scatter_record: scatter one record (all Q, columns c0..c0+ncol1-1 of
! the packed tensor) into the run-time layouts
!######################################################################
  subroutine scatter_record(c0, ncol1, buf)

    integer(is), intent(in) :: c0, ncol1
    real(dp), intent(in)    :: buf(naux, ncol1)

    integer(is), parameter  :: ktile = 32
    integer(is) :: kt, k1, jc, ij, p, q, a, b, i, j, k, iq, blk, Ql, x, y, tmp
    integer(is) :: nk

    nk = merge(naux, nplane, vv_full)

    ! vir-vir planes, tiled over planes
    do kt = 1, nk, ktile
       k1 = min(kt + ktile - 1, nk)
       !$omp parallel do private(jc, ij, p, q, a, b, k, tmp) schedule(static, 16)
       do jc = 1, ncol1
          ij = c0 + jc - 1
          p  = pair_p(ij)
          q  = pair_q(ij)
          a  = Pinv(p)
          b  = Pinv(q)
          if (a == 0 .or. b == 0) cycle
          if (a < b) then
             tmp = a; a = b; b = tmp
          endif
          if (vv_full) then
             if (store_sp) then
                do k = kt, k1
                   Bvv_sp(a,b,k) = real(buf(k,jc), sp)
                   Bvv_sp(b,a,k) = real(buf(k,jc), sp)
                enddo
             else
                do k = kt, k1
                   Bvv(a,b,k) = buf(k,jc)
                   Bvv(b,a,k) = buf(k,jc)
                enddo
             endif
          else
             if (store_sp) then
                do k = kt, k1
                   Bvv_sp(a,b,k) = real(buf(2*k-1,jc), sp)
                   if (2*k <= naux .and. a /= b) Bvv_sp(b,a,k) = real(buf(2*k,jc), sp)
                enddo
             else
                do k = kt, k1
                   Bvv(a,b,k) = buf(2*k-1,jc)
                   if (2*k <= naux .and. a /= b) Bvv(b,a,k) = buf(2*k,jc)
                enddo
             endif
          endif
       enddo
       !$omp end parallel do
    enddo

    ! diagonals, hole-hole block, SOMO slices, particle-core block
    !$omp parallel do private(jc, ij, p, q, i, j, iq, blk, Ql, x, y, a, b) schedule(static, 16)
    do jc = 1, ncol1
       ij = c0 + jc - 1
       p  = pair_p(ij)
       q  = pair_q(ij)
       if (p == q) Dall(:,p) = buf(:,jc)
       i = Hinv(p)
       j = Hinv(q)
       if (i > 0 .and. j > 0) then
          do iq = 1, naux
             blk = (iq-1)/nQ + 1
             Ql  = iq - (blk-1)*nQ
             if (store_sp) then
                Boo_sp(Ql,i,j,blk) = real(buf(iq,jc), sp)
                Boo_sp(Ql,j,i,blk) = real(buf(iq,jc), sp)
             else
                Boo(Ql,i,j,blk) = buf(iq,jc)
                Boo(Ql,j,i,blk) = buf(iq,jc)
             endif
          enddo
       endif
       x = somo_index(p)
       y = somo_index(q)
       ! core-SOMO
       if (x > 0 .and. j > 0 .and. j <= nC) Bco(:,j,x) = buf(:,jc)
       if (y > 0 .and. i > 0 .and. i <= nC) Bco(:,i,y) = buf(:,jc)
       ! virtual-SOMO
       if (x > 0 .and. Pinv(q) > 2) Bvo(:,Pinv(q)-2,x) = buf(:,jc)
       if (y > 0 .and. Pinv(p) > 2) Bvo(:,Pinv(p)-2,y) = buf(:,jc)
       ! particle-core (extended)
       if (extended) then
          a = Pinv(p)
          if (a > 0 .and. j > 0 .and. j <= nC) Bcv(a,j,:) = buf(:,jc)
          b = Pinv(q)
          if (b > 0 .and. i > 0 .and. i <= nC) Bcv(b,i,:) = buf(:,jc)
       endif
    enddo
    !$omp end parallel do

  end subroutine scatter_record

!######################################################################
! load_bcv_pass: second streaming pass over the file collecting only
! the particle-core block (integrals already resident from a
! non-extended section)
!######################################################################
  subroutine load_bcv_pass(fname)

    character(len=*), intent(in) :: fname

    integer(is)           :: unit, dims(2), nrec, cpr, irec, c0, c1, ncol1
    integer(is)           :: n_ij, p, q, ij, jc, a, b, i, j
    integer(is), allocatable :: ptab(:), qtab(:)
    real(dp), allocatable :: buf(:,:)
    logical               :: exists
    real(dp)              :: t0

    t0 = wall_time()
    if (.not. ints_loaded) call mrsf_error('load_bcv_pass: integrals not loaded')
    if (store_sp) call mrsf_error('the extended method requires double-precision integrals')
    inquire(file=trim(fname), exist=exists)
    if (.not. exists) call mrsf_error('integral file not found: '//trim(fname))

    if (allocated(Bcv)) deallocate(Bcv)
    allocate(Bcv(nvirb,max(nC,1_is),naux), source=0.0_dp)

    n_ij = nmo*(nmo+1)/2
    allocate(ptab(n_ij), qtab(n_ij))
    ij = 0
    do p = 1, nmo
       do q = 1, p
          ij = ij + 1
          ptab(ij) = p
          qtab(ij) = q
       enddo
    enddo

    call freeunit(unit)
    open(unit, file=trim(fname), form='unformatted', status='old')
    read(unit) dims(1)
    read(unit) dims(2)
    if (dims(1) /= naux .or. dims(2) /= n_ij) then
       close(unit)
       call mrsf_error('load_bcv_pass: integral file does not match the loaded integrals')
    endif
    read(unit) nrec
    read(unit) cpr
    allocate(buf(naux,cpr))
    do irec = 1, nrec
       c0    = (irec-1)*cpr + 1
       c1    = min(irec*cpr, n_ij)
       ncol1 = c1 - c0 + 1
       read(unit) buf(1:naux, 1:ncol1)
       !$omp parallel do private(ij, p, q, a, b, i, j)
       do jc = 1, ncol1
          ij = c0 + jc - 1
          p = ptab(ij); q = qtab(ij)
          a = Pinv(p); j = Hinv(q)
          if (a > 0 .and. j > 0 .and. j <= nC) Bcv(a,j,:) = buf(:,jc)
          b = Pinv(q); i = Hinv(p)
          if (b > 0 .and. i > 0 .and. i <= nC) Bcv(b,i,:) = buf(:,jc)
       enddo
       !$omp end parallel do
    enddo
    close(unit)
    deallocate(buf, ptab, qtab)
    if (verbose) write(6,'(2x,a,f10.2,a)') 'particle-core block loaded in ', wall_time() - t0, ' s'

  end subroutine load_bcv_pass

!######################################################################
! somo_index: 1 for O1, 2 for O2, 0 otherwise
!######################################################################
  pure function somo_index(p) result(x)

    integer(is), intent(in) :: p
    integer(is)             :: x

    x = 0
    if (p == iO1) x = 1
    if (p == iO2) x = 2

  end function somo_index

!######################################################################
! build_pairing_blocks: G, M, N, H and K12 from the SOMO slices
!######################################################################
  subroutine build_pairing_blocks()

    integer(is) :: x, y, blk, nQl

    if (allocated(Gp)) deallocate(Gp, Mp, Np, Hp)
    allocate(Gp(max(nC,1_is),max(nC,1_is),2,2), source=0.0_dp)
    allocate(Mp(max(nV,1_is),max(nV,1_is),2,2), source=0.0_dp)
    allocate(Np(max(nC,1_is),max(nV,1_is),2,2), source=0.0_dp)
    allocate(Hp(max(nC,1_is),max(nV,1_is)), source=0.0_dp)

    do x = 1, 2
       do y = 1, 2
          if (nC > 0) call dgemm('T','N', nC, nC, naux, 1.0_dp, Bco(1,1,x), naux, &
               Bco(1,1,y), naux, 0.0_dp, Gp(1,1,x,y), nC)
          if (nV > 0) call dgemm('T','N', nV, nV, naux, 1.0_dp, Bvo(1,1,x), naux, &
               Bvo(1,1,y), naux, 0.0_dp, Mp(1,1,x,y), nV)
          if (nC > 0 .and. nV > 0) call dgemm('T','N', nC, nV, naux, 1.0_dp, &
               Bco(1,1,x), naux, Bvo(1,1,y), naux, 0.0_dp, Np(1,1,x,y), nC)
       enddo
    enddo
    Hp = Np(:,:,1,2) - Np(:,:,2,1)

    ! K12 = (O1 O2|O1 O2) from the hole-hole block
    K12 = 0.0_dp
    do blk = 1, nblk
       nQl = min(nQ, naux - (blk-1)*nQ)
       if (store_sp) then
          K12 = K12 + sum(real(Boo_sp(1:nQl,nC+1,nC+2,blk),dp)**2)
       else
          K12 = K12 + sum(Boo(1:nQl,nC+1,nC+2,blk)**2)
       endif
    enddo

  end subroutine build_pairing_blocks

!######################################################################
! build_diag0: multiplicity-independent diagonal
!   diag0(a,i) = Fb(a,a) - Fa(i,i) - c (ii|aa)
!######################################################################
  subroutine build_diag0()

    real(dp), allocatable :: Dpart(:,:), Dhole(:,:), J(:,:)
    integer(is)           :: a, i

    allocate(Dpart(naux,nvirb), Dhole(naux,nocca), J(nvirb,nocca))
    do a = 1, nvirb
       Dpart(:,a) = Dall(:,Pmap(a))
    enddo
    do i = 1, nocca
       Dhole(:,i) = Dall(:,Hmap(i))
    enddo
    call dgemm('T','N', nvirb, nocca, naux, 1.0_dp, Dpart, naux, Dhole, naux, 0.0_dp, J, nvirb)

    if (allocated(diag0)) deallocate(diag0)
    allocate(diag0(nvirb,nocca))
    do i = 1, nocca
       do a = 1, nvirb
          diag0(a,i) = FbPP(a,a) - FaHH(i,i) - chf * J(a,i)
       enddo
    enddo

    deallocate(Dpart, Dhole, J)

  end subroutine build_diag0

!######################################################################
! free_ints
!######################################################################
  subroutine free_ints()

    if (allocated(Bvv)) deallocate(Bvv)
    if (allocated(Bvv_sp)) deallocate(Bvv_sp)
    if (allocated(Dall)) deallocate(Dall)
    if (allocated(Boo)) deallocate(Boo)
    if (allocated(Boo_sp)) deallocate(Boo_sp)
    if (allocated(Bco)) deallocate(Bco)
    if (allocated(Bvo)) deallocate(Bvo)
    if (allocated(Bcv)) deallocate(Bcv)
    if (allocated(Gp)) deallocate(Gp, Mp, Np, Hp)
    if (allocated(diag0)) deallocate(diag0)
    if (allocated(Twork)) deallocate(Twork)
    if (allocated(plane_scr)) deallocate(plane_scr)
    if (allocated(Boo_scr)) deallocate(Boo_scr)
    if (allocated(loaded_occ)) deallocate(loaded_occ)
    call etensor_free()
    loaded_file = ''
    nvmax = 0
    ints_loaded = .false.

  end subroutine free_ints

end module mrsf_integrals
