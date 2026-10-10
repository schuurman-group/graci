!**********************************************************************
! mrsf_interface: C-bound entry points of the MRSF-TDDFT library
!**********************************************************************
module mrsf_interface

  use iso_c_binding
  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_space
  use mrsf_integrals
  use mrsf_sigma
  use mrsf_extended
  use mrsf_xcgrid, only: time_xc, nxc_calls, nxc_vecs, time_probe, nprobe_calls
  use mrsf_davidson
  use mrsf_etensor
  use mrsf_density

  implicit none

contains

!######################################################################
! mrsf_set_extended: select the extended method (must precede
! mrsf_initialise, which lays out the response vector)
!######################################################################
  subroutine mrsf_set_extended(flag1) bind(c, name='mrsf_set_extended')

    logical(c_bool), intent(in) :: flag1

    extended = logical(flag1)

  end subroutine mrsf_set_extended

!######################################################################
! mrsf_set_frozen_core: number of doubly occupied MOs (the lowest by
! orbital energy) excluded from the response space (must precede
! mrsf_initialise)
!######################################################################
  subroutine mrsf_set_frozen_core(nfc1) bind(c, name='mrsf_set_frozen_core')

    integer(is), intent(in) :: nfc1

    nfc_req = nfc1

  end subroutine mrsf_set_frozen_core

!######################################################################
! mrsf_set_exchange: 0 auto (explicit tensor when it fits mem_budget),
! 1 explicit (ij|ab) tensor, 2 DF plane sweep (before mrsf_int_initialise)
!######################################################################
  subroutine mrsf_set_exchange(mode1) bind(c, name='mrsf_set_exchange')

    integer(is), intent(in) :: mode1

    exchange_mode = mode1

  end subroutine mrsf_set_exchange

!######################################################################
! mrsf_initialise
!######################################################################
  subroutine mrsf_initialise(nmo1, nel1, imult1, mosym1, occ1, moen1, fa1, fb1, &
       chf1, spc1, escf1, ipg1, label1, verbose1) bind(c, name='mrsf_initialise')

    integer(is), intent(in)            :: nmo1, nel1, imult1, ipg1
    integer(ib), intent(in)            :: mosym1(nmo1)
    real(dp), intent(in)               :: occ1(nmo1), moen1(nmo1)
    real(dp), intent(in)               :: fa1(nmo1,nmo1), fb1(nmo1,nmo1)
    real(dp), intent(in)               :: chf1, spc1(3), escf1
    character(kind=c_char), intent(in) :: label1(*)
    logical(c_bool), intent(in)        :: verbose1
    integer(is)                        :: a, b, i, j

    ! release the orbital data of a previous initialisation; the
    ! integrals are kept and reused if their key matches
    call free_orbital_data()

    nmo       = nmo1
    nel       = nel1
    imult_ref = imult1
    ipg       = ipg1
    chf       = chf1
    spc       = spc1
    escf      = escf1
    verbose   = logical(verbose1)
    call c2fstr(label1, label)

    if (imult_ref /= 3) call mrsf_error('the reference must be a triplet')

    allocate(occ(nmo), moen(nmo), mosym(nmo))
    occ   = occ1
    moen  = moen1
    mosym = int(mosym1, is)

    call setup_space()

    ! Fock matrix blocks
    allocate(FaHH(nocca,nocca), FbPP(nvirb,nvirb))
    do j = 1, nocca
       do i = 1, nocca
          FaHH(i,j) = fa1(Hmap(i),Hmap(j))
       enddo
    enddo
    do b = 1, nvirb
       do a = 1, nvirb
          FbPP(a,b) = fb1(Pmap(a),Pmap(b))
       enddo
    enddo

    ! the work buffer is laid out for the current column count
    if (allocated(Twork)) deallocate(Twork)
    nvmax = 0

    init_done = .true.

    if (verbose) then
       write(6,'(/,2x,a)') 'MRSF-TDDFT library initialised'
       write(6,'(2x,a,i0,a,i0,a,i0,a,i0)') 'nmo = ', nmo, ', nC = ', nC, &
            ', nV = ', nV, ', xdim = ', xdim
       if (extended) write(6,'(2x,a,i0,a,i0)') 'extended method: CV slots = ', ncv, &
            ', xdim_tot = ', xdim_tot
       write(6,'(2x,a,i0,a,i0)') 'SOMOs (MO indices): ', iO1, ', ', iO2
       write(6,'(2x,a,f8.4)') 'HF exchange fraction: ', chf
       if (nfc > 0) write(6,'(2x,a,i0,a,i0,a)') 'frozen core: ', nfc, &
            ' doubly occupied MOs excluded (', nocca_act, ' active holes)'
    endif

  end subroutine mrsf_initialise

!######################################################################
! mrsf_int_initialise
!######################################################################
  subroutine mrsf_int_initialise(method1, prec1, erifile1, vvstore1, membudget1) &
       bind(c, name='mrsf_int_initialise')

    character(kind=c_char), intent(in) :: method1(*), prec1(*), erifile1(*), vvstore1(*)
    real(dp), intent(in)               :: membudget1
    character(len=255)                 :: method, prec, erifile, vvstore

    call c2fstr(method1, method)
    call c2fstr(prec1, prec)
    call c2fstr(erifile1, erifile)
    call c2fstr(vvstore1, vvstore)

    if (trim(method) /= 'df') call mrsf_error('only density-fitted integrals are supported')
    mem_budget = membudget1 * 1.0e9_dp

    call load_df_file(erifile, prec, vvstore, mem_budget)

  end subroutine mrsf_int_initialise

!######################################################################
! mrsf_ext_initialise: extended method data (after mrsf_int_initialise
! and, if the kernel is used, after the grid initialisation)
!   fdft1(nmo,nmo): closed-shell KS matrix of G in the MO basis
!######################################################################
  subroutine mrsf_ext_initialise(nmo1, fdft1, ccp1, use_kernel1, fprime1, erifile1) &
       bind(c, name='mrsf_ext_initialise')

    integer(is), intent(in)            :: nmo1
    real(dp), intent(in)               :: fdft1(nmo1,nmo1)
    real(dp), intent(in)               :: ccp1
    logical(c_bool), intent(in)        :: use_kernel1
    integer(is), intent(in)            :: fprime1
    character(kind=c_char), intent(in) :: erifile1(*)
    character(len=255)                 :: erifile

    if (nmo1 /= nmo) call mrsf_error('mrsf_ext_initialise: inconsistent nmo')
    if (.not. extended) call mrsf_error('mrsf_ext_initialise: mrsf_set_extended not set')
    call c2fstr(erifile1, erifile)
    if (.not. allocated(Bcv)) call load_bcv_pass(erifile)
    call ext_initialise(fdft1, ccp1, logical(use_kernel1), fprime1)

  end subroutine mrsf_ext_initialise

!######################################################################
! mrsf_get_dims: xdim1 is the full vector length (xdim_tot)
!######################################################################
  subroutine mrsf_get_dims(nocca1, nvirb1, xdim1, naux1) bind(c, name='mrsf_get_dims')

    integer(is), intent(out) :: nocca1, nvirb1, xdim1, naux1

    nocca1 = nocca
    nvirb1 = nvirb
    xdim1  = xdim_tot
    naux1  = naux

  end subroutine mrsf_get_dims

!######################################################################
! mrsf_get_ext_dims: column count and number of CV slots
!######################################################################
  subroutine mrsf_get_ext_dims(ncol1, ncv1) bind(c, name='mrsf_get_ext_dims')

    integer(is), intent(out) :: ncol1, ncv1

    ncol1 = ncol
    ncv1  = ncv

  end subroutine mrsf_get_ext_dims

!######################################################################
! mrsf_diag: Davidson diagonalisation for one (irrep, mult) block
!######################################################################
  subroutine mrsf_diag(irrep1, imult1, nroots1, nextra1, maxvec1, maxiter1, tol1, &
       precond1, inner1, xdim1, ener1, xvec1, niter1, iconv1) bind(c, name='mrsf_diag')

    integer(is), intent(in)    :: irrep1, imult1, nextra1, maxvec1, maxiter1, xdim1
    integer(is), intent(in)    :: precond1, inner1
    integer(is), intent(inout) :: nroots1
    real(dp), intent(in)       :: tol1
    real(dp), intent(out)      :: ener1(*), xvec1(*)
    integer(is), intent(out)   :: niter1, iconv1
    real(dp)                   :: t0

    if (xdim1 /= xdim_tot) call mrsf_error('mrsf_diag: inconsistent xdim')
    if (imult1 /= 1 .and. imult1 /= 3) call mrsf_error('mrsf_diag: mult must be 1 or 3')
    if (extended .and. .not. ext_ready) call mrsf_error('mrsf_diag: mrsf_ext_initialise not called')
    if (precond1 /= 1 .and. precond1 /= 2) call mrsf_error('mrsf_diag: precond must be 1 or 2')

    t0 = wall_time()
    call davidson_solve(imult1, irrep1, nroots1, nextra1, maxvec1, maxiter1, tol1, &
         precond1, max(inner1, 1_is), ener1, xvec1, niter1, iconv1)
    if (verbose) write(6,'(2x,a,f10.2,a,/)') 'Davidson time: ', wall_time()-t0, ' s'

  end subroutine mrsf_diag

!######################################################################
! mrsf_sigma: A x for a batch of compressed vectors (validation)
!######################################################################
  subroutine mrsf_sigma_c(imult1, irrep1, nvec1, xdim1, x1, ax1) bind(c, name='mrsf_sigma')

    integer(is), intent(in) :: imult1, irrep1, nvec1, xdim1
    real(dp), intent(in)    :: x1(xdim1,nvec1)
    real(dp), intent(out)   :: ax1(xdim1,nvec1)

    if (xdim1 /= xdim_tot) call mrsf_error('mrsf_sigma: inconsistent xdim')
    if (extended .and. .not. ext_ready) call mrsf_error('mrsf_sigma: mrsf_ext_initialise not called')
    call sigma_batch(nvec1, imult1, irrep1, x1, ax1)

  end subroutine mrsf_sigma_c

!######################################################################
! mrsf_diagonal
!######################################################################
  subroutine mrsf_diagonal(imult1, xdim1, d1) bind(c, name='mrsf_diagonal')

    integer(is), intent(in) :: imult1, xdim1
    real(dp), intent(out)   :: d1(xdim1)

    if (xdim1 /= xdim_tot) call mrsf_error('mrsf_diagonal: inconsistent xdim')
    call diagonal(imult1, d1)

  end subroutine mrsf_diagonal

!######################################################################
! mrsf_density: state 1-RDMs (spin-summed, MO basis); ncol1 = nocca
! (standard) or nocca + nC (extended); arrays passed by address
!######################################################################
  subroutine mrsf_density_c(nmo1, occ1, imult1, ncol1, xdim1, nroots1, xvec1, dmat1) &
       bind(c, name='mrsf_density')

    integer(is), intent(in) :: nmo1, imult1, ncol1, xdim1, nroots1
    real(dp), intent(in)    :: occ1(nmo1), xvec1(xdim1,nroots1)
    real(dp), intent(out)   :: dmat1(nmo1,nmo1,nroots1)
    integer(is), allocatable :: ipairs(:,:)
    integer(is)              :: k

    allocate(ipairs(2,nroots1))
    do k = 1, nroots1
       ipairs(:,k) = k
    enddo
    call tdm_pairs(nmo1, occ1, imult1, ncol1, xdim1, nroots1, nroots1, nroots1, ipairs, &
         xvec1, xvec1, dmat1)
    deallocate(ipairs)

  end subroutine mrsf_density_c

!######################################################################
! mrsf_tdm: 1-TDMs <bra|E_pq|ket> for pairs of states
!######################################################################
  subroutine mrsf_tdm(nmo1, occ1, imult1, ncol1, xdim1, npairs1, nb1, nk1, ipairs1, xb1, xk1, &
       rho1) bind(c, name='mrsf_tdm')

    integer(is), intent(in) :: nmo1, imult1, ncol1, xdim1, npairs1, nb1, nk1
    real(dp), intent(in)    :: occ1(nmo1)
    integer(is), intent(in) :: ipairs1(2,npairs1)
    real(dp), intent(in)    :: xb1(xdim1,nb1), xk1(xdim1,nk1)
    real(dp), intent(out)   :: rho1(nmo1,nmo1,npairs1)

    call tdm_pairs(nmo1, occ1, imult1, ncol1, xdim1, npairs1, nb1, nk1, ipairs1, xb1, xk1, rho1)

  end subroutine mrsf_tdm

!######################################################################
! mrsf_finalise
!######################################################################
  subroutine mrsf_finalise() bind(c, name='mrsf_finalise')

    call mrsf_finalise_f()

  end subroutine mrsf_finalise

!######################################################################
! mrsf_report_timings: prints the accumulated timings of the current
! section (if verbose) and resets the counters, leaving the loaded
! integrals and orbital data in place
!######################################################################
  subroutine mrsf_report_timings() bind(c, name='mrsf_report_timings')

    call report_timings()

  end subroutine mrsf_report_timings

  subroutine free_orbital_data()

    if (allocated(occ)) deallocate(occ)
    if (allocated(moen)) deallocate(moen)
    if (allocated(mosym)) deallocate(mosym)
    if (allocated(Hmap)) deallocate(Hmap)
    if (allocated(Pmap)) deallocate(Pmap)
    if (allocated(Hinv)) deallocate(Hinv)
    if (allocated(Pinv)) deallocate(Pinv)
    if (allocated(slot_irrep)) deallocate(slot_irrep)
    if (allocated(FaHH)) deallocate(FaHH)
    if (allocated(FbPP)) deallocate(FbPP)
    if (allocated(frozen_hole)) deallocate(frozen_hole)
    if (allocated(hact)) deallocate(hact)
    if (allocated(cact)) deallocate(cact)
    nfc = 0
    call ext_free()
    init_done = .false.

  end subroutine free_orbital_data

  subroutine mrsf_finalise_f()

    call report_timings()
    call free_ints()
    call free_orbital_data()
    extended = .false.
    nfc_req = 0

  end subroutine mrsf_finalise_f

  subroutine report_timings()

    real(dp) :: gflop

    if (verbose .and. (nsigma_calls > 0 .or. nxc_calls > 0)) then
       write(6,'(/,2x,a)') 'MRSF timings'
       if (nsigma_calls > 0) then
          write(6,'(2x,a,f10.2,a)') 'integral ingestion : ', time_load, ' s'
          write(6,'(2x,a,f10.2,a,i0,a,i0,a)') 'sigma vectors      : ', time_sigma, &
               ' s (', nsigma_calls, ' calls, ', nsigma_vecs, ' vectors)'
          if (ninner_vecs > 0) write(6,'(2x,a,f10.2,a,i0,a)') 'inner corrections  : ', &
               time_inner, ' s (', ninner_vecs, ' kernel-free vectors)'
          write(6,'(2x,a,f10.2,a)') '  exchange term    : ', time_exch, ' s'
          if (time_exch > 0.0_dp) then
             ! sweep: 2 naux nvirb^2 ncol + step 2: 2 naux nvirb nocca ncol per vector
             gflop = 2.0_dp * real(naux,dp) * real(nvirb,dp) * real(ncol,dp) &
                  * real(nvirb + nocca,dp) * real(nsigma_vecs,dp) / 1.0e9_dp
             write(6,'(2x,a,f10.1,a,f10.2,a)') '  exchange kernel  : ', gflop, &
                  ' GFLOP, ', gflop / time_exch, ' GFLOP/s'
          endif
          if (netens_vecs > 0) then
             write(6,'(2x,a,f10.2,a,f10.2,a,i0,a,f10.1,a,f10.2,a)') '  exchange tensor  : build ', &
                  time_ebuild, ' s, sigma ', time_etens, ' s (', netens_vecs, ' vectors, ', &
                  flop_etens/1.0e9_dp, ' GFLOP, ', flop_etens/1.0e9_dp/max(time_etens,1.0e-12_dp), ' GFLOP/s)'
          endif
          if (extended) then
             write(6,'(2x,a,f10.2,a)') '  extended terms   : ', time_ext, ' s (Coulomb pass, couplings)'
             if (nxcs_vecs > 0) then
                ! psiL (2 nao nC_act ncomp) + three projection dgemms (6 nC_act nao) per point
                gflop = (2.0_dp*real(nao_g_report(),dp)*real(nC_act,dp)*real(ncomp_report(),dp) &
                     + 6.0_dp*real(nC_act,dp)*real(nao_g_report(),dp)) * real(ngrid_report(),dp) &
                     * real(nxcs_vecs,dp) / 1.0e9_dp
                write(6,'(2x,a,f10.2,a,i0,a,f10.1,a,f10.2,a)') '  CV kernel (grid) : ', time_xcs, &
                     ' s (', nxcs_vecs, ' vectors, ', gflop, ' GFLOP dgemm, ', gflop/time_xcs, ' GFLOP/s)'
             endif
          endif
       endif
       if (nxc_calls > 0 .and. .not. extended) write(6,'(2x,a,f10.2,a,i0,a,i0,a)') &
            'xc kernel (grid)   : ', time_xc, ' s (', nxc_calls, ' calls, ', nxc_vecs, ' vectors)'
       if (nprobe_calls > 0) write(6,'(2x,a,f10.2,a,i0,a)') 'xc probe (grid)    : ', time_probe, &
            ' s (', nprobe_calls, ' calls)'
    endif

    time_load = 0.0_dp; time_sigma = 0.0_dp; time_exch = 0.0_dp
    time_ext = 0.0_dp; time_xcs = 0.0_dp; nxcs_vecs = 0
    nsigma_calls = 0; nsigma_vecs = 0
    time_ebuild = 0.0_dp; time_etens = 0.0_dp; flop_etens = 0.0_dp; netens_vecs = 0
    time_inner = 0.0_dp; ninner_vecs = 0
    time_xc = 0.0_dp; nxc_calls = 0; nxc_vecs = 0; time_probe = 0.0_dp; nprobe_calls = 0

  end subroutine report_timings

  ! grid dimensions for the kernel flop count (kept out of the use list
  ! of mrsf_xcgrid to avoid name clashes)
  function nao_g_report() result(n)
    use mrsf_xcgrid, only: nao_g
    integer(is) :: n
    n = nao_g
  end function nao_g_report

  function ncomp_report() result(n)
    use mrsf_xcgrid, only: ncomp
    integer(is) :: n
    n = ncomp
  end function ncomp_report

  function ngrid_report() result(n)
    use mrsf_xcgrid, only: ngrid
    integer(is) :: n
    n = ngrid
  end function ngrid_report

end module mrsf_interface
