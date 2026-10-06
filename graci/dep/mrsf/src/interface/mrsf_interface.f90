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
  use mrsf_davidson
  use mrsf_density

  implicit none

contains

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

    init_done = .true.

    if (verbose) then
       write(6,'(/,2x,a)') 'MRSF-TDDFT library initialised'
       write(6,'(2x,a,i0,a,i0,a,i0,a,i0)') 'nmo = ', nmo, ', nC = ', nC, &
            ', nV = ', nV, ', xdim = ', xdim
       write(6,'(2x,a,i0,a,i0)') 'SOMOs (MO indices): ', iO1, ', ', iO2
       write(6,'(2x,a,f8.4)') 'HF exchange fraction: ', chf
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
! mrsf_get_dims
!######################################################################
  subroutine mrsf_get_dims(nocca1, nvirb1, xdim1, naux1) bind(c, name='mrsf_get_dims')

    integer(is), intent(out) :: nocca1, nvirb1, xdim1, naux1

    nocca1 = nocca
    nvirb1 = nvirb
    xdim1  = xdim
    naux1  = naux

  end subroutine mrsf_get_dims

!######################################################################
! mrsf_diag: Davidson diagonalisation for one (irrep, mult) block
!######################################################################
  subroutine mrsf_diag(irrep1, imult1, nroots1, nextra1, maxvec1, maxiter1, tol1, &
       xdim1, ener1, xvec1, niter1, iconv1) bind(c, name='mrsf_diag')

    integer(is), intent(in)    :: irrep1, imult1, nextra1, maxvec1, maxiter1, xdim1
    integer(is), intent(inout) :: nroots1
    real(dp), intent(in)       :: tol1
    real(dp), intent(out)      :: ener1(*), xvec1(*)
    integer(is), intent(out)   :: niter1, iconv1
    real(dp)                   :: t0

    if (xdim1 /= xdim) call mrsf_error('mrsf_diag: inconsistent xdim')
    if (imult1 /= 1 .and. imult1 /= 3) call mrsf_error('mrsf_diag: mult must be 1 or 3')

    t0 = wall_time()
    call davidson_solve(imult1, irrep1, nroots1, nextra1, maxvec1, maxiter1, tol1, &
         ener1, xvec1, niter1, iconv1)
    if (verbose) write(6,'(2x,a,f10.2,a,/)') 'Davidson time: ', wall_time()-t0, ' s'

  end subroutine mrsf_diag

!######################################################################
! mrsf_sigma: A x for a batch of compressed vectors (validation)
!######################################################################
  subroutine mrsf_sigma_c(imult1, irrep1, nvec1, xdim1, x1, ax1) bind(c, name='mrsf_sigma')

    integer(is), intent(in) :: imult1, irrep1, nvec1, xdim1
    real(dp), intent(in)    :: x1(xdim1,nvec1)
    real(dp), intent(out)   :: ax1(xdim1,nvec1)

    if (xdim1 /= xdim) call mrsf_error('mrsf_sigma: inconsistent xdim')
    call sigma_batch(nvec1, imult1, irrep1, x1, ax1)

  end subroutine mrsf_sigma_c

!######################################################################
! mrsf_diagonal
!######################################################################
  subroutine mrsf_diagonal(imult1, xdim1, d1) bind(c, name='mrsf_diagonal')

    integer(is), intent(in) :: imult1, xdim1
    real(dp), intent(out)   :: d1(xdim1)

    if (xdim1 /= xdim) call mrsf_error('mrsf_diagonal: inconsistent xdim')
    call diagonal(imult1, d1)

  end subroutine mrsf_diagonal

!######################################################################
! mrsf_density: state 1-RDMs (spin-summed, MO basis)
!######################################################################
  subroutine mrsf_density_c(nmo1, occ1, imult1, xdim1, nroots1, xvec1, dmat1) &
       bind(c, name='mrsf_density')

    integer(is), intent(in) :: nmo1, imult1, xdim1, nroots1
    real(dp), intent(in)    :: occ1(nmo1), xvec1(xdim1,nroots1)
    real(dp), intent(out)   :: dmat1(nmo1,nmo1,nroots1)
    integer(is), allocatable :: ipairs(:,:)
    integer(is)              :: k

    allocate(ipairs(2,nroots1))
    do k = 1, nroots1
       ipairs(:,k) = k
    enddo
    call tdm_pairs(nmo1, occ1, imult1, xdim1, nroots1, nroots1, nroots1, ipairs, &
         xvec1, xvec1, dmat1, .true.)
    deallocate(ipairs)

  end subroutine mrsf_density_c

!######################################################################
! mrsf_tdm: 1-TDMs <bra|E_pq|ket> for pairs of states
!######################################################################
  subroutine mrsf_tdm(nmo1, occ1, imult1, xdim1, npairs1, nb1, nk1, ipairs1, xb1, xk1, &
       rho1) bind(c, name='mrsf_tdm')

    integer(is), intent(in) :: nmo1, imult1, xdim1, npairs1, nb1, nk1
    real(dp), intent(in)    :: occ1(nmo1)
    integer(is), intent(in) :: ipairs1(2,npairs1)
    real(dp), intent(in)    :: xb1(xdim1,nb1), xk1(xdim1,nk1)
    real(dp), intent(out)   :: rho1(nmo1,nmo1,npairs1)

    call tdm_pairs(nmo1, occ1, imult1, xdim1, npairs1, nb1, nk1, ipairs1, xb1, xk1, &
         rho1, .false.)

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
    init_done = .false.

  end subroutine free_orbital_data

  subroutine mrsf_finalise_f()

    call report_timings()
    call free_ints()
    call free_orbital_data()

  end subroutine mrsf_finalise_f

  subroutine report_timings()

    real(dp) :: gflop

    if (verbose .and. nsigma_calls > 0) then
       write(6,'(/,2x,a)') 'MRSF timings'
       write(6,'(2x,a,f10.2,a)') 'integral ingestion : ', time_load, ' s'
       write(6,'(2x,a,f10.2,a,i0,a,i0,a)') 'sigma vectors      : ', time_sigma, &
            ' s (', nsigma_calls, ' calls, ', nsigma_vecs, ' vectors)'
       write(6,'(2x,a,f10.2,a)') '  exchange term    : ', time_exch, ' s'
       if (time_exch > 0.0_dp) then
          gflop = 2.0_dp * real(naux,dp) * real(nvirb,dp) * real(nocca,dp) &
               * real(nvirb + nocca,dp) * real(nsigma_vecs,dp) / 1.0e9_dp
          write(6,'(2x,a,f10.1,a,f10.2,a)') '  exchange kernel  : ', gflop, &
               ' GFLOP, ', gflop / time_exch, ' GFLOP/s'
       endif
    endif

    time_load = 0.0_dp; time_sigma = 0.0_dp; time_exch = 0.0_dp
    nsigma_calls = 0; nsigma_vecs = 0

  end subroutine report_timings

end module mrsf_interface
