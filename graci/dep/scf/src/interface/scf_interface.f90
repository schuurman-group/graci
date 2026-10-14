!**********************************************************************
! scf_interface: bind(C) entry points of libscf
!   scf_set_basis     basis (PySCF _atm, _bas, _env, ao_loc, cart2sph)
!   scf_set_grid      quadrature grid, screening, AO cache
!   scf_set_operator  h, S, E_nuc, the DF tensor (address kept), the
!                     functional, the symmetry-adapted orbitals
!   scf_run           one SCF (ROKS/ROHF for nalpha > nbeta, else RKS/RHF)
!   scf_fock_eval     Fock matrices, energy, Roothaan matrix and gradient
!                     norm of given orbitals (tests)
!   scf_gto_eval      AO values of all shells (tests)
!   scf_timings, scf_free
! Arrays are passed by address (F-ordered float64 from Python's 'dptr'
! or copied int32 arrays); ierr = 0 on success.
!**********************************************************************
module scf_interface

  use iso_c_binding
  use scf_constants
  use scf_io
  use scf_gto
  use scf_xc
  use scf_grid
  use scf_df
  use scf_fock
  use scf_solver

  implicit none

  character(len=255) :: basis_stamp = ''
  logical            :: op_set = .false.

contains

!######################################################################
! scf_set_basis: atm(6,natm), bas(8,nbas), env(nenv) (0-based pointers),
! aoloc(nbas+1) (0-based offsets), c2sm(nc2s): cart2sph(l) column-major
! for l = 0..lmax starting at c2soff(l)+1; the basis is kept until the
! stamp changes (a new basis also discards the grid)
!######################################################################
  subroutine scf_set_basis(natm, nbas, nenv, atm, bas, env, nao, aoloc, cart, lmax, nc2s, c2sm, c2soff, &
       stamp1, ierr) bind(c, name='scf_set_basis')

    integer(is), intent(in)            :: natm, nbas, nenv, nao, lmax, nc2s
    integer(is), intent(in)            :: atm(6,natm), bas(8,nbas), aoloc(nbas+1), c2soff(0:lmax)
    real(dp), intent(in)               :: env(nenv), c2sm(nc2s)
    logical(c_bool), intent(in)        :: cart
    character(kind=c_char), intent(in) :: stamp1(*)
    integer(is), intent(out)           :: ierr
    character(len=255)                 :: stamp
    logical                            :: lcart

    ierr = 0
    call c2fstr(stamp1, stamp)
    if (stamp == basis_stamp .and. nbas_b == nbas .and. nao_b == nao .and. len_trim(stamp) > 0) return
    lcart = cart
    call grid_free()
    call gto_set(natm, nbas, nenv, atm, bas, env, nao, aoloc, lcart, lmax, nc2s, c2sm, c2soff, ierr)
    if (ierr == 0) then
       basis_stamp = stamp
    else
       basis_stamp = ''
    endif

  end subroutine scf_set_basis

!######################################################################
! scf_set_grid: xyz(3,ng) (PySCF's C-ordered coordinates), weights,
! ncomp 4 (GGA) or 1 (LDA), budget (GB) for the AO cache; reused when
! the stamp is unchanged; ierr 3: no basis
!######################################################################
  subroutine scf_set_grid(ng, xyz, w, ncomp, budget, stamp1, ierr) bind(c, name='scf_set_grid')

    integer(is), intent(in)            :: ng, ncomp
    real(dp), intent(in)               :: xyz(3,ng), w(ng), budget
    character(kind=c_char), intent(in) :: stamp1(*)
    integer(is), intent(out)           :: ierr
    character(len=255)                 :: stamp

    ierr = 0
    if (nbas_b == 0) then
       ierr = 3
       return
    endif
    call c2fstr(stamp1, stamp)
    if (stamp == gstamp .and. ngrid == ng .and. ncomp_g == ncomp .and. len_trim(stamp) > 0) return
    call grid_set(ng, xyz, w, ncomp, budget*1.0e9_dp, stamp, ierr)

  end subroutine scf_set_grid

!######################################################################
! scf_set_operator: h, S (nao,nao), E_nuc, the DF tensor B(npair,naux)
! (the C-ordered (naux, npair) PySCF array; its address is kept and must
! stay valid while the library uses it), c_x, the functional (libxc ids,
! families 1 LDA / 2 GGA, coefficients; nfunc = 0 for HF), nspin (2 for
! ROKS, 1 for RKS), the symmetry-adapted orbitals so(nao,nao) in irrep
! blocks of nso(i) columns with PySCF irrep ids irrid(i)
! ierr 4: nao differs from the basis; 5: functional; 6: symmetry blocks
!######################################################################
  subroutine scf_set_operator(nao, naux, h, S, enuc, B, cx, nfunc, ids, fam, coefs, nspin, &
       nirr1, nso, irrid, so, ierr) bind(c, name='scf_set_operator')

    integer(is), intent(in)       :: nao, naux, nfunc, nspin, nirr1
    integer(is), intent(in)       :: ids(*), fam(*), nso(nirr1), irrid(nirr1)
    real(dp), intent(in)          :: h(nao,nao), S(nao,nao), enuc, cx, coefs(*)
    real(dp), intent(in), target  :: B(nao*(nao+1)/2, naux)
    real(dp), intent(in)          :: so(nao,nao)
    integer(is), intent(out)      :: ierr
    integer(is)                   :: ie

    ierr = 0
    op_set = .false.
    if (nbas_b > 0 .and. nao /= nao_b) then
       ierr = 4
       return
    endif
    call fock_set(nao, h, S, enuc, cx, nfunc > 0)
    call df_set(nao, naux, c_loc(B))
    call xc_setup(nfunc, ids, fam, coefs, nspin, ie)
    if (ie /= 0) then
       ierr = 5
       return
    endif
    call solver_set_symmetry(nao, nirr1, nso, irrid, so, S, ie)
    if (ie /= 0) then
       ierr = 6
       return
    endif
    op_set = .true.

  end subroutine scf_set_operator

!######################################################################
! checks shared by scf_run and scf_fock_eval; ierr 7: no operator, 8: XC
! without a grid, 9: electron counts, 10: functional spin
!######################################################################
  subroutine check_ready(nao, na, nb, ierr)
    integer(is), intent(in)  :: nao, na, nb
    integer(is), intent(out) :: ierr
    ierr = 0
    if (.not. op_set .or. nao /= nao_f) then
       ierr = 7
    else if (have_xc .and. ngrid == 0) then
       ierr = 8
    else if (na < nb .or. nb < 0 .or. na > nao) then
       ierr = 9
    else if (have_xc .and. ((na > nb .and. xc_nspin /= 2) .or. (na == nb .and. xc_nspin /= 1))) then
       ierr = 10
    endif
  end subroutine check_ready

!######################################################################
! scf_run: one SCF from the start factor X0(nao,k0a) (alpha density
! X0 X0^T, beta density from its first k0b columns; RKS: D = 2 X0 X0^T)
!   opts: conv_tol, conv_tol_grad, max_cycle, diis_space, diis_start,
!         level_shift, exact_until, verbose
!   outputs: C, mo_energy, mo_occ, mo_irrep (PySCF's order), Fa, Fb (the
!         spin Fock matrices of the final density; RKS: Fa = Fb = F),
!   ener: E, E_1e, E_coul, E_exch, E_xc, E_nuc, |g|, dE of the last cycle
!   info: converged, cycles, shortcut cycles, occupation rule ok
!######################################################################
  subroutine scf_run(nao, na, nb, k0a, k0b, X0, opts, C, e, occ, irr, Fa, Fb, ener, info, ierr) &
       bind(c, name='scf_run')

    integer(is), intent(in)  :: nao, na, nb, k0a, k0b
    real(dp), intent(in)     :: X0(nao,k0a), opts(8)
    real(dp), intent(out)    :: C(nao,nao), e(nao), occ(nao), Fa(nao,nao), Fb(nao,nao), ener(8)
    integer(is), intent(out) :: irr(nao), info(8), ierr
    real(dp)                 :: Etot

    info = 0
    ener = 0.0_dp
    call check_ready(nao, na, nb, ierr)
    if (ierr /= 0) return
    if (k0b > k0a .or. k0a < 1) then
       ierr = 9
       return
    endif
    conv_tol      = opts(1)
    conv_tol_grad = opts(2)
    max_cycle     = max(1, nint(opts(3)))
    diis_space    = max(1, nint(opts(4)))
    diis_start    = nint(opts(5))
    level_shift   = opts(6)
    exact_until   = nint(opts(7))
    verbose       = opts(8) > 0.5_dp

    if (na > nb) then
       call scf_iterate_roks(nao, na, nb, k0a, k0b, X0, C, e, occ, irr, Fa, Fb, Etot)
    else
       call scf_iterate_rks(nao, na, k0a, X0, C, e, occ, irr, Fa, Etot)
       Fb = Fa
    endif
    ener(1:6) = ecomp
    ener(7) = g_final
    ener(8) = de_final
    info(1) = merge(1, 0, converged)
    info(2) = ncycle
    info(3) = nshort
    info(4) = merge(1, 0, occ_rule_ok)

  end subroutine scf_run

!######################################################################
! scf_fock_eval: for orbitals C with occupations occ (2/1/0) and irrep
! ids irr: mode 0 ROKS exact (Fa, Fb, Fc, dCO = (Fa - Fb) C_O, R from the
! projector formula, exact gradient), mode 1 ROKS shortcut (Fc, dCO, R
! and gradient from F_c and dCO), mode 2 RKS (Fa = F, R = F)
!   ener: E, E_1e, E_coul, E_exch, E_xc, E_nuc, |g|
!######################################################################
  subroutine scf_fock_eval(mode, nao, C, occ, irr, ener, Fa, Fb, Fc, dCO, R, ierr) bind(c, name='scf_fock_eval')

    integer(is), intent(in)  :: mode, nao, irr(nao)
    real(dp), intent(in)     :: C(nao,nao), occ(nao)
    real(dp), intent(out)    :: ener(8), Fa(nao,nao), Fb(nao,nao), Fc(nao,nao), dCO(nao,2), R(nao,nao)
    integer(is), intent(out) :: ierr

    real(dp), allocatable :: X(:,:), Da(:,:), Db(:,:), fcmo(:,:), dmo(:,:)
    integer(is) :: na, nb, i, k
    real(dp)    :: Etot, gn

    ener = 0.0_dp
    nb = count(occ > 1.5_dp)
    na = nb + count(abs(occ - 1.0_dp) < 1.0e-6_dp)
    call check_ready(nao, na, nb, ierr)
    if (ierr /= 0) return
    if (na - nb > 2) then
       ierr = 9
       return
    endif
    allocate(X(nao,max(na,1)), Da(nao,nao), Db(nao,nao), fcmo(nao,nao), dmo(nao,2))
    k = 0
    do i = 1, nao
       if (occ(i) > 1.5_dp) then
          k = k + 1
          X(:,k) = C(:,i)
       endif
    enddo
    do i = 1, nao
       if (abs(occ(i) - 1.0_dp) < 1.0e-6_dp) then
          k = k + 1
          X(:,k) = C(:,i)
       endif
    enddo
    Fa = 0.0_dp; Fb = 0.0_dp; Fc = 0.0_dp; dCO = 0.0_dp; R = 0.0_dp
    select case (mode)
    case (0)
       call build_roks(nao, na, nb, X, 0, Etot, Fa, Fb, Fc, dCO, Da, Db)
       call roothaan_ao(nao, Fa, Fb, Da, Db, R)
       gn = gnorm_roks_exact(nao, C, occ, irr, Fa, Fb)
    case (1)
       call build_roks(nao, na, nb, X, 1, Etot, Fa, Fb, Fc, dCO, Da, Db)
       call roothaan_short(nao, C, occ, Fc, dCO, R, fcmo, dmo)
       gn = gnorm_roks_short(nao, occ, irr, fcmo, dmo)
    case default
       call build_rks(nao, nb, X, Etot, Fa, Da)
       R = Fa
       gn = gnorm_rks(nao, C, occ, irr, Fa)
    end select
    ener(1:6) = ecomp
    ener(7) = gn
    deallocate(X, Da, Db, fcmo, dmo)

  end subroutine scf_fock_eval

!######################################################################
! scf_gto_eval: values (deriv 0) or values and gradients (deriv 1) of all
! AOs at the points xyz(3,np): ao(np, nao, 1+3*deriv)
!######################################################################
  subroutine scf_gto_eval(np, xyz, deriv, ao, ierr) bind(c, name='scf_gto_eval')

    integer(is), intent(in)  :: np, deriv
    real(dp), intent(in)     :: xyz(3,np)
    real(dp), intent(out)    :: ao(np,nao_b,*)
    integer(is), intent(out) :: ierr

    integer(is), allocatable :: shl(:)
    real(dp), allocatable    :: x(:), y(:), z(:)
    integer(is) :: i

    ierr = 0
    if (nbas_b == 0) then
       ierr = 3
       return
    endif
    allocate(shl(nbas_b), x(np), y(np), z(np))
    do i = 1, nbas_b
       shl(i) = i
    enddo
    x = xyz(1,:); y = xyz(2,:); z = xyz(3,:)
    call gto_eval(np, x, y, z, nbas_b, shl, deriv, np, nao_b, ao)
    deallocate(shl, x, y, z)

  end subroutine scf_gto_eval

!######################################################################
! scf_timings: J, K (s), K flops, XC passes (s), AO cache fill (s), flops
! of one V product over the grid, significant-AO fraction, V cost
! fraction, diagonalisations, matrix work, SCF total, grid points,
! blocks, AO cache used (1/0), max significant AOs per block
!######################################################################
  subroutine scf_timings(t) bind(c, name='scf_timings')
    real(dp), intent(out) :: t(16)
    t = 0.0_dp
    t(1)  = time_j
    t(2)  = time_k
    t(3)  = flop_k
    t(4)  = time_xcpass
    t(5)  = time_aofill
    t(6)  = flop_xc
    t(7)  = frac_sig
    t(8)  = cost_sig
    t(9)  = time_diag
    t(10) = time_mat
    t(11) = time_total
    t(12) = real(ngrid, dp)
    t(13) = real(nblk, dp)
    t(14) = merge(1.0_dp, 0.0_dp, cached)
    t(15) = real(nsigmax, dp)
    t(16) = time_y
  end subroutine scf_timings

  subroutine scf_xc_subtimers(t) bind(c, name='scf_xc_subtimers')
    real(dp), intent(out) :: t(6)
    t = tsub
  end subroutine scf_xc_subtimers

  subroutine scf_reset_timings() bind(c, name='scf_reset_timings')
    time_j = 0.0_dp; time_k = 0.0_dp; flop_k = 0.0_dp; time_y = 0.0_dp
    time_xcpass = 0.0_dp; time_diag = 0.0_dp; time_mat = 0.0_dp; time_total = 0.0_dp
    tsub = 0.0_dp
  end subroutine scf_reset_timings

!######################################################################
! scf_free: release everything
!######################################################################
  subroutine scf_free() bind(c, name='scf_free')
    call grid_free()
    call gto_free()
    call df_free()
    call fock_free()
    call xc_release()
    call solver_free()
    basis_stamp = ''
    op_set = .false.
  end subroutine scf_free

end module scf_interface
