!**********************************************************************
! scf_solver: the SCF iterations (ROKS/ROHF and RKS/RHF)
!   - per-irrep S-orthonormal bases from the symmetry-adapted orbitals
!     of PySCF (one block for C1); diagonalisation per irrep, orbitals
!     irrep-pure by construction
!   - CDIIS on the Roothaan (ROKS) or Fock (RKS) matrix, error
!     X^T (F D S - S D F) X in the orthonormal basis, D = D_a + D_b
!   - occupations: PySCF's ROHF rule (closed: the nb lowest Roothaan
!     eigenvalues; open: the nopen lowest alpha energies among the rest)
!     while the occupation pattern settles (exact Fock builds), then the
!     open orbitals followed by maximum overlap and the cheaper shortcut
!     build; the rule is checked exactly at convergence
!   - convergence: |dE| < conv_tol and |g| < conv_tol_grad (PySCF's
!     orbital gradient); final step: exact build at the final orbitals,
!     canonicalisation within (irrep x closed/open/virtual) blocks, PySCF's
!     ordering (C1: ascending eigenvalue; symmetry: closed | open |
!     virtual, each sorted by the eigenvalue rounded to 1e-9)
!**********************************************************************
module scf_solver

  use scf_constants
  use scf_io
  use scf_df
  use scf_grid
  use scf_fock

  implicit none

  integer(is)              :: nirr = 0
  integer(is), allocatable :: irr_n(:), irr_id(:), irr_off(:)
  real(dp), allocatable    :: Xorth(:,:)

  ! options (set before scf_iterate)
  real(dp)    :: conv_tol = 1.0e-8_dp, conv_tol_grad = 1.0e-4_dp, level_shift = 0.0_dp
  integer(is) :: max_cycle = 50, diis_space = 14, diis_start = 1, exact_until = 3

  ! results
  logical     :: converged = .false., occ_rule_ok = .true.
  integer(is) :: ncycle = 0, nshort = 0
  real(dp)    :: g_final = 0.0_dp, de_final = 0.0_dp, time_diag = 0.0_dp, time_total = 0.0_dp

  external :: dsyevd, dgesv

contains

!######################################################################
! solver_set_symmetry: so(nao,nao) = the symmetry-adapted orbitals of
! the irreps one after the other (n(i) columns each, PySCF irrep ids
! ids(i)); the identity and one irrep without symmetry
!######################################################################
  subroutine solver_set_symmetry(nao, nir, n, ids, so, S, ierr)

    integer(is), intent(in)  :: nao, nir, n(nir), ids(nir)
    real(dp), intent(in)     :: so(nao,nao), S(nao,nao)
    integer(is), intent(out) :: ierr

    real(dp), allocatable :: s_ir(:,:), w(:), T(:,:)
    integer(is) :: i, k, m, j

    ierr = 0
    if (allocated(irr_n)) deallocate(irr_n, irr_id, irr_off, Xorth)
    nirr = nir
    allocate(irr_n(nir), irr_id(nir), irr_off(nir+1), Xorth(nao,nao))
    irr_n = n
    irr_id = ids
    irr_off(1) = 0
    do i = 1, nir
       irr_off(i+1) = irr_off(i) + n(i)
    enddo
    if (irr_off(nir+1) /= nao) then
       ierr = 1
       return
    endif
    do i = 1, nir
       m = n(i)
       if (m == 0) cycle
       k = irr_off(i)
       allocate(T(nao,m), s_ir(m,m), w(m))
       call dgemm('N', 'N', nao, m, nao, 1.0_dp, S, nao, so(1,k+1), nao, 0.0_dp, T, nao)
       call dgemm('T', 'N', m, m, nao, 1.0_dp, so(1,k+1), nao, T, nao, 0.0_dp, s_ir, m)
       call eigh(m, s_ir, w, ierr)
       if (ierr /= 0) return
       if (minval(w) <= 0.0_dp) then
          ierr = 2
          return
       endif
       do j = 1, m
          s_ir(:,j) = s_ir(:,j)/sqrt(w(j))
       enddo
       call dgemm('N', 'N', nao, m, m, 1.0_dp, so(1,k+1), nao, s_ir, m, 0.0_dp, Xorth(1,k+1), nao)
       deallocate(T, s_ir, w)
    enddo

  end subroutine solver_set_symmetry

  subroutine solver_free()
    if (allocated(irr_n)) deallocate(irr_n, irr_id, irr_off, Xorth)
    nirr = 0
  end subroutine solver_free

!######################################################################
! eigh: symmetric eigenproblem in place (eigenvectors in A)
!######################################################################
  subroutine eigh(n, A, w, ierr)
    integer(is), intent(in)  :: n
    real(dp), intent(inout)  :: A(n,n)
    real(dp), intent(out)    :: w(n)
    integer(is), intent(out) :: ierr
    real(dp), allocatable    :: work(:)
    integer(is), allocatable :: iwork(:)
    integer(is) :: lwork, liwork
    real(dp)    :: wq(1)
    integer(is) :: iq(1)
    call dsyevd('V', 'U', n, A, n, w, wq, -1, iq, -1, ierr)
    lwork = int(wq(1)); liwork = iq(1)
    allocate(work(lwork), iwork(liwork))
    call dsyevd('V', 'U', n, A, n, w, work, lwork, iwork, liwork, ierr)
    deallocate(work, iwork)
  end subroutine eigh

!######################################################################
! diag_irrep: eigenpairs of F per irrep block, sorted by eigenvalue
! (stable); irr = PySCF irrep id of each orbital
!######################################################################
  subroutine diag_irrep(nao, F, e, C, irr)

    integer(is), intent(in) :: nao
    real(dp), intent(in)    :: F(nao,nao)
    real(dp), intent(out)   :: e(nao), C(nao,nao)
    integer(is), intent(out):: irr(nao)

    real(dp), allocatable    :: T(:,:), Fp(:,:), w(:), e0(:), C0(:,:)
    integer(is), allocatable :: ir0(:), perm(:)
    integer(is) :: i, k, m, ierr
    real(dp) :: t0

    t0 = wall_time()
    allocate(e0(nao), C0(nao,nao), ir0(nao), perm(nao))
    do i = 1, nirr
       m = irr_n(i)
       if (m == 0) cycle
       k = irr_off(i)
       allocate(T(nao,m), Fp(m,m), w(m))
       call dgemm('N', 'N', nao, m, nao, 1.0_dp, F, nao, Xorth(1,k+1), nao, 0.0_dp, T, nao)
       call dgemm('T', 'N', m, m, nao, 1.0_dp, Xorth(1,k+1), nao, T, nao, 0.0_dp, Fp, m)
       call eigh(m, Fp, w, ierr)
       call dgemm('N', 'N', nao, m, m, 1.0_dp, Xorth(1,k+1), nao, Fp, m, 0.0_dp, C0(1,k+1), nao)
       e0(k+1:k+m) = w
       ir0(k+1:k+m) = irr_id(i)
       deallocate(T, Fp, w)
    enddo
    call stable_sort(nao, e0, perm)
    do i = 1, nao
       e(i) = e0(perm(i))
       C(:,i) = C0(:,perm(i))
       irr(i) = ir0(perm(i))
    enddo
    deallocate(e0, C0, ir0, perm)
    time_diag = time_diag + wall_time() - t0

  end subroutine diag_irrep

!######################################################################
! stable_sort: perm such that key(perm) is ascending, equal keys kept in
! their original order (merge sort)
!######################################################################
  subroutine stable_sort(n, key, perm)
    integer(is), intent(in)  :: n
    real(dp), intent(in)     :: key(n)
    integer(is), intent(out) :: perm(n)
    integer(is), allocatable :: tmp(:)
    integer(is) :: width, lo, mid, hi, i, j, k
    allocate(tmp(n))
    do i = 1, n
       perm(i) = i
    enddo
    width = 1
    do while (width < n)
       lo = 1
       do while (lo <= n)
          mid = min(lo + width - 1, n)
          hi  = min(lo + 2*width - 1, n)
          i = lo; j = mid + 1; k = lo
          do while (i <= mid .and. j <= hi)
             if (key(perm(j)) < key(perm(i))) then
                tmp(k) = perm(j); j = j + 1
             else
                tmp(k) = perm(i); i = i + 1
             endif
             k = k + 1
          enddo
          do while (i <= mid)
             tmp(k) = perm(i); i = i + 1; k = k + 1
          enddo
          do while (j <= hi)
             tmp(k) = perm(j); j = j + 1; k = k + 1
          enddo
          lo = lo + 2*width
       enddo
       perm = tmp
       width = 2*width
    enddo
    deallocate(tmp)
  end subroutine stable_sort

!######################################################################
! DIIS: store (F, err) and return the extrapolated matrix
!######################################################################
  subroutine diis_step(nao, F, D, hist_F, hist_e, nh, Fx)

    integer(is), intent(in)    :: nao
    real(dp), intent(in)       :: F(nao,nao), D(nao,nao)
    real(dp), intent(inout)    :: hist_F(nao,nao,*), hist_e(nao,nao,*)
    integer(is), intent(inout) :: nh
    real(dp), intent(out)      :: Fx(nao,nao)

    real(dp), allocatable    :: T1(:,:), T2(:,:), Bm(:,:), c(:), Bv(:,:), w(:), g(:)
    integer(is), allocatable :: ipiv(:)
    integer(is) :: i, k, n, info

    allocate(T1(nao,nao), T2(nao,nao))
    ! err = X^T (F D S - S D F) X
    call dgemm('N', 'N', nao, nao, nao, 1.0_dp, F, nao, D, nao, 0.0_dp, T1, nao)
    call dgemm('N', 'N', nao, nao, nao, 1.0_dp, T1, nao, Sm, nao, 0.0_dp, T2, nao)
    T1 = T2 - transpose(T2)
    call dgemm('N', 'N', nao, nao, nao, 1.0_dp, T1, nao, Xorth, nao, 0.0_dp, T2, nao)
    if (nh == diis_space) then
       do i = 1, nh - 1
          hist_F(:,:,i) = hist_F(:,:,i+1)
          hist_e(:,:,i) = hist_e(:,:,i+1)
       enddo
       nh = nh - 1
    endif
    nh = nh + 1
    hist_F(:,:,nh) = F
    call dgemm('T', 'N', nao, nao, nao, 1.0_dp, Xorth, nao, T2, nao, 0.0_dp, hist_e(1,1,nh), nao)
    n = nh
    allocate(Bm(n+1,n+1), c(n+1), ipiv(n+1))
    Bm = 0.0_dp
    Bm(1,2:n+1) = 1.0_dp
    Bm(2:n+1,1) = 1.0_dp
    do i = 1, n
       do k = 1, i
          Bm(i+1,k+1) = sum(hist_e(:,:,i)*hist_e(:,:,k))
          Bm(k+1,i+1) = Bm(i+1,k+1)
       enddo
    enddo
    ! PySCF's rule: solve, or the pseudo-inverse without the eigenvalues
    ! below 1e-14 in magnitude
    allocate(Bv(n+1,n+1), w(n+1), g(n+1))
    Bv = Bm
    call eigh(n+1, Bv, w, info)
    g = 0.0_dp
    g(1) = 1.0_dp
    if (any(abs(w) < 1.0e-14_dp) .or. info /= 0) then
       c = 0.0_dp
       do i = 1, n+1
          if (abs(w(i)) > 1.0e-14_dp) c = c + Bv(:,i)*(Bv(1,i)/w(i))
       enddo
    else
       c = g
       call dgesv(n+1, 1, Bm, n+1, ipiv, c, n+1, info)
    endif
    deallocate(Bv, w, g)
    Fx = 0.0_dp
    do i = 1, n
       Fx = Fx + c(i+1)*hist_F(:,:,i)
    enddo
    deallocate(T1, T2, Bm, c, ipiv)

  end subroutine diis_step

!######################################################################
! rohf_occupations: PySCF's rule; ea = alpha energies of the orbitals
!######################################################################
  subroutine rohf_occupations(nao, na, nb, e, ea, occ)
    integer(is), intent(in) :: nao, na, nb
    real(dp), intent(in)    :: e(nao), ea(nao)
    real(dp), intent(out)   :: occ(nao)
    integer(is), allocatable :: order(:), rest(:), o2(:)
    real(dp), allocatable    :: er(:)
    integer(is) :: i, nr
    allocate(order(nao), rest(nao), o2(nao), er(nao))
    call stable_sort(nao, e, order)
    occ = 0.0_dp
    do i = 1, nb
       occ(order(i)) = 2.0_dp
    enddo
    nr = nao - nb
    do i = 1, nr
       rest(i) = order(nb+i)
       er(i) = ea(rest(i))
    enddo
    call stable_sort(nr, er(1:nr), o2(1:nr))
    do i = 1, na - nb
       occ(rest(o2(i))) = 1.0_dp
    enddo
    deallocate(order, rest, o2, er)
  end subroutine rohf_occupations

!######################################################################
! scf_iterate_roks: the ROKS/ROHF SCF from the start factor X0 (nao,ka0)
! (alpha density X0 X0^T, beta density X0(:,1:kb0) X0(:,1:kb0)^T);
! outputs C, mo_energy, mo_occ, mo_irrep (ordered), Fa, Fb, E
!######################################################################
  subroutine scf_iterate_roks(nao, na, nb, ka0, kb0, X0, C, e, occ, irr, Fa, Fb, Etot)

    integer(is), intent(in)  :: nao, na, nb, ka0, kb0
    real(dp), intent(in)     :: X0(nao,ka0)
    real(dp), intent(out)    :: C(nao,nao), e(nao), occ(nao), Fa(nao,nao), Fb(nao,nao), Etot
    integer(is), intent(out) :: irr(nao)

    real(dp), allocatable    :: X(:,:), Fc(:,:), dCO(:,:), Da(:,:), Db(:,:), R(:,:), Rx(:,:), D(:,:)
    real(dp), allocatable    :: hF(:,:,:), he(:,:,:), Cn(:,:), en(:), ea(:), occn(:), fcmo(:,:), dmo(:,:)
    real(dp), allocatable    :: SCo(:,:), ov(:,:), T(:,:), sc(:)
    integer(is), allocatable :: irn(:), pat(:), patl(:), order(:), rest(:), o2(:)
    real(dp)    :: E_last, dE, gn, t0
    integer(is) :: cyc, nh, i, k, no, stable, nr, attempt
    logical     :: exact, have_orb, last_exact

    t0 = wall_time()
    no = na - nb
    allocate(X(nao,na), Fc(nao,nao), dCO(nao,max(no,1)), Da(nao,nao), Db(nao,nao), R(nao,nao), &
         Rx(nao,nao), D(nao,nao), hF(nao,nao,diis_space), he(nao,nao,diis_space), Cn(nao,nao), &
         en(nao), ea(nao), occn(nao), irn(nao), fcmo(nao,nao), dmo(nao,max(no,1)), &
         pat(2*max(no,1)), patl(2*max(no,1)), order(nao), rest(nao), o2(nao), SCo(nao,max(no,1)), &
         ov(max(no,1),nao), T(nao,nao), sc(nao))

    ! first Fock build from the start factor (exact)
    call build_roks(nao, ka0, kb0, X0, 0, Etot, Fa, Fb, Fc, dCO, Da, Db)
    call roothaan_ao(nao, Fa, Fb, Da, Db, R)
    E_last = Etot
    dE = 0.0_dp
    last_exact = .true.
    nh = 0
    exact = .true.
    stable = 0
    patl = -1
    have_orb = .false.
    converged = .false.
    nshort = 0
    ncycle = 0
    occ = 0.0_dp

    do attempt = 1, 2
       do cyc = ncycle + 1, max_cycle
          ncycle = cyc
          D = Da + Db
          if (cyc >= diis_start) then
             call diis_step(nao, R, D, hF, he, nh, Rx)
          else
             Rx = R
          endif
          if (level_shift > 0.0_dp) then
             ! PySCF's level shift with D/2: virtual by s, open by s/2
             call dgemm('N', 'N', nao, nao, nao, 1.0_dp, Sm, nao, D, nao, 0.0_dp, T, nao)
             call dgemm('N', 'N', nao, nao, nao, -0.5_dp*level_shift, T, nao, Sm, nao, 1.0_dp, Rx, nao)
             Rx = Rx + level_shift*Sm
          endif
          call diag_irrep(nao, Rx, en, Cn, irn)
          ! occupations
          if (exact .or. .not. have_orb) then
             call dgemm('N', 'N', nao, nao, nao, 1.0_dp, Fa, nao, Cn, nao, 0.0_dp, T, nao)
             do i = 1, nao
                ea(i) = dot_product(Cn(:,i), T(:,i))
             enddo
             call rohf_occupations(nao, na, nb, en, ea, occn)
          else
             ! closed: the nb lowest eigenvalues; open: largest overlap with the previous open orbitals
             call stable_sort(nao, en, order)
             occn = 0.0_dp
             do i = 1, nb
                occn(order(i)) = 2.0_dp
             enddo
             nr = nao - nb
             call dgemm('T', 'N', no, nao, nao, 1.0_dp, SCo, nao, Cn, nao, 0.0_dp, ov, max(no,1))
             do i = 1, nr
                rest(i) = order(nb+i)
                sc(i) = -sum(ov(1:no,rest(i))**2)
             enddo
             call stable_sort(nr, sc, o2)
             do i = 1, no
                occn(rest(o2(i))) = 1.0_dp
             enddo
          endif
          ! pattern of the open orbitals (index and irrep)
          k = 0
          do i = 1, nao
             if (abs(occn(i) - 1.0_dp) < 1.0e-6_dp) then
                k = k + 1
                pat(k) = i
                pat(no+k) = irn(i)
             endif
          enddo
          if (all(pat(1:2*no) == patl(1:2*no))) then
             stable = stable + 1
          else
             stable = 0
          endif
          patl = pat
          C = Cn
          occ = occn
          irr = irn
          have_orb = .true.
          if (exact .and. stable >= exact_until) exact = .false.
          ! factor [closed | open]
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
          if (no > 0) then
             call dgemm('N', 'N', nao, no, nao, 1.0_dp, Sm, nao, X(1,nb+1), nao, 0.0_dp, SCo, nao)
          endif
          if (exact) then
             call build_roks(nao, na, nb, X, 0, Etot, Fa, Fb, Fc, dCO, Da, Db)
             call roothaan_ao(nao, Fa, Fb, Da, Db, R)
             gn = gnorm_roks_exact(nao, C, occ, irr, Fa, Fb)
          else
             call build_roks(nao, na, nb, X, 1, Etot, Fa, Fb, Fc, dCO, Da, Db)
             call roothaan_short(nao, C, occ, Fc, dCO, R, fcmo, dmo)
             gn = gnorm_roks_short(nao, occ, irr, fcmo, dmo)
             nshort = nshort + 1
          endif
          last_exact = exact
          dE = Etot - E_last
          E_last = Etot
          if (verbose) then
             write(6,'(2x,a,i4,a,f20.12,a,es10.2,a,es10.2,2x,a)') 'cycle ', cyc, '  E = ', Etot, '  dE = ', dE, &
                  '  |g| = ', gn, merge('exact   ', 'shortcut', exact)
             flush(6)
          endif
          if (abs(dE) < conv_tol .and. gn < conv_tol_grad) then
             converged = .true.
             exit
          endif
       enddo

       ! final step: exact build (unless the last cycle made one),
       ! canonicalisation, occupation check
       if (.not. last_exact) then
          call build_roks(nao, na, nb, X, 0, Etot, Fa, Fb, Fc, dCO, Da, Db)
          call roothaan_ao(nao, Fa, Fb, Da, Db, R)
          last_exact = .true.
       endif
       call canonicalise(nao, C, occ, irr, R, e)
       call dgemm('N', 'N', nao, nao, nao, 1.0_dp, Fa, nao, C, nao, 0.0_dp, T, nao)
       do i = 1, nao
          ea(i) = dot_product(C(:,i), T(:,i))
       enddo
       call rohf_occupations(nao, na, nb, e, ea, occn)
       occ_rule_ok = all(abs(occn - occ) < 1.0e-6_dp)
       if (occ_rule_ok .or. .not. converged .or. ncycle >= max_cycle) exit
       ! the rule picks other orbitals: continue with exact builds
       exact = .true.
       stable = 0
       converged = .false.
    enddo

    g_final = gnorm_roks_exact(nao, C, occ, irr, Fa, Fb)
    de_final = dE
    call order_orbitals(nao, C, e, occ, irr)
    time_total = time_total + wall_time() - t0
    deallocate(X, Fc, dCO, Da, Db, R, Rx, D, hF, he, Cn, en, ea, occn, irn, fcmo, dmo, pat, patl, &
         order, rest, o2, SCo, ov, T, sc)

  end subroutine scf_iterate_roks

!######################################################################
! scf_iterate_rks: closed-shell SCF from the start factor X0 (D = 2 X0 X0^T)
!######################################################################
  subroutine scf_iterate_rks(nao, nocc, k0, X0, C, e, occ, irr, F, Etot)

    integer(is), intent(in)  :: nao, nocc, k0
    real(dp), intent(in)     :: X0(nao,k0)
    real(dp), intent(out)    :: C(nao,nao), e(nao), occ(nao), F(nao,nao), Etot
    integer(is), intent(out) :: irr(nao)

    real(dp), allocatable    :: X(:,:), D(:,:), Fx(:,:), hF(:,:,:), he(:,:,:), en(:)
    integer(is), allocatable :: order(:)
    real(dp)    :: E_last, dE, gn, t0
    integer(is) :: cyc, nh, i, k

    t0 = wall_time()
    allocate(X(nao,nocc), D(nao,nao), Fx(nao,nao), hF(nao,nao,diis_space), he(nao,nao,diis_space), &
         en(nao), order(nao))
    call build_rks(nao, k0, X0, Etot, F, D)
    E_last = Etot
    dE = 0.0_dp
    nh = 0
    converged = .false.
    ncycle = 0
    nshort = 0
    do cyc = 1, max_cycle
       ncycle = cyc
       if (cyc >= diis_start) then
          call diis_step(nao, F, D, hF, he, nh, Fx)
       else
          Fx = F
       endif
       call diag_irrep(nao, Fx, e, C, irr)
       call stable_sort(nao, e, order)
       occ = 0.0_dp
       do i = 1, nocc
          occ(order(i)) = 2.0_dp
       enddo
       k = 0
       do i = 1, nao
          if (occ(i) > 1.5_dp) then
             k = k + 1
             X(:,k) = C(:,i)
          endif
       enddo
       call build_rks(nao, nocc, X, Etot, F, D)
       gn = gnorm_rks(nao, C, occ, irr, F)
       dE = Etot - E_last
       E_last = Etot
       if (verbose) then
          write(6,'(2x,a,i4,a,f20.12,a,es10.2,a,es10.2)') 'cycle ', cyc, '  E = ', Etot, '  dE = ', dE, '  |g| = ', gn
          flush(6)
       endif
       if (abs(dE) < conv_tol .and. gn < conv_tol_grad) then
          converged = .true.
          exit
       endif
    enddo
    call canonicalise(nao, C, occ, irr, F, e)
    g_final = gnorm_rks(nao, C, occ, irr, F)
    de_final = dE
    occ_rule_ok = .true.
    call order_orbitals(nao, C, e, occ, irr)
    time_total = time_total + wall_time() - t0
    deallocate(X, D, Fx, hF, he, en, order)

  end subroutine scf_iterate_rks

!######################################################################
! canonicalise: diagonalise R within the (irrep x occupation class)
! blocks of the orbitals C; e = block eigenvalues
!######################################################################
  subroutine canonicalise(nao, C, occ, irr, R, e)

    integer(is), intent(in) :: nao, irr(nao)
    real(dp), intent(inout) :: C(nao,nao)
    real(dp), intent(in)    :: occ(nao), R(nao,nao)
    real(dp), intent(out)   :: e(nao)

    real(dp), allocatable    :: Cb(:,:), T(:,:), Rb(:,:), w(:)
    integer(is), allocatable :: idx(:)
    integer(is) :: cls, i, j, n, ierr
    real(dp)    :: oc

    allocate(idx(nao))
    e = 0.0_dp
    do cls = 0, 2
       oc = real(cls, dp)
       do i = 1, nirr
          n = 0
          do j = 1, nao
             if (abs(occ(j) - oc) < 1.0e-6_dp .and. irr(j) == irr_id(i)) then
                n = n + 1
                idx(n) = j
             endif
          enddo
          if (n == 0) cycle
          allocate(Cb(nao,n), T(nao,n), Rb(n,n), w(n))
          Cb = C(:,idx(1:n))
          call dgemm('N', 'N', nao, n, nao, 1.0_dp, R, nao, Cb, nao, 0.0_dp, T, nao)
          call dgemm('T', 'N', n, n, nao, 1.0_dp, Cb, nao, T, nao, 0.0_dp, Rb, n)
          call eigh(n, Rb, w, ierr)
          call dgemm('N', 'N', nao, n, n, 1.0_dp, Cb, nao, Rb, n, 0.0_dp, T, nao)
          C(:,idx(1:n)) = T
          e(idx(1:n)) = w
          deallocate(Cb, T, Rb, w)
       enddo
    enddo
    deallocate(idx)

  end subroutine canonicalise

!######################################################################
! order_orbitals: PySCF's final ordering
!######################################################################
  subroutine order_orbitals(nao, C, e, occ, irr)

    integer(is), intent(in)    :: nao
    real(dp), intent(inout)    :: C(nao,nao), e(nao), occ(nao)
    integer(is), intent(inout) :: irr(nao)

    real(dp), allocatable    :: C0(:,:), e0(:), o0(:), key(:)
    integer(is), allocatable :: i0(:), perm(:), idx(:), p2(:)
    integer(is) :: i, n, m, cls

    allocate(C0(nao,nao), e0(nao), o0(nao), i0(nao), perm(nao), key(nao), idx(nao), p2(nao))
    C0 = C; e0 = e; o0 = occ; i0 = irr
    if (nirr <= 1) then
       call stable_sort(nao, e0, perm)
    else
       m = 0
       do cls = 2, 0, -1
          n = 0
          do i = 1, nao
             if (abs(o0(i) - real(cls,dp)) < 1.0e-6_dp) then
                n = n + 1
                idx(n) = i
                key(n) = anint(e0(i)*1.0e9_dp)
             endif
          enddo
          call stable_sort(n, key(1:n), p2(1:n))
          do i = 1, n
             perm(m+i) = idx(p2(i))
          enddo
          m = m + n
       enddo
    endif
    do i = 1, nao
       C(:,i) = C0(:,perm(i))
       e(i) = e0(perm(i))
       occ(i) = o0(perm(i))
       irr(i) = i0(perm(i))
    enddo
    deallocate(C0, e0, o0, i0, perm, key, idx, p2)

  end subroutine order_orbitals

end module scf_solver
