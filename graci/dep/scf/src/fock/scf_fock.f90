!**********************************************************************
! scf_fock: Fock matrices, energies, the ROKS Roothaan matrix and the
! orbital gradients in PySCF's conventions:
!   ROKS: F_s = h + J[D_a + D_b] - c_x K[D_s] + V_xc^s,
!         E = E_nuc + Tr(h D) + 1/2 Tr(D J) - 1/2 c_x sum_s Tr(D_s K_s) + E_xc
!         Roothaan matrix (MO blocks): F_c on closed-closed, open-open,
!         virtual-virtual and closed-virtual, F_b on closed-open, F_a on
!         open-virtual (PySCF's projector formula); gradient elements
!         virtual-closed F_a + F_b, virtual-open F_a, open-closed F_b
!   RKS:  F = h + J[D] - 1/2 c_x K[D] + V_xc, gradient 2 F(virtual, occ)
! Shortcut (orbital densities): F_a - F_b is needed only in the rows and
! columns of the open orbitals, so the iterations form only F_c and
! (F_a - F_b) C_O.
!**********************************************************************
module scf_fock

  use scf_constants
  use scf_io
  use scf_df
  use scf_grid
  use scf_xc

  implicit none

  integer(is)           :: nao_f = 0
  real(dp), allocatable :: Hc(:,:), Sm(:,:)
  real(dp)              :: enuc_f = 0.0_dp, cx_f = 0.0_dp
  logical               :: have_xc = .false.
  real(dp)              :: ecomp(6) = 0.0_dp      ! E_tot, E_1e, E_coul, E_exch(HF), E_xc, E_nuc
  real(dp)              :: time_mat = 0.0_dp

  ! dgemm (scf_grid) and dsyrk (scf_df) are declared external there

contains

  subroutine fock_set(nao, h, S, enuc, cx, xc)
    integer(is), intent(in) :: nao
    real(dp), intent(in)    :: h(nao,nao), S(nao,nao), enuc, cx
    logical, intent(in)     :: xc
    if (allocated(Hc)) deallocate(Hc, Sm)
    allocate(Hc(nao,nao), Sm(nao,nao))
    nao_f = nao
    Hc = h
    Sm = S
    enuc_f = enuc
    cx_f = cx
    have_xc = xc
  end subroutine fock_set

  subroutine fock_free()
    if (allocated(Hc)) deallocate(Hc, Sm)
    nao_f = 0
  end subroutine fock_free

!######################################################################
! density from a factor: D = f X(:,1:k) X(:,1:k)^T
!######################################################################
  subroutine density(nao, k, X, f, D)
    integer(is), intent(in) :: nao, k
    real(dp), intent(in)    :: X(nao,*), f
    real(dp), intent(out)   :: D(nao,nao)
    D = 0.0_dp
    if (k > 0) then
       call dsyrk('U', 'N', nao, k, f, X, nao, 0.0_dp, D, nao)
       call fill_lower(nao, D)
    endif
  end subroutine density

!######################################################################
! build_roks: densities D_a = X X^T, D_b = X(:,1:kb) X(:,1:kb)^T; mode 0
! (exact): F_a, F_b and F_c; mode 1 (shortcut): F_c and dCO = (F_a -
! F_b) X(:,kb+1:ka); returns the energy
!######################################################################
  subroutine build_roks(nao, ka, kb, X, mode, E, Fa, Fb, Fc, dCO, Da, Db)

    integer(is), intent(in) :: nao, ka, kb, mode
    real(dp), intent(in)    :: X(nao,ka)
    real(dp), intent(out)   :: E, Fa(nao,nao), Fb(nao,nao), Fc(nao,nao), dCO(nao,*), Da(nao,nao), Db(nao,nao)

    real(dp), allocatable :: D(:,:), J(:,:), Kb_(:,:), Ko(:,:), V1(:,:), V2(:,:), dV(:,:)
    real(dp) :: exc, e1, ec, ex
    integer(is) :: no, wk

    no = ka - kb
    allocate(D(nao,nao), J(nao,nao), Kb_(nao,nao), Ko(nao,nao), V1(nao,nao), V2(nao,nao), dV(nao,max(no,1)))
    call density(nao, ka, X, 1.0_dp, Da)
    call density(nao, kb, X, 1.0_dp, Db)
    D = Da + Db
    wk = 0
    if (cx_f /= 0.0_dp) wk = 1
    call df_jk(nao, D, ka, kb, X, wk, J, Kb_, Ko)
    exc = 0.0_dp
    V1 = 0.0_dp; V2 = 0.0_dp; dV = 0.0_dp
    if (have_xc) call grid_pass_roks(nao, ka, kb, X, mode, exc, V1, V2, dV)
    e1 = sum(Hc*D)
    ec = 0.5_dp*sum(D*J)
    ex = -0.5_dp*cx_f*(sum(Da*(Kb_ + Ko)) + sum(Db*Kb_))
    E  = enuc_f + e1 + ec + ex + exc
    ecomp = [E, e1, ec, ex, exc, enuc_f]
    if (mode == 0) then
       Fa = Hc + J - cx_f*(Kb_ + Ko) + V1
       Fb = Hc + J - cx_f*Kb_ + V2
       Fc = 0.5_dp*(Fa + Fb)
       if (no > 0) dCO(1:nao,1:no) = matmul(Fa - Fb, X(:,kb+1:ka))
    else
       Fc = Hc + J - 0.5_dp*cx_f*(2.0_dp*Kb_ + Ko) + V1
       if (no > 0) then
          call dgemm('N', 'N', nao, no, nao, -cx_f, Ko, nao, X(1,kb+1), nao, 0.0_dp, dCO, nao)
          dCO(1:nao,1:no) = dCO(1:nao,1:no) + dV(:,1:no)
       endif
    endif
    deallocate(D, J, Kb_, Ko, V1, V2, dV)

  end subroutine build_roks

!######################################################################
! build_rks: D = 2 X X^T; F; returns the energy
!######################################################################
  subroutine build_rks(nao, k, X, E, F, D)

    integer(is), intent(in) :: nao, k
    real(dp), intent(in)    :: X(nao,k)
    real(dp), intent(out)   :: E, F(nao,nao), D(nao,nao)

    real(dp), allocatable :: J(:,:), K1(:,:), K2(:,:), V(:,:)
    real(dp) :: exc, e1, ec, ex
    integer(is) :: wk

    allocate(J(nao,nao), K1(nao,nao), K2(nao,nao), V(nao,nao))
    call density(nao, k, X, 2.0_dp, D)
    wk = 0
    if (cx_f /= 0.0_dp) wk = 2
    call df_jk(nao, D, k, k, X, wk, J, K1, K2)
    exc = 0.0_dp
    V = 0.0_dp
    if (have_xc) call grid_pass_rks(nao, k, X, exc, V)
    ! K[D] = 2 K1
    e1 = sum(Hc*D)
    ec = 0.5_dp*sum(D*J)
    ex = -0.5_dp*cx_f*sum(D*K1)
    E  = enuc_f + e1 + ec + ex + exc
    ecomp = [E, e1, ec, ex, exc, enuc_f]
    F = Hc + J - cx_f*K1 + V
    deallocate(J, K1, K2, V)

  end subroutine build_rks

!######################################################################
! roothaan_ao: PySCF's projector formula (any densities)
!######################################################################
  subroutine roothaan_ao(nao, Fa, Fb, Da, Db, R)

    integer(is), intent(in) :: nao
    real(dp), intent(in)    :: Fa(nao,nao), Fb(nao,nao), Da(nao,nao), Db(nao,nao)
    real(dp), intent(out)   :: R(nao,nao)

    real(dp), allocatable :: fc(:,:), pc(:,:), po(:,:), pv(:,:), T(:,:)
    integer(is) :: i
    real(dp) :: t0

    t0 = wall_time()
    allocate(fc(nao,nao), pc(nao,nao), po(nao,nao), pv(nao,nao), T(nao,nao))
    fc = 0.5_dp*(Fa + Fb)
    call dgemm('N', 'N', nao, nao, nao, 1.0_dp, Db, nao, Sm, nao, 0.0_dp, pc, nao)
    call dgemm('N', 'N', nao, nao, nao, 1.0_dp, Da - Db, nao, Sm, nao, 0.0_dp, po, nao)
    call dgemm('N', 'N', nao, nao, nao, -1.0_dp, Da, nao, Sm, nao, 0.0_dp, pv, nao)
    do i = 1, nao
       pv(i,i) = pv(i,i) + 1.0_dp
    enddo
    R = 0.0_dp
    call sandwich(pc, fc, pc, 0.5_dp)
    call sandwich(po, fc, po, 0.5_dp)
    call sandwich(pv, fc, pv, 0.5_dp)
    call sandwich(po, Fb, pc, 1.0_dp)
    call sandwich(po, Fa, pv, 1.0_dp)
    call sandwich(pv, fc, pc, 1.0_dp)
    T = transpose(R)
    R = R + T
    deallocate(fc, pc, po, pv, T)
    time_mat = time_mat + wall_time() - t0

  contains

    ! R += f A^T F B
    subroutine sandwich(A, Fm, Bm, f)
      real(dp), intent(in) :: A(nao,nao), Fm(nao,nao), Bm(nao,nao), f
      call dgemm('N', 'N', nao, nao, nao, 1.0_dp, Fm, nao, Bm, nao, 0.0_dp, T, nao)
      call dgemm('T', 'N', nao, nao, nao, f, A, nao, T, nao, 1.0_dp, R, nao)
    end subroutine sandwich

  end subroutine roothaan_ao

!######################################################################
! roothaan_short: the same matrix for orbital densities from F_c and
! dCO = (F_a - F_b) C_O; also returns fcmo = C^T F_c C and dmo = C^T dCO
! for the gradient
!######################################################################
  subroutine roothaan_short(nao, C, occ, Fc, dCO, R, fcmo, dmo)

    integer(is), intent(in) :: nao
    real(dp), intent(in)    :: C(nao,nao), occ(nao), Fc(nao,nao), dCO(nao,*)
    real(dp), intent(out)   :: R(nao,nao), fcmo(nao,nao), dmo(nao,*)

    real(dp), allocatable :: T(:,:), SC(:,:), RM(:,:)
    integer(is), allocatable :: op(:)
    integer(is) :: no, p, k
    real(dp) :: t0

    t0 = wall_time()
    no = count(abs(occ - 1.0_dp) < 1.0e-6_dp)
    allocate(T(nao,nao), SC(nao,nao), RM(nao,nao), op(max(no,1)))
    k = 0
    do p = 1, nao
       if (abs(occ(p) - 1.0_dp) < 1.0e-6_dp) then
          k = k + 1
          op(k) = p
       endif
    enddo
    call dgemm('N', 'N', nao, nao, nao, 1.0_dp, Fc, nao, C, nao, 0.0_dp, T, nao)
    call dgemm('T', 'N', nao, nao, nao, 1.0_dp, C, nao, T, nao, 0.0_dp, fcmo, nao)
    if (no > 0) call dgemm('T', 'N', nao, no, nao, 1.0_dp, C, nao, dCO, nao, 0.0_dp, dmo, nao)
    RM = fcmo
    do k = 1, no
       do p = 1, nao
          if (abs(occ(p) - 2.0_dp) < 1.0e-6_dp) then
             RM(p,op(k)) = RM(p,op(k)) - 0.5_dp*dmo(p,k)
             RM(op(k),p) = RM(op(k),p) - 0.5_dp*dmo(p,k)
          else if (abs(occ(p)) < 1.0e-6_dp) then
             RM(p,op(k)) = RM(p,op(k)) + 0.5_dp*dmo(p,k)
             RM(op(k),p) = RM(op(k),p) + 0.5_dp*dmo(p,k)
          endif
       enddo
    enddo
    call dgemm('N', 'N', nao, nao, nao, 1.0_dp, Sm, nao, C, nao, 0.0_dp, SC, nao)
    call dgemm('N', 'N', nao, nao, nao, 1.0_dp, SC, nao, RM, nao, 0.0_dp, T, nao)
    call dgemm('N', 'T', nao, nao, nao, 1.0_dp, T, nao, SC, nao, 0.0_dp, R, nao)
    deallocate(T, SC, RM, op)
    time_mat = time_mat + wall_time() - t0

  end subroutine roothaan_short

!######################################################################
! gradient norms (PySCF's get_grad; symmetry-forbidden elements zeroed)
!######################################################################
  function gnorm_roks_exact(nao, C, occ, irr, Fa, Fb) result(gn)

    integer(is), intent(in) :: nao, irr(nao)
    real(dp), intent(in)    :: C(nao,nao), occ(nao), Fa(nao,nao), Fb(nao,nao)
    real(dp)                :: gn
    real(dp), allocatable   :: T(:,:), fam(:,:), fbm(:,:)
    integer(is) :: p, q
    real(dp) :: g

    allocate(T(nao,nao), fam(nao,nao), fbm(nao,nao))
    call dgemm('N', 'N', nao, nao, nao, 1.0_dp, Fa, nao, C, nao, 0.0_dp, T, nao)
    call dgemm('T', 'N', nao, nao, nao, 1.0_dp, C, nao, T, nao, 0.0_dp, fam, nao)
    call dgemm('N', 'N', nao, nao, nao, 1.0_dp, Fb, nao, C, nao, 0.0_dp, T, nao)
    call dgemm('T', 'N', nao, nao, nao, 1.0_dp, C, nao, T, nao, 0.0_dp, fbm, nao)
    gn = 0.0_dp
    do q = 1, nao
       do p = 1, nao
          if (irr(p) /= irr(q)) cycle
          g = 0.0_dp
          if (occ(p) < 0.5_dp .and. occ(q) > 0.5_dp) g = g + fam(p,q)
          if (occ(p) < 1.5_dp .and. occ(q) > 1.5_dp) g = g + fbm(p,q)
          gn = gn + g*g
       enddo
    enddo
    gn = sqrt(gn)
    deallocate(T, fam, fbm)

  end function gnorm_roks_exact

  function gnorm_roks_short(nao, occ, irr, fcmo, dmo) result(gn)

    integer(is), intent(in) :: nao, irr(nao)
    real(dp), intent(in)    :: occ(nao), fcmo(nao,nao), dmo(nao,*)
    real(dp)                :: gn
    integer(is), allocatable :: oslot(:)
    integer(is) :: p, q, k
    real(dp) :: g

    allocate(oslot(nao))
    oslot = 0
    k = 0
    do p = 1, nao
       if (abs(occ(p) - 1.0_dp) < 1.0e-6_dp) then
          k = k + 1
          oslot(p) = k
       endif
    enddo
    gn = 0.0_dp
    do q = 1, nao
       do p = 1, nao
          if (irr(p) /= irr(q)) cycle
          g = 0.0_dp
          if (occ(p) < 0.5_dp .and. occ(q) > 1.5_dp) then
             g = 2.0_dp*fcmo(p,q)                              ! virtual-closed: F_a + F_b
          else if (occ(p) < 0.5_dp .and. oslot(q) > 0) then
             g = fcmo(p,q) + 0.5_dp*dmo(p,oslot(q))            ! virtual-open: F_a
          else if (oslot(p) > 0 .and. occ(q) > 1.5_dp) then
             g = fcmo(p,q) - 0.5_dp*dmo(q,oslot(p))            ! open-closed: F_b
          endif
          gn = gn + g*g
       enddo
    enddo
    gn = sqrt(gn)
    deallocate(oslot)

  end function gnorm_roks_short

  function gnorm_rks(nao, C, occ, irr, F) result(gn)

    integer(is), intent(in) :: nao, irr(nao)
    real(dp), intent(in)    :: C(nao,nao), occ(nao), F(nao,nao)
    real(dp)                :: gn
    real(dp), allocatable   :: T(:,:), fm(:,:)
    integer(is) :: p, q

    allocate(T(nao,nao), fm(nao,nao))
    call dgemm('N', 'N', nao, nao, nao, 1.0_dp, F, nao, C, nao, 0.0_dp, T, nao)
    call dgemm('T', 'N', nao, nao, nao, 1.0_dp, C, nao, T, nao, 0.0_dp, fm, nao)
    gn = 0.0_dp
    do q = 1, nao
       do p = 1, nao
          if (irr(p) /= irr(q)) cycle
          if (occ(p) < 0.5_dp .and. occ(q) > 0.5_dp) gn = gn + 4.0_dp*fm(p,q)**2
       enddo
    enddo
    gn = sqrt(gn)
    deallocate(T, fm)

  end function gnorm_rks

end module scf_fock
