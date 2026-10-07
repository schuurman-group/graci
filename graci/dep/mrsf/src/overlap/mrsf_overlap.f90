!**********************************************************************
! mrsf_overlap: overlaps <Psi_I|Psi'_J> between MRSF-TDDFT states of two
! geometries (bra, ket), by determinant factorisation (Lee, Kim, Lee,
! Choi, JCTC 15, 882 (2019)) with the two-index determinants evaluated
! either exactly (adjugate / Schur-complement / Jacobi identities of the
! MO overlap matrix) or by the truncated Leibniz formula TLF(0/1/2).
!
! Determinants are taken in the canonical form (alpha MOs ascending,
! then beta MOs ascending). With R+ = |C Cbar O1 O2|, R- = |C Cbar O1bar
! O2bar|, Phi+_{ia} = a+_{a beta} a_{i alpha} R+, Phi-_{ia} = a+_{a alpha}
! a_{i beta} R- (reduced to canonical form with the fermionic phases
! c+_{ia} = (-1)^(posH(i) + nC + 1 + qC(a)), c-_{ia} = -c+_{ia}, where
! posH(i) is the rank of hole i among the hole MOs and qC(a) the number
! of core MOs below particle a), the state of an expanded amplitude
! matrix X~ (particle a, hole i) is
!   |Psi> = sum_{non-OO} X~_{ai} (Phi+_{ia} + s' Phi-_{ia})/sqrt2
!         + sum_{OO slots} X~_{ai} Phi+_{ia}
! with s' = -1 (singlets), +1 (triplets); verified against brute-force
! Loewdin overlaps, S^2 and the library's densities (Phase-0 prototype
! mrsf_overlap_ref.py). With the plain minors of M = Cb^T Sao Ck,
!   Sa_{ij} = det M[H\i, H'\j],  Sb_{ab} = det M[C+a, C'+b],
!   T_{ib}  = det M[H\i, C'+b],  U_{aj}  = det M[C+a, H'\j]
! (rows/columns in ascending MO order) and the phase-absorbed families
! A+_{ai} = c+ X~/sqrt2 (OO: c+ X~), A-_{ai} = s' c- X~/sqrt2 (OO: 0):
!   <Psi|Psi'> = sum (A+ A'+ + A- A'-) Sa Sb + sum (A+ A'- + A- A'+) T U,
! evaluated as Z = A Sa, Y = Sb A', W = A T, Q = U A'^T and dgemm over
! all state pairs.
!**********************************************************************
module mrsf_overlap

  use mrsf_constants
  use mrsf_io
  use mrsf_space, only: classify_orbitals
  use mrsf_density, only: expand_one

  implicit none

  real(dp)    :: time_ovl = 0.0_dp
  integer(is) :: novl_calls = 0

contains

!######################################################################
! state_overlap: entry point. occ_b/occ_k: reference occupations (2/1/0)
! of the bra/ket MOs; method: 0 exact, 1/2/3 = TLF(0/1/2); xb(xdim,nb),
! xk(xdim,nk): compressed amplitudes (particle fastest); S(nb,nk) out.
! ierr: 0 ok, 1 class mismatch, 2 singular core/hole block (explicit
! minors used), 3 unknown method
!######################################################################
  subroutine state_overlap(nao_b, nao_k, nmo, occ_b, occ_k, mult, method, nb, nk, &
       Cb, Ck, Sao, xb, xk, S, ierr)

    integer(is), intent(in)  :: nao_b, nao_k, nmo, mult, method, nb, nk
    real(dp), intent(in)     :: occ_b(nmo), occ_k(nmo), Cb(nao_b,nmo), Ck(nao_k,nmo)
    real(dp), intent(in)     :: Sao(nao_b,nao_k), xb(*), xk(*)
    real(dp), intent(out)    :: S(nb,nk)
    integer(is), intent(out) :: ierr

    integer(is), allocatable :: hb(:), pb(:), hk(:), pk(:)
    integer(is) :: ncb, nvb, nck, nvk, io1b, io2b, io1k, io2k, nocca, nvirb
    real(dp) :: t0

    t0 = wall_time()
    ierr = 0
    S = 0.0_dp
    call classify_orbitals(nmo, occ_b, ncb, nvb, hb, pb, io1b, io2b)
    call classify_orbitals(nmo, occ_k, nck, nvk, hk, pk, io1k, io2k)
    if (ncb /= nck .or. nvb /= nvk) then
       ierr = 1
       return
    endif
    if (method < 0 .or. method > 3) then
       ierr = 3
       return
    endif
    nocca = ncb + 2
    nvirb = nvb + 2
    call overlap_core(nao_b, nao_k, nmo, ncb, nocca, nvirb, hb, pb, hk, pk, mult, method, &
         nb, nk, Cb, Ck, Sao, xb, xk, S, ierr)
    deallocate(hb, pb, hk, pk)
    time_ovl = time_ovl + wall_time() - t0
    novl_calls = novl_calls + 1

  end subroutine state_overlap

!######################################################################
! overlap_core
!######################################################################
  subroutine overlap_core(nao_b, nao_k, nmo, nc, nocca, nvirb, hb, pb, hk, pk, mult, method, &
       nb, nk, Cb, Ck, Sao, xb, xk, S, ierr)

    integer(is), intent(in)    :: nao_b, nao_k, nmo, nc, nocca, nvirb, mult, method, nb, nk
    integer(is), intent(in)    :: hb(nocca), pb(nvirb), hk(nocca), pk(nvirb)
    real(dp), intent(in)       :: Cb(nao_b,nmo), Ck(nao_k,nmo), Sao(nao_b,nao_k)
    real(dp), intent(in)       :: xb(nvirb*nocca,nb), xk(nvirb*nocca,nk)
    real(dp), intent(out)      :: S(nb,nk)
    integer(is), intent(inout) :: ierr

    real(dp), allocatable :: M(:,:), tmp(:,:), Sa(:,:), Sb(:,:), T(:,:), U(:,:)
    real(dp), allocatable :: Ap(:,:,:), Am(:,:,:), Bp(:,:,:), Bm(:,:,:)
    real(dp), allocatable :: Zp(:,:,:), Zm(:,:,:), Yp(:,:,:), Ym(:,:,:)
    real(dp), allocatable :: Wp(:,:,:), Wm(:,:,:), Qp(:,:,:), Qm(:,:,:)
    integer(is) :: v, nz, nw

    ! MO overlap matrix M = Cb^T Sao Ck (bra rows, ket columns)
    allocate(M(nmo,nmo), tmp(nao_b,nmo))
    call dgemm('N','N', nao_b, nmo, nao_k, 1.0_dp, Sao, nao_b, Ck, nao_k, 0.0_dp, tmp, nao_b)
    call dgemm('T','N', nmo, nmo, nao_b, 1.0_dp, Cb, nao_b, tmp, nao_b, 0.0_dp, M, nmo)
    deallocate(tmp)

    ! two-index blocks
    allocate(Sa(nocca,nocca), Sb(nvirb,nvirb), T(nocca,nvirb), U(nvirb,nocca))
    if (method == 0) then
       call blocks_exact(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, Sa, Sb, T, U, ierr)
    else
       call blocks_tlf(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, method-1, Sa, Sb, T, U)
    endif

    ! phase-absorbed amplitude families
    allocate(Ap(nvirb,nocca,nb), Am(nvirb,nocca,nb), Bp(nvirb,nocca,nk), Bm(nvirb,nocca,nk))
    do v = 1, nb
       call families(nmo, nc, nocca, nvirb, hb, pb, mult, xb(1,v), Ap(1,1,v), Am(1,1,v))
    enddo
    do v = 1, nk
       call families(nmo, nc, nocca, nvirb, hk, pk, mult, xk(1,v), Bp(1,1,v), Bm(1,1,v))
    enddo

    ! direct term: Z = A Sa (nvirb x nocca), Y = Sb B (nvirb x nocca), S += <Z, Y>
    nz = nvirb*nocca
    allocate(Zp(nvirb,nocca,nb), Zm(nvirb,nocca,nb), Yp(nvirb,nocca,nk), Ym(nvirb,nocca,nk))
    do v = 1, nb
       call dgemm('N','N', nvirb, nocca, nocca, 1.0_dp, Ap(1,1,v), nvirb, Sa, nocca, 0.0_dp, Zp(1,1,v), nvirb)
       call dgemm('N','N', nvirb, nocca, nocca, 1.0_dp, Am(1,1,v), nvirb, Sa, nocca, 0.0_dp, Zm(1,1,v), nvirb)
    enddo
    do v = 1, nk
       call dgemm('N','N', nvirb, nocca, nvirb, 1.0_dp, Sb, nvirb, Bp(1,1,v), nvirb, 0.0_dp, Yp(1,1,v), nvirb)
       call dgemm('N','N', nvirb, nocca, nvirb, 1.0_dp, Sb, nvirb, Bm(1,1,v), nvirb, 0.0_dp, Ym(1,1,v), nvirb)
    enddo
    call dgemm('T','N', nb, nk, nz, 1.0_dp, Zp, nz, Yp, nz, 0.0_dp, S, nb)
    call dgemm('T','N', nb, nk, nz, 1.0_dp, Zm, nz, Ym, nz, 1.0_dp, S, nb)
    deallocate(Zp, Zm, Yp, Ym)

    ! cross term: W = A T (nvirb x nvirb), Q = U B^T (nvirb x nvirb), S += <W+, Q-> + <W-, Q+>
    nw = nvirb*nvirb
    allocate(Wp(nvirb,nvirb,nb), Wm(nvirb,nvirb,nb), Qp(nvirb,nvirb,nk), Qm(nvirb,nvirb,nk))
    do v = 1, nb
       call dgemm('N','N', nvirb, nvirb, nocca, 1.0_dp, Ap(1,1,v), nvirb, T, nocca, 0.0_dp, Wp(1,1,v), nvirb)
       call dgemm('N','N', nvirb, nvirb, nocca, 1.0_dp, Am(1,1,v), nvirb, T, nocca, 0.0_dp, Wm(1,1,v), nvirb)
    enddo
    do v = 1, nk
       call dgemm('N','T', nvirb, nvirb, nocca, 1.0_dp, U, nvirb, Bp(1,1,v), nvirb, 0.0_dp, Qp(1,1,v), nvirb)
       call dgemm('N','T', nvirb, nvirb, nocca, 1.0_dp, U, nvirb, Bm(1,1,v), nvirb, 0.0_dp, Qm(1,1,v), nvirb)
    enddo
    call dgemm('T','N', nb, nk, nw, 1.0_dp, Wp, nw, Qm, nw, 1.0_dp, S, nb)
    call dgemm('T','N', nb, nk, nw, 1.0_dp, Wm, nw, Qp, nw, 1.0_dp, S, nb)
    deallocate(Wp, Wm, Qp, Qm, Ap, Am, Bp, Bm, Sa, Sb, T, U, M)

  end subroutine overlap_core

!######################################################################
! families: A+ and A- (nvirb x nocca) of one compressed vector
!######################################################################
  subroutine families(nmo, nc, nocca, nvirb, hmap, pmap, mult, x, Ap, Am)

    integer(is), intent(in) :: nmo, nc, nocca, nvirb, hmap(nocca), pmap(nvirb), mult
    real(dp), intent(in)    :: x(nvirb,nocca)
    real(dp), intent(out)   :: Ap(nvirb,nocca), Am(nvirb,nocca)

    real(dp), allocatable :: Xe(:,:)
    integer(is), allocatable :: posH(:), qC(:)
    integer(is) :: il, al
    real(dp) :: cp, spair

    allocate(Xe(nvirb,nocca), posH(nocca), qC(nvirb))
    call expand_one(nc, nvirb, nocca, mult, x, Xe)
    call positions(nc, nocca, nvirb, hmap, pmap, posH, qC)
    spair = -1.0_dp
    if (mult /= 1) spair = 1.0_dp
    do il = 1, nocca
       do al = 1, nvirb
          cp = real((-1)**(posH(il) + nc + 1 + qC(al)), dp)
          if (al <= 2 .and. il > nc) then
             Ap(al,il) = Xe(al,il)*cp
             Am(al,il) = 0.0_dp
          else
             Ap(al,il) = Xe(al,il)*cp*isqrt2
             Am(al,il) = -spair*Xe(al,il)*cp*isqrt2
          endif
       enddo
    enddo
    deallocate(Xe, posH, qC)

  end subroutine families

!######################################################################
! positions: posH(il) = rank (0-based) of hole il among the hole MOs,
! qC(al) = number of core MOs below particle al
!######################################################################
  subroutine positions(nc, nocca, nvirb, hmap, pmap, posH, qC)

    integer(is), intent(in)  :: nc, nocca, nvirb, hmap(nocca), pmap(nvirb)
    integer(is), intent(out) :: posH(nocca), qC(nvirb)
    integer(is) :: il, kl, al

    do il = 1, nocca
       posH(il) = 0
       do kl = 1, nocca
          if (hmap(kl) < hmap(il)) posH(il) = posH(il) + 1
       enddo
    enddo
    do al = 1, nvirb
       qC(al) = 0
       do kl = 1, nc
          if (hmap(kl) < pmap(al)) qC(al) = qC(al) + 1
       enddo
    enddo

  end subroutine positions

!######################################################################
! sort_list: ascending insertion sort of a short integer list
!######################################################################
  subroutine sort_list(n, a)
    integer(is), intent(in)    :: n
    integer(is), intent(inout) :: a(n)
    integer(is) :: i, j, t
    do i = 2, n
       t = a(i); j = i - 1
       do while (j >= 1)
          if (a(j) <= t) exit
          a(j+1) = a(j); j = j - 1
       enddo
       a(j+1) = t
    enddo
  end subroutine sort_list

!######################################################################
! det_sub: determinant of M[rows, cols] (n x n) by LU
!######################################################################
  function det_sub(nmo, M, n, rows, cols) result(d)
    integer(is), intent(in) :: nmo, n, rows(n), cols(n)
    real(dp), intent(in)    :: M(nmo,nmo)
    real(dp) :: d
    real(dp), allocatable :: A(:,:)
    integer(is), allocatable :: ipiv(:)
    integer(is) :: i, j, info
    if (n == 0) then
       d = 1.0_dp
       return
    endif
    allocate(A(n,n), ipiv(n))
    do j = 1, n
       do i = 1, n
          A(i,j) = M(rows(i), cols(j))
       enddo
    enddo
    call dgetrf(n, n, A, n, ipiv, info)
    d = 1.0_dp
    if (info > 0) then
       d = 0.0_dp
    else
       do i = 1, n
          d = d*A(i,i)
          if (ipiv(i) /= i) d = -d
       enddo
    endif
    deallocate(A, ipiv)
  end function det_sub

!######################################################################
! lu_inverse: determinant and inverse of A (n x n); info /= 0 on failure
!######################################################################
  subroutine lu_inverse(n, A, detA, Ainv, info)
    integer(is), intent(in)  :: n
    real(dp), intent(in)     :: A(n,n)
    real(dp), intent(out)    :: detA, Ainv(n,n)
    integer(is), intent(out) :: info
    integer(is), allocatable :: ipiv(:)
    real(dp), allocatable :: work(:)
    integer(is) :: i
    detA = 1.0_dp
    info = 0
    if (n == 0) return
    allocate(ipiv(n), work(n*n))
    Ainv = A
    call dgetrf(n, n, Ainv, n, ipiv, info)
    if (info == 0) then
       do i = 1, n
          detA = detA*Ainv(i,i)
          if (ipiv(i) /= i) detA = -detA
       enddo
       call dgetri(n, Ainv, n, ipiv, work, n*n, info)
    endif
    deallocate(ipiv, work)
  end subroutine lu_inverse

!######################################################################
! blocks_exact: Sa, Sb, T, U from the inverses of the hole block
! M[Hs,Hs'] and the core block M[Cs,Cs'] (sorted MO lists)
!######################################################################
  subroutine blocks_exact(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, Sa, Sb, T, U, ierr)

    integer(is), intent(in)    :: nmo, nc, nocca, nvirb, hb(nocca), pb(nvirb), hk(nocca), pk(nvirb)
    real(dp), intent(in)       :: M(nmo,nmo)
    real(dp), intent(out)      :: Sa(nocca,nocca), Sb(nvirb,nvirb), T(nocca,nvirb), U(nvirb,nocca)
    integer(is), intent(inout) :: ierr

    integer(is), allocatable :: Hsb(:), Hsk(:), Csb(:), Csk(:), posHb(:), posHk(:), qCb(:), qCk(:)
    integer(is), allocatable :: qCOb(:), qCOk(:)
    real(dp), allocatable :: A(:,:), Ainv(:,:), Cc(:,:), Cinv(:,:), X1(:,:), X2(:,:)
    real(dp), allocatable :: MPC(:,:), MCP(:,:), MPP(:,:), MHP(:,:), MPH(:,:), W(:,:), Wk(:,:), rowv(:), colv(:)
    real(dp), allocatable :: m2b(:,:), m2k(:,:), tmpT(:,:), tmpU(:,:)
    real(dp) :: detA, detC, sg, d2
    integer(is) :: i, j, il, jl, al, bl, kl, ll, info, po1b, po2b, po1k, po2k, x, oo, ri, rk, cj, cl, i1, i2
    integer(is) :: sgn_r, sgn_c
    real(dp), parameter :: sing_tol = 1.0e-12_dp

    allocate(Hsb(nocca), Hsk(nocca), Csb(nc), Csk(nc), posHb(nocca), posHk(nocca), qCb(nvirb), qCk(nvirb))
    allocate(qCOb(nocca), qCOk(nocca))
    Hsb = hb; Hsk = hk; call sort_list(nocca, Hsb); call sort_list(nocca, Hsk)
    Csb = hb(1:nc); Csk = hk(1:nc); call sort_list(nc, Csb); call sort_list(nc, Csk)
    call positions(nc, nocca, nvirb, hb, pb, posHb, qCb)
    call positions(nc, nocca, nvirb, hk, pk, posHk, qCk)
    ! rank of hole il among the core MOs (for the SOMO rows/columns of the bordered blocks)
    do il = 1, nocca
       qCOb(il) = 0; qCOk(il) = 0
       do kl = 1, nc
          if (hb(kl) < hb(il)) qCOb(il) = qCOb(il) + 1
          if (hk(kl) < hk(il)) qCOk(il) = qCOk(il) + 1
       enddo
    enddo
    po1b = posHb(nc+1); po2b = posHb(nc+2)
    po1k = posHk(nc+1); po2k = posHk(nc+2)

    ! hole block and core block with their inverses
    allocate(A(nocca,nocca), Ainv(nocca,nocca), Cc(nc,nc), Cinv(nc,nc))
    do j = 1, nocca
       do i = 1, nocca
          A(i,j) = M(Hsb(i), Hsk(j))
       enddo
    enddo
    do j = 1, nc
       do i = 1, nc
          Cc(i,j) = M(Csb(i), Csk(j))
       enddo
    enddo
    call lu_inverse(nocca, A, detA, Ainv, info)
    if (info == 0) call lu_inverse(nc, Cc, detC, Cinv, info)
    if (info /= 0 .or. abs(detA) < sing_tol .or. abs(detC) < sing_tol) then
       ierr = 2
       call blocks_minors(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, Sa, Sb, T, U)
       deallocate(Hsb, Hsk, Csb, Csk, posHb, posHk, qCb, qCk, qCOb, qCOk, A, Ainv, Cc, Cinv)
       return
    endif

    ! Sa: cofactors of the hole block
    do jl = 1, nocca
       cj = posHk(jl)
       do il = 1, nocca
          ri = posHb(il)
          Sa(il,jl) = real((-1)**(ri + cj), dp)*detA*Ainv(cj+1, ri+1)
       enddo
    enddo

    ! gathered blocks
    allocate(MPC(nvirb,nc), MCP(nc,nvirb), MPP(nvirb,nvirb), MHP(nocca,nvirb), MPH(nvirb,nocca))
    do kl = 1, nc
       do al = 1, nvirb
          MPC(al,kl) = M(pb(al), Csk(kl))
       enddo
       do bl = 1, nvirb
          MCP(kl,bl) = M(Csb(kl), pk(bl))
       enddo
    enddo
    do bl = 1, nvirb
       do al = 1, nvirb
          MPP(al,bl) = M(pb(al), pk(bl))
       enddo
       do il = 1, nocca
          MHP(il,bl) = M(hb(il), pk(bl))
       enddo
    enddo
    do jl = 1, nocca
       do al = 1, nvirb
          MPH(al,jl) = M(pb(al), hk(jl))
       enddo
    enddo

    ! Sb: bordered core block, X1 = Cinv MCP, Schur = MPP - MPC X1
    allocate(X1(nc,nvirb), X2(nvirb,nc))
    if (nc > 0) then
       call dgemm('N','N', nc, nvirb, nc, 1.0_dp, Cinv, nc, MCP, nc, 0.0_dp, X1, nc)
       call dgemm('N','N', nvirb, nc, nc, 1.0_dp, MPC, nvirb, Cinv, nc, 0.0_dp, X2, nvirb)
       call dgemm('N','N', nvirb, nvirb, nc, -1.0_dp, MPC, nvirb, X1, nc, 1.0_dp, MPP, nvirb)
    endif
    do bl = 1, nvirb
       sgn_c = (-1)**(nc - qCk(bl))
       do al = 1, nvirb
          sgn_r = (-1)**(nc - qCb(al))
          Sb(al,bl) = real(sgn_r*sgn_c, dp)*detC*MPP(al,bl)
       enddo
    enddo

    ! T, U: SOMO rows/columns by the bordered form
    allocate(rowv(nvirb), colv(nvirb))
    do x = 1, 2
       il = nc + x                       ! hole O_x (bra); the other SOMO is the row
       oo = nc + 3 - x
       sgn_r = (-1)**(nc - qCOb(oo))
       ! rowv(b) = M(O_oo, b) - M(O_oo, Csk) Cinv M(Csb, b) = M(O_oo,b) - M(O_oo,Csk) X1(:,b)
       do bl = 1, nvirb
          rowv(bl) = M(hb(oo), pk(bl))
          do kl = 1, nc
             rowv(bl) = rowv(bl) - M(hb(oo), Csk(kl))*X1(kl,bl)
          enddo
          sgn_c = (-1)**(nc - qCk(bl))
          T(il,bl) = real(sgn_r*sgn_c, dp)*detC*rowv(bl)
       enddo
       jl = nc + x                       ! hole O_x (ket); the other SOMO is the column
       oo = nc + 3 - x
       sgn_c = (-1)**(nc - qCOk(oo))
       do al = 1, nvirb
          colv(al) = M(pb(al), hk(oo))
          do kl = 1, nc
             colv(al) = colv(al) - X2(al,kl)*M(Csb(kl), hk(oo))
          enddo
          sgn_r = (-1)**(nc - qCb(al))
          U(al,jl) = real(sgn_r*sgn_c, dp)*detC*colv(al)
       enddo
    enddo

    ! T, U: core holes by Laplace expansion with the Jacobi second minors
    ! m2b(il,kl) = det A[rows\{il,kl}, cols\{O1',O2'}], m2k(jl,ll) = det A[rows\{O1,O2}, cols\{jl,ll}]
    allocate(m2b(nocca,nocca), m2k(nocca,nocca), W(nocca,nocca), Wk(nocca,nocca), tmpT(nocca,nvirb), tmpU(nvirb,nocca))
    W = 0.0_dp; Wk = 0.0_dp
    do il = 1, nocca
       ri = posHb(il)
       do kl = 1, nocca
          if (kl == il) cycle
          rk = posHb(kl)
          i1 = min(ri, rk); i2 = max(ri, rk)
          d2 = Ainv(po1k+1, i1+1)*Ainv(po2k+1, i2+1) - Ainv(po1k+1, i2+1)*Ainv(po2k+1, i1+1)
          if (po1k > po2k) d2 = -d2
          m2b(il,kl) = detA*real((-1)**(ri + rk + po1k + po2k), dp)*d2
          ! rank of k in sorted(H\i)
          rk = posHb(kl)
          if (posHb(kl) > posHb(il)) rk = rk - 1
          W(il,kl) = real((-1)**rk, dp)*m2b(il,kl)
       enddo
    enddo
    do jl = 1, nocca
       cj = posHk(jl)
       do ll = 1, nocca
          if (ll == jl) cycle
          cl = posHk(ll)
          i1 = min(cj, cl); i2 = max(cj, cl)
          ! det Ainv[rows {i1,i2} (ket positions), cols {po1b,po2b} (bra positions)]
          d2 = Ainv(i1+1, po1b+1)*Ainv(i2+1, po2b+1) - Ainv(i1+1, po2b+1)*Ainv(i2+1, po1b+1)
          if (po1b > po2b) d2 = -d2
          m2k(jl,ll) = detA*real((-1)**(cj + cl + po1b + po2b), dp)*d2
          cl = posHk(ll)
          if (posHk(ll) > posHk(jl)) cl = cl - 1
          Wk(jl,ll) = real((-1)**cl, dp)*m2k(jl,ll)
       enddo
    enddo
    call dgemm('N','N', nocca, nvirb, nocca, 1.0_dp, W, nocca, MHP, nocca, 0.0_dp, tmpT, nocca)
    call dgemm('N','T', nvirb, nocca, nocca, 1.0_dp, MPH, nvirb, Wk, nocca, 0.0_dp, tmpU, nvirb)
    do bl = 1, nvirb
       sgn_c = (-1)**qCk(bl)
       do il = 1, nc
          T(il,bl) = real(sgn_c, dp)*tmpT(il,bl)
       enddo
    enddo
    do jl = 1, nc
       do al = 1, nvirb
          sgn_r = (-1)**qCb(al)
          U(al,jl) = real(sgn_r, dp)*tmpU(al,jl)
       enddo
    enddo

    deallocate(Hsb, Hsk, Csb, Csk, posHb, posHk, qCb, qCk, qCOb, qCOk)
    deallocate(A, Ainv, Cc, Cinv, X1, X2, MPC, MCP, MPP, MHP, MPH, rowv, colv, m2b, m2k, W, Wk, tmpT, tmpU)

  end subroutine blocks_exact

!######################################################################
! blocks_minors: the same blocks by explicit LU determinants (fallback)
!######################################################################
  subroutine blocks_minors(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, Sa, Sb, T, U)

    integer(is), intent(in) :: nmo, nc, nocca, nvirb, hb(nocca), pb(nvirb), hk(nocca), pk(nvirb)
    real(dp), intent(in)    :: M(nmo,nmo)
    real(dp), intent(out)   :: Sa(nocca,nocca), Sb(nvirb,nvirb), T(nocca,nvirb), U(nvirb,nocca)

    integer(is), allocatable :: rows(:), cols(:)
    integer(is) :: il, jl, al, bl, n1

    n1 = nocca - 1
    allocate(rows(n1), cols(n1))
    do jl = 1, nocca
       call remove_one(nocca, hk, jl, cols)
       call sort_list(n1, cols)
       do il = 1, nocca
          call remove_one(nocca, hb, il, rows)
          call sort_list(n1, rows)
          Sa(il,jl) = det_sub(nmo, M, n1, rows, cols)
       enddo
       do al = 1, nvirb
          rows(1:nc) = hb(1:nc); rows(n1) = pb(al)
          call sort_list(n1, rows)
          U(al,jl) = det_sub(nmo, M, n1, rows, cols)
       enddo
    enddo
    do bl = 1, nvirb
       cols(1:nc) = hk(1:nc); cols(n1) = pk(bl)
       call sort_list(n1, cols)
       do al = 1, nvirb
          rows(1:nc) = hb(1:nc); rows(n1) = pb(al)
          call sort_list(n1, rows)
          Sb(al,bl) = det_sub(nmo, M, n1, rows, cols)
       enddo
       do il = 1, nocca
          call remove_one(nocca, hb, il, rows)
          call sort_list(n1, rows)
          T(il,bl) = det_sub(nmo, M, n1, rows, cols)
       enddo
    enddo
    deallocate(rows, cols)

  end subroutine blocks_minors

  subroutine remove_one(n, a, k, b)
    integer(is), intent(in)  :: n, a(n), k
    integer(is), intent(out) :: b(n-1)
    integer(is) :: i, j
    j = 0
    do i = 1, n
       if (i == k) cycle
       j = j + 1
       b(j) = a(i)
    enddo
  end subroutine remove_one

!######################################################################
! blocks_tlf: truncated Leibniz expansions of the same minors, keeping
! the permutations with at most k off-diagonal factors
!######################################################################
  subroutine blocks_tlf(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, k, Sa, Sb, T, U)

    integer(is), intent(in) :: nmo, nc, nocca, nvirb, hb(nocca), pb(nvirb), hk(nocca), pk(nvirb), k
    real(dp), intent(in)    :: M(nmo,nmo)
    real(dp), intent(out)   :: Sa(nocca,nocca), Sb(nvirb,nvirb), T(nocca,nvirb), U(nvirb,nocca)

    integer(is), allocatable :: rows(:), cols(:)
    integer(is) :: il, jl, al, bl, n1

    n1 = nocca - 1
    allocate(rows(n1), cols(n1))
    !$omp parallel private(rows, cols, il, jl, al, bl)
    !$omp do
    do jl = 1, nocca
       call remove_one(nocca, hk, jl, cols)
       call sort_list(n1, cols)
       do il = 1, nocca
          call remove_one(nocca, hb, il, rows)
          call sort_list(n1, rows)
          Sa(il,jl) = tlf_minor(nmo, M, n1, rows, cols, k)
       enddo
       do al = 1, nvirb
          rows(1:nc) = hb(1:nc); rows(n1) = pb(al)
          call sort_list(n1, rows)
          U(al,jl) = tlf_minor(nmo, M, n1, rows, cols, k)
       enddo
    enddo
    !$omp end do
    !$omp do
    do bl = 1, nvirb
       cols(1:nc) = hk(1:nc); cols(n1) = pk(bl)
       call sort_list(n1, cols)
       do al = 1, nvirb
          rows(1:nc) = hb(1:nc); rows(n1) = pb(al)
          call sort_list(n1, rows)
          Sb(al,bl) = tlf_minor(nmo, M, n1, rows, cols, k)
       enddo
       do il = 1, nocca
          call remove_one(nocca, hb, il, rows)
          call sort_list(n1, rows)
          T(il,bl) = tlf_minor(nmo, M, n1, rows, cols, k)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    deallocate(rows, cols)

  end subroutine blocks_tlf

!######################################################################
! tlf_minor: det M[rows, cols] (sorted MO lists, length n) truncated to
! permutations with at most k off-diagonal factors (bra MO /= ket MO)
!######################################################################
  function tlf_minor(nmo, M, n, rows, cols, k) result(total)

    integer(is), intent(in) :: nmo, n, rows(n), cols(n), k
    real(dp), intent(in)    :: M(nmo,nmo)
    real(dp) :: total

    integer(is) :: common(n), extra_r(2), extra_c(2), ncom, mext, mr, mc, i, j, a, b, kk, p
    integer(is) :: map(n)         ! map(i) = ket MO assigned to row i
    logical :: found

    total = 0.0_dp
    if (n == 0) then
       total = 1.0_dp
       return
    endif
    ! common elements and the extra rows/columns
    ncom = 0; mr = 0; mc = 0
    do i = 1, n
       found = .false.
       do j = 1, n
          if (cols(j) == rows(i)) then
             found = .true.; exit
          endif
       enddo
       if (found) then
          ncom = ncom + 1; common(ncom) = rows(i)
       else
          mr = mr + 1
          if (mr > 2) return
          extra_r(mr) = rows(i)
       endif
    enddo
    do j = 1, n
       found = .false.
       do i = 1, n
          if (rows(i) == cols(j)) then
             found = .true.; exit
          endif
       enddo
       if (.not. found) then
          mc = mc + 1
          if (mc > 2) return
          extra_c(mc) = cols(j)
       endif
    enddo
    mext = mr
    if (mext > k) return
    ! base mapping: common elements to themselves
    do i = 1, n
       map(i) = rows(i)
    enddo
    if (mext == 0) then
       total = total + term(map)
       if (k >= 2) then
          do a = 1, ncom
             do b = a + 1, ncom
                call swap_map(map, common(a), common(b))
                total = total + term(map)
                call swap_map(map, common(a), common(b))
             enddo
          enddo
       endif
    else if (mext == 1) then
       call set_map(map, extra_r(1), extra_c(1))
       total = total + term(map)
       if (k >= 2) then
          do kk = 1, ncom
             call set_map(map, extra_r(1), common(kk))
             call set_map(map, common(kk), extra_c(1))
             total = total + term(map)
             call set_map(map, common(kk), common(kk))
          enddo
          call set_map(map, extra_r(1), extra_c(1))
       endif
    else
       call set_map(map, extra_r(1), extra_c(1)); call set_map(map, extra_r(2), extra_c(2))
       total = total + term(map)
       call set_map(map, extra_r(1), extra_c(2)); call set_map(map, extra_r(2), extra_c(1))
       total = total + term(map)
    endif

  contains

    subroutine set_map(mp, r, c)
      integer(is), intent(inout) :: mp(n)
      integer(is), intent(in)    :: r, c
      integer(is) :: q
      do q = 1, n
         if (rows(q) == r) then
            mp(q) = c; return
         endif
      enddo
    end subroutine set_map

    subroutine swap_map(mp, r1, r2)
      integer(is), intent(inout) :: mp(n)
      integer(is), intent(in)    :: r1, r2
      integer(is) :: q1, q2, t
      q1 = 0; q2 = 0
      do p = 1, n
         if (rows(p) == r1) q1 = p
         if (rows(p) == r2) q2 = p
      enddo
      t = mp(q1); mp(q1) = mp(q2); mp(q2) = t
    end subroutine swap_map

    function term(mp) result(val)
      integer(is), intent(in) :: mp(n)
      real(dp) :: val
      integer(is) :: perm(n), q, r, s, ninv
      ! column positions of the assigned ket MOs
      do q = 1, n
         perm(q) = 0
         do r = 1, n
            if (cols(r) == mp(q)) then
               perm(q) = r; exit
            endif
         enddo
      enddo
      ninv = 0
      do q = 1, n
         do r = q + 1, n
            if (perm(q) > perm(r)) ninv = ninv + 1
         enddo
      enddo
      val = 1.0_dp
      if (mod(ninv, 2) == 1) val = -1.0_dp
      do s = 1, n
         val = val*M(rows(s), mp(s))
      enddo
    end function term

  end function tlf_minor

end module mrsf_overlap

!**********************************************************************
! C-bound entry point
!**********************************************************************
module mrsf_overlap_interface

  use iso_c_binding
  use mrsf_constants
  use mrsf_overlap

  implicit none

contains

  subroutine mrsf_state_overlap_c(nao_b, nao_k, nmo, occ_b, occ_k, mult, method, nb, nk, &
       Cb, Ck, Sao, xb, xk, S, ierr) bind(c, name='mrsf_state_overlap')
    integer(is), intent(in)  :: nao_b, nao_k, nmo, mult, method, nb, nk
    real(dp), intent(in)     :: occ_b(*), occ_k(*), Cb(*), Ck(*), Sao(*), xb(*), xk(*)
    real(dp), intent(out)    :: S(*)
    integer(is), intent(out) :: ierr
    call state_overlap(nao_b, nao_k, nmo, occ_b, occ_k, mult, method, nb, nk, Cb, Ck, Sao, xb, xk, S, ierr)
  end subroutine mrsf_state_overlap_c

end module mrsf_overlap_interface
