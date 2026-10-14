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
! all state pairs. Extended (EMRSF) states add the CV configuration
! state functions on G = |C Cbar O1 O1bar|; their blocks (CV-CV,
! MRSF-CV) are evaluated by the same kind of identities on the G block
! of M (see extended_terms below).
!**********************************************************************
module mrsf_overlap

  use mrsf_constants
  use mrsf_io
  use mrsf_space, only: classify_orbitals
  use mrsf_density_mod, only: expand_one

  implicit none

  real(dp)    :: time_ovl = 0.0_dp
  integer(is) :: novl_calls = 0

contains

!######################################################################
! state_overlap: entry point. occ_b/occ_k: reference occupations (2/1/0)
! of the bra/ket MOs; method: 0 exact, 1/2/3 = TLF(0/1/2); align: 1 =
! align the ket MOs to the bra MOs within each orbital class (C, O, V)
! by orthogonal Procrustes rotations of the MO overlap blocks (the ket
! amplitudes are transformed accordingly, so the exact overlaps are
! unchanged and the TLF expansions become valid when the orbitals are
! reordered, sign-flipped or mixed between the two calculations; for
! extended objects the SOMO pair is only sign-aligned, see align_ket_ext);
! ncolb/ncolk: columns of the bra/ket vectors, nocca (standard MRSF; 0
! is accepted) or nocca + nC (extended MRSF: CV amplitudes in the last
! nC columns, SOMO rows masked); xb(nvirb*ncolb,nb), xk(nvirb*ncolk,nk):
! compressed amplitudes (particle fastest); S(nb,nk) out. ierr: 0 ok,
! 1 class mismatch, 2 singular core/hole/G block (explicit minors used),
! 3 unknown method, 4 SVD failure, 5 SOMO pair rotated or swapped
! between the geometries (extended objects, alignment restricted to
! sign flips; exact overlaps unaffected), 6 TLF closed forms of the CV
! blocks need class-ordered MOs (core, O1, O2, virtuals in index order),
! 7 = 2 and 5, 8 invalid column count
!######################################################################
  subroutine state_overlap(nao_b, nao_k, nmo, occ_b, occ_k, mult, method, align, nb, nk, ncolb, ncolk, &
       Cb, Ck, Sao, xb, xk, S, ierr)

    integer(is), intent(in)  :: nao_b, nao_k, nmo, mult, method, align, nb, nk, ncolb, ncolk
    real(dp), intent(in)     :: occ_b(nmo), occ_k(nmo), Cb(nao_b,nmo), Ck(nao_k,nmo)
    real(dp), intent(in)     :: Sao(nao_b,nao_k), xb(*), xk(*)
    real(dp), intent(out)    :: S(nb,nk)
    integer(is), intent(out) :: ierr

    integer(is), allocatable :: hb(:), pb(:), hk(:), pk(:)
    integer(is) :: ncb, nvb, nck, nvk, io1b, io2b, io1k, io2k, nocca, nvirb, nclb, nclk
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
    nclb = ncolb
    nclk = ncolk
    if (nclb == 0) nclb = nocca
    if (nclk == 0) nclk = nocca
    if ((nclb /= nocca .and. nclb /= nocca + ncb) .or. (nclk /= nocca .and. nclk /= nocca + ncb)) then
       ierr = 8
       return
    endif
    call overlap_core(nao_b, nao_k, nmo, ncb, nocca, nvirb, hb, pb, hk, pk, mult, method, align, &
         nb, nk, nclb, nclk, Cb, Ck, Sao, xb, xk, S, ierr)
    deallocate(hb, pb, hk, pk)
    time_ovl = time_ovl + wall_time() - t0
    novl_calls = novl_calls + 1

  end subroutine state_overlap

!######################################################################
! overlap_core
!######################################################################
  subroutine overlap_core(nao_b, nao_k, nmo, nc, nocca, nvirb, hb, pb, hk, pk, mult, method, align, &
       nb, nk, ncolb, ncolk, Cb, Ck, Sao, xb, xk, S, ierr)

    integer(is), intent(in)    :: nao_b, nao_k, nmo, nc, nocca, nvirb, mult, method, align, nb, nk, ncolb, ncolk
    integer(is), intent(in)    :: hb(nocca), pb(nvirb), hk(nocca), pk(nvirb)
    real(dp), intent(in)       :: Cb(nao_b,nmo), Ck(nao_k,nmo), Sao(nao_b,nao_k)
    real(dp), intent(in)       :: xb(nvirb*ncolb,nb), xk(nvirb*ncolk,nk)
    real(dp), intent(out)      :: S(nb,nk)
    integer(is), intent(inout) :: ierr

    real(dp), allocatable :: M(:,:), tmp(:,:), Sa(:,:), Sb(:,:), T(:,:), U(:,:)
    real(dp), allocatable :: Ap(:,:,:), Am(:,:,:), Bp(:,:,:), Bm(:,:,:)
    real(dp), allocatable :: Zp(:,:,:), Zm(:,:,:), Yp(:,:,:), Ym(:,:,:)
    real(dp), allocatable :: Wp(:,:,:), Wm(:,:,:), Qp(:,:,:), Qm(:,:,:)
    real(dp), allocatable :: Uh(:,:), Up(:,:), Xe(:,:), Xt(:,:), Yb(:,:,:), Yk(:,:,:), Yt(:,:)
    real(dp) :: sgn_o, offdiag
    integer(is) :: v, nz, nw, info, nv, jl, i0
    logical :: extb, extk, ext

    extb = ncolb > nocca
    extk = ncolk > nocca
    ext = extb .or. extk
    nv = nvirb - 2
    if (ext .and. method /= 0) then
       if (.not. natural_order(nc, nocca, hb, hk)) then
          ierr = 6
          return
       endif
    endif

    ! MO overlap matrix M = Cb^T Sao Ck (bra rows, ket columns)
    allocate(M(nmo,nmo), tmp(nao_b,nmo))
    call dgemm('N','N', nao_b, nmo, nao_k, 1.0_dp, Sao, nao_b, Ck, nao_k, 0.0_dp, tmp, nao_b)
    call dgemm('T','N', nmo, nmo, nao_b, 1.0_dp, Cb, nao_b, tmp, nao_b, 0.0_dp, M, nmo)
    deallocate(tmp)

    ! alignment of the ket MOs to the bra MOs within the classes: M -> M U,
    ! ket amplitudes X~ -> det(U_O) U_P^T X~ U_H (CV amplitudes Y -> U_V^T Y U_C)
    allocate(Uh(nocca,nocca), Up(nvirb,nvirb))
    sgn_o = 1.0_dp
    offdiag = 0.0_dp
    if (align == 1) then
       if (ext) then
          call align_ket_ext(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, Uh, Up, sgn_o, offdiag, info)
       else
          call align_ket(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, Uh, Up, sgn_o, info)
       endif
       if (info /= 0) then
          ierr = 4
          deallocate(M, Uh, Up)
          return
       endif
       if (ext .and. offdiag > 0.2_dp) ierr = 5
    endif

    ! two-index blocks
    allocate(Sa(nocca,nocca), Sb(nvirb,nvirb), T(nocca,nvirb), U(nvirb,nocca))
    if (method == 0) then
       call blocks_exact(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, Sa, Sb, T, U, ierr)
    else
       call blocks_tlf(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, method-1, Sa, Sb, T, U)
    endif

    ! phase-absorbed amplitude families (the first nocca columns of each vector)
    allocate(Ap(nvirb,nocca,nb), Am(nvirb,nocca,nb), Bp(nvirb,nocca,nk), Bm(nvirb,nocca,nk))
    do v = 1, nb
       call families(nmo, nc, nocca, nvirb, hb, pb, mult, xb(1,v), Ap(1,1,v), Am(1,1,v))
    enddo
    allocate(Xe(nvirb,nocca), Xt(nvirb,nocca))
    do v = 1, nk
       call expand_one(nc, nvirb, nocca, mult, xk(1,v), Xe)
       if (align == 1) then
          call dgemm('T','N', nvirb, nocca, nvirb, sgn_o, Up, nvirb, Xe, nvirb, 0.0_dp, Xt, nvirb)
          call dgemm('N','N', nvirb, nocca, nocca, 1.0_dp, Xt, nvirb, Uh, nocca, 0.0_dp, Xe, nvirb)
       endif
       call families_expanded(nc, nocca, nvirb, hk, pk, mult, Xe, Bp(1,1,v), Bm(1,1,v))
    enddo
    deallocate(Xe, Xt)

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
    deallocate(Wp, Wm, Qp, Qm)

    ! extended objects: CV amplitudes (virtual rows of the last nc columns)
    if (ext .and. nc > 0 .and. nv > 0) then
       allocate(Yb(nv,nc,nb), Yk(nv,nc,nk))
       Yb = 0.0_dp
       Yk = 0.0_dp
       if (extb) then
          do v = 1, nb
             do jl = 1, nc
                i0 = nvirb*(nocca + jl - 1)
                Yb(:,jl,v) = xb(i0+3:i0+nvirb, v)
             enddo
          enddo
       endif
       if (extk) then
          allocate(Yt(nv,nc))
          do v = 1, nk
             do jl = 1, nc
                i0 = nvirb*(nocca + jl - 1)
                Yk(:,jl,v) = xk(i0+3:i0+nvirb, v)
             enddo
             if (align == 1) then
                call dgemm('T','N', nv, nc, nv, 1.0_dp, Up(3,3), nvirb, Yk(1,1,v), nv, 0.0_dp, Yt, nv)
                call dgemm('N','N', nv, nc, nc, 1.0_dp, Yt, nv, Uh, nocca, 0.0_dp, Yk(1,1,v), nv)
             endif
          enddo
          deallocate(Yt)
       endif
       call extended_terms(nmo, nc, nocca, nvirb, hb, pb, hk, pk, mult, method, nb, nk, extb, extk, &
            M, Sa, Sb, T, U, Ap, Am, Bp, Bm, Yb, Yk, S, ierr)
       deallocate(Yb, Yk)
    endif
    deallocate(Ap, Am, Bp, Bm, Sa, Sb, T, U, M, Uh, Up)

  end subroutine overlap_core

!######################################################################
! families: A+ and A- (nvirb x nocca) of one compressed vector
!######################################################################
  subroutine families(nmo, nc, nocca, nvirb, hmap, pmap, mult, x, Ap, Am)

    integer(is), intent(in) :: nmo, nc, nocca, nvirb, hmap(nocca), pmap(nvirb), mult
    real(dp), intent(in)    :: x(nvirb,nocca)
    real(dp), intent(out)   :: Ap(nvirb,nocca), Am(nvirb,nocca)

    real(dp), allocatable :: Xe(:,:)

    allocate(Xe(nvirb,nocca))
    call expand_one(nc, nvirb, nocca, mult, x, Xe)
    call families_expanded(nc, nocca, nvirb, hmap, pmap, mult, Xe, Ap, Am)
    deallocate(Xe)

  end subroutine families

  subroutine families_expanded(nc, nocca, nvirb, hmap, pmap, mult, Xe, Ap, Am)

    integer(is), intent(in) :: nc, nocca, nvirb, hmap(nocca), pmap(nvirb), mult
    real(dp), intent(in)    :: Xe(nvirb,nocca)
    real(dp), intent(out)   :: Ap(nvirb,nocca), Am(nvirb,nocca)

    integer(is), allocatable :: posH(:), qC(:)
    integer(is) :: il, al
    real(dp) :: cp, spair

    allocate(posH(nocca), qC(nvirb))
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
    deallocate(posH, qC)

  end subroutine families_expanded

!######################################################################
! align_ket: orthogonal Procrustes alignment of the ket MOs to the bra
! MOs within the core, SOMO and virtual classes. For each class block
! B = M[bra X, ket X] = V S W^T, U_X = W V^T makes B U_X symmetric
! positive definite; M(:, ket X) <- M(:, ket X) U_X. Uh/Up: the hole /
! particle rotations in the local orderings; sgn_o = det(U_O) (global
! sign of the rotated reference determinants)
!######################################################################
  subroutine align_ket(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, Uh, Up, sgn_o, info)

    integer(is), intent(in)    :: nmo, nc, nocca, nvirb, hb(nocca), pb(nvirb), hk(nocca), pk(nvirb)
    real(dp), intent(inout)    :: M(nmo,nmo)
    real(dp), intent(out)      :: Uh(nocca,nocca), Up(nvirb,nvirb), sgn_o
    integer(is), intent(out)   :: info

    real(dp), allocatable :: Uc(:,:), Uo(:,:), Uv(:,:)
    integer(is) :: nv

    nv = nvirb - 2
    allocate(Uc(nc,nc), Uo(2,2), Uv(nv,nv))
    info = 0
    call procrustes(nmo, nc, hb(1:nc), hk(1:nc), M, Uc, info)
    if (info == 0) call procrustes(nmo, 2_is, hb(nc+1:nc+2), hk(nc+1:nc+2), M, Uo, info)
    if (info == 0) call procrustes(nmo, nv, pb(3:nvirb), pk(3:nvirb), M, Uv, info)
    Uh = 0.0_dp; Up = 0.0_dp
    Uh(1:nc,1:nc) = Uc
    Uh(nc+1:nc+2,nc+1:nc+2) = Uo
    Up(1:2,1:2) = Uo
    Up(3:nvirb,3:nvirb) = Uv
    sgn_o = Uo(1,1)*Uo(2,2) - Uo(1,2)*Uo(2,1)
    sgn_o = sign(1.0_dp, sgn_o)
    deallocate(Uc, Uo, Uv)

  end subroutine align_ket

  subroutine procrustes(nmo, n, rows, cols, M, U, info)

    integer(is), intent(in)    :: nmo, n, rows(n), cols(n)
    real(dp), intent(inout)    :: M(nmo,nmo)
    real(dp), intent(out)      :: U(n,n)
    integer(is), intent(out)   :: info

    real(dp), allocatable :: B(:,:), Usv(:,:), VT(:,:), sv(:), work(:), tmp(:,:), tmp2(:,:)
    integer(is) :: i, j, lwork

    info = 0
    if (n == 0) return
    allocate(B(n,n), Usv(n,n), VT(n,n), sv(n), tmp(nmo,n), tmp2(nmo,n))
    do j = 1, n
       do i = 1, n
          B(i,j) = M(rows(i), cols(j))
       enddo
    enddo
    lwork = max(1, 5*n)*n + 64*n
    allocate(work(lwork))
    call dgesvd('A', 'A', n, n, B, n, sv, Usv, n, VT, n, work, lwork, info)
    if (info /= 0) then
       deallocate(B, Usv, VT, sv, work, tmp, tmp2)
       return
    endif
    ! U = W V^T = VT^T Usv^T
    call dgemm('T','T', n, n, n, 1.0_dp, VT, n, Usv, n, 0.0_dp, U, n)
    ! M(:, cols) <- M(:, cols) U
    do j = 1, n
       tmp(:,j) = M(:,cols(j))
    enddo
    call dgemm('N','N', nmo, n, n, 1.0_dp, tmp, nmo, U, n, 0.0_dp, tmp2, nmo)
    do j = 1, n
       M(:,cols(j)) = tmp2(:,j)
    enddo
    deallocate(B, Usv, VT, sv, work, tmp, tmp2)

  end subroutine procrustes

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

!######################################################################
! ---------------- extended (EMRSF) states: CV contributions -----------
! An extended state adds sum_{j in C, b in V} y_bj chi_jb to the MRSF
! expansion, chi_jb = (D^a_jb + t D^b_jb)/sqrt2 with the canonical
! determinants D^a_jb = (sorted((G\j)+b); G), D^b_jb = (G; sorted((G\j)+b)),
! G = |C Cbar O1 O1bar| (the G configuration of the extended method),
! t = +1 (singlets) / -1 (triplets), and the coefficients
!   ca(j,b) = sG (-1)^(pG(j) + nb(j,b)) / sqrt2,  cb = t ca,
! sG = (-1)^(posH(O2) + nC + 1 + qC(O1)) (phase of the MRSF slot O2->O1
! which is the same determinant), pG(j) the rank of j in G, nb(j,b) the
! number of G\j below b. All overlaps are products of plain minors of M
! (CV-CV: det M[(G\j)+b, (G'\j')+b'] det M[G,G'] and the Cramer-type
! minors det M[(G\j)+b, G'] det M[G, (G'\j')+b']; MRSF-CV: the minors
! X1 = det M[H\i, (G'\j')+b'], X2 = det M[C+a, G'], X3 = det M[H\i, G'],
! X4 = det M[C+a, (G'\j')+b'] with <Phi+|D'^a> = X1 X2, <Phi+|D'^b> =
! X3 X4, <Phi-|D'^a> = X4 X3, <Phi-|D'^b> = X2 X1), evaluated in closed
! form (one LU of N = M[G,G'] plus the hole/core LUs: Jacobi, Schur and
! Cramer identities) or by the truncated Leibniz formula (closed forms
! of the same separable structure; class-ordered MOs required), and
! contracted with the amplitudes by dgemm over all state pairs. The
! bra-CV/ket-MRSF terms use the same cross routine on M^T. Validated
! against brute-force Loewdin overlaps (emrsf_overlap_ref.py).
! Implementation note: every table that is combined elementwise with
! the CV amplitudes Y(b,j) is stored in the same (virtual, core) layout
! and the contractions over Y go through BLAS; loops that read a
! (core, virtual) table against Y(b,j,v) were miscompiled by ifx 2025.3
! with -ipo (Y addressed with the stride of the table), see the
! synthetic reproducer tests/dbg_synth_tlf.py of the development tree.
!######################################################################

!######################################################################
! extended_terms: adds the CV-CV, MRSF-CV and CV-MRSF contributions to S
!######################################################################
  subroutine extended_terms(nmo, nc, nocca, nvirb, hb, pb, hk, pk, mult, method, nb, nk, extb, extk, &
       M, Sa, Sb, T, U, Ap, Am, Bp, Bm, Yb, Yk, S, ierr)

    integer(is), intent(in)    :: nmo, nc, nocca, nvirb, hb(nocca), pb(nvirb), hk(nocca), pk(nvirb)
    integer(is), intent(in)    :: mult, method, nb, nk
    logical, intent(in)        :: extb, extk
    real(dp), intent(in)       :: M(nmo,nmo), Sa(nocca,nocca), Sb(nvirb,nvirb), T(nocca,nvirb), U(nvirb,nocca)
    real(dp), intent(in)       :: Ap(nvirb,nocca,nb), Am(nvirb,nocca,nb), Bp(nvirb,nocca,nk), Bm(nvirb,nocca,nk)
    real(dp), intent(in)       :: Yb(nvirb-2,nc,nb), Yk(nvirb-2,nc,nk)
    real(dp), intent(inout)    :: S(nb,nk)
    integer(is), intent(inout) :: ierr

    real(dp), allocatable :: MT(:,:), R(:,:), X2(:), X3(:)
    integer(is) :: kt
    real(dp) :: sGb, sGk, tcsf

    kt = method - 1
    tcsf = 1.0_dp
    if (mult /= 1) tcsf = -1.0_dp
    sGb = g_phase(nc, nocca, nvirb, hb, pb)
    sGk = g_phase(nc, nocca, nvirb, hk, pk)

    if (extb .and. extk) then
       call cv_contract(nmo, nc, nocca, nvirb, hb, pb, hk, pk, nb, nk, M, Yb, Yk, sGb, sGk, tcsf, method, kt, S, ierr)
    endif
    ! bra MRSF, ket CV
    if (extk) then
       allocate(R(nb,nk), X2(nvirb), X3(nocca))
       X2 = Sb(:,1)
       X3 = Sa(:,nc+2)
       call cross_contract(M, nmo, nc, nocca, nvirb, hb, pb, hk, pk, nb, nk, Ap, Am, X2, X3, Yk, sGk, tcsf, &
            method, kt, R, ierr)
       S = S + R
       deallocate(R, X2, X3)
    endif
    ! bra CV, ket MRSF: the same routine on M^T with the roles exchanged
    if (extb) then
       allocate(MT(nmo,nmo), R(nk,nb), X2(nvirb), X3(nocca))
       MT = transpose(M)
       X2 = T(nc+2,:)
       X3 = U(1,:)
       call cross_contract(MT, nmo, nc, nocca, nvirb, hk, pk, hb, pb, nk, nb, Bp, Bm, X2, X3, Yb, sGb, tcsf, &
            method, kt, R, ierr)
       S = S + transpose(R)
       deallocate(MT, R, X2, X3)
    endif

  end subroutine extended_terms

!######################################################################
! g_phase: sG = (-1)^(posH(O2) + nc + 1 + qC(O1)), the phase of the
! MRSF slot (O1 <- O2) = sG |G> in the canonical form
!######################################################################
  function g_phase(nc, nocca, nvirb, hmap, pmap) result(sG)

    integer(is), intent(in) :: nc, nocca, nvirb, hmap(nocca), pmap(nvirb)
    real(dp) :: sG
    integer(is), allocatable :: posH(:), qC(:)

    allocate(posH(nocca), qC(nvirb))
    call positions(nc, nocca, nvirb, hmap, pmap, posH, qC)
    sG = real((-1)**(posH(nc+2) + nc + 1 + qC(1)), dp)
    deallocate(posH, qC)

  end function g_phase

!######################################################################
! natural_order: the orbital classes are in index order (core 1..nc,
! O1 = nc+1, O2 = nc+2, virtuals above) on both sides, as required by
! the TLF closed forms of the CV blocks
!######################################################################
  function natural_order(nc, nocca, hb, hk) result(ok)

    integer(is), intent(in) :: nc, nocca, hb(nocca), hk(nocca)
    logical :: ok

    ok = (hb(nc+1) == nc + 1 .and. hb(nc+2) == nc + 2 .and. hk(nc+1) == nc + 1 .and. hk(nc+2) == nc + 2)
    if (nc > 0) ok = ok .and. (hb(nc) == nc .and. hk(nc) == nc)

  end function natural_order

!######################################################################
! perm_sign_list: sign of the permutation sorting a list of distinct
! integers
!######################################################################
  function perm_sign_list(n, a) result(s)

    integer(is), intent(in) :: n, a(n)
    integer(is) :: s, i, j

    s = 1
    do i = 1, n
       do j = i + 1, n
          if (a(i) > a(j)) s = -s
       enddo
    enddo

  end function perm_sign_list

!######################################################################
! cv_tables: per-side data of the CV amplitudes of one side: the sorted
! G list, the rank pG(jl) (1-based) of core jl in it, the coefficient
! table ca(bl,jl) and the sort signs sr(bl,jl) = sign((G\j)+b -> sorted)
! (all tables in the (virtual, core) layout of the CV amplitudes: every
! elementwise combination with Y keeps one memory layout)
!######################################################################
  subroutine cv_tables(nmo, nc, nocca, nvirb, hmap, pmap, sG, Gs, pG, ca, sr)

    integer(is), intent(in)  :: nmo, nc, nocca, nvirb, hmap(nocca), pmap(nvirb)
    real(dp), intent(in)     :: sG
    integer(is), intent(out) :: Gs(nc+1), pG(nc)
    real(dp), intent(out)    :: ca(nvirb-2,nc), sr(nvirb-2,nc)

    integer(is) :: nG, nv, jl, bl, i, nbelow, j, b

    nG = nc + 1
    nv = nvirb - 2
    Gs(1:nc) = hmap(1:nc)
    Gs(nG) = hmap(nc+1)
    call sort_list(nG, Gs)
    do jl = 1, nc
       pG(jl) = 0
       do i = 1, nG
          if (Gs(i) == hmap(jl)) pG(jl) = i
       enddo
    enddo
    do jl = 1, nc
       j = hmap(jl)
       do bl = 1, nv
          b = pmap(2+bl)
          nbelow = 0
          do i = 1, nG
             if (Gs(i) /= j .and. Gs(i) < b) nbelow = nbelow + 1
          enddo
          ca(bl,jl) = sG*real((-1)**(pG(jl) - 1 + nbelow), dp)*isqrt2
          sr(bl,jl) = real((-1)**(nG - 1 - nbelow), dp)
       enddo
    enddo

  end subroutine cv_tables

!######################################################################
! cv_contract: CV-CV contributions. With the signed amplitudes the
! overlap of two CV parts is
!   sGb sGk { dN [ sum_i <S_i^T Y K_i, Y'> + (Y.Qx)(Px.Y') ] + t dR (Y.Qr)(Pr.Y') }
! exact: dN = dR = detN^2, K = N^-1 (core block, transposed indexing),
! S = Schur complement of N in M_VV', Qx = Qr = M_VG' N^-1, Px = Pr =
! (N^-1 M_GV')^T; TLF(k): the truncated-Leibniz closed forms
! (cv_tlf_pieces). Px, Pr, Qx, Qr and the amplitudes share the (virtual,
! core) layout.
!######################################################################
  subroutine cv_contract(nmo, nc, nocca, nvirb, hb, pb, hk, pk, nb, nk, M, Yb, Yk, sGb, sGk, t, method, kt, S, ierr)

    integer(is), intent(in)    :: nmo, nc, nocca, nvirb, hb(nocca), pb(nvirb), hk(nocca), pk(nvirb), nb, nk
    integer(is), intent(in)    :: method, kt
    real(dp), intent(in)       :: M(nmo,nmo), Yb(nvirb-2,nc,nb), Yk(nvirb-2,nc,nk), sGb, sGk, t
    real(dp), intent(inout)    :: S(nb,nk)
    integer(is), intent(inout) :: ierr

    integer(is), allocatable :: Gsb(:), Gsk(:), pGb(:), pGk(:)
    real(dp), allocatable :: cab(:,:), srb(:,:), cak(:,:), srk(:,:)
    real(dp), allocatable :: N(:,:), Ninv(:,:), MVG(:,:), MGV(:,:), Sch(:,:), Q(:,:), P(:,:)
    real(dp), allocatable :: K(:,:,:), Sv(:,:,:), Qx(:,:), Px(:,:), Qr(:,:), Pr(:,:)
    real(dp), allocatable :: Ub(:,:,:), tmp(:,:), qxs(:), qrs(:), pxs(:), prs(:), R(:,:), dSv(:,:)
    logical :: diagS(3)
    integer(is) :: nv, nG, i, j, b, jl, jl2, info, nkp, ip, v, bl
    real(dp) :: detN, dN, dR, fac
    real(dp), parameter :: sing_tol = 1.0e-12_dp

    nv = nvirb - 2
    nG = nc + 1
    if (nc == 0 .or. nv == 0) return
    allocate(Gsb(nG), Gsk(nG), pGb(nc), pGk(nc), cab(nv,nc), srb(nv,nc), cak(nv,nc), srk(nv,nc))
    call cv_tables(nmo, nc, nocca, nvirb, hb, pb, sGb, Gsb, pGb, cab, srb)
    call cv_tables(nmo, nc, nocca, nvirb, hk, pk, sGk, Gsk, pGk, cak, srk)

    allocate(K(nc,nc,3), Sv(nv,nv,3), Qx(nv,nc), Px(nv,nc), Qr(nv,nc), Pr(nv,nc))
    K = 0.0_dp; Sv = 0.0_dp; Qx = 0.0_dp; Px = 0.0_dp; Qr = 0.0_dp; Pr = 0.0_dp
    diagS = .false.
    nkp = 0
    if (method == 0) then
       allocate(N(nG,nG), Ninv(nG,nG))
       do j = 1, nG
          do i = 1, nG
             N(i,j) = M(Gsb(i), Gsk(j))
          enddo
       enddo
       call lu_inverse(nG, N, detN, Ninv, info)
       if (info /= 0 .or. abs(detN) < sing_tol) then
          if (ierr == 5) then
             ierr = 7
          else if (ierr /= 7) then
             ierr = 2
          endif
          call cv_fallback(nmo, nc, nocca, nvirb, nv, nG, hb, pb, hk, pk, nb, nk, M, Yb, Yk, Gsb, Gsk, pGb, pGk, &
               cab, srb, cak, srk, t, S)
          deallocate(N, Ninv, K, Sv, Qx, Px, Qr, Pr, Gsb, Gsk, pGb, pGk, cab, srb, cak, srk)
          return
       endif
       allocate(MVG(nv,nG), MGV(nG,nv), Sch(nv,nv), Q(nv,nG), P(nG,nv))
       do i = 1, nG
          do b = 1, nv
             MVG(b,i) = M(pb(2+b), Gsk(i))
             MGV(i,b) = M(Gsb(i), pk(2+b))
          enddo
       enddo
       do bl = 1, nv
          do b = 1, nv
             Sch(b,bl) = M(pb(2+b), pk(2+bl))
          enddo
       enddo
       call dgemm('N','N', nv, nG, nG, 1.0_dp, MVG, nv, Ninv, nG, 0.0_dp, Q, nv)
       call dgemm('N','N', nG, nv, nG, 1.0_dp, Ninv, nG, MGV, nG, 0.0_dp, P, nG)
       call dgemm('N','N', nv, nv, nG, -1.0_dp, Q, nv, MGV, nG, 1.0_dp, Sch, nv)
       nkp = 1
       do jl2 = 1, nc
          do jl = 1, nc
             K(jl,jl2,1) = Ninv(pGk(jl2), pGb(jl))
          enddo
       enddo
       Sv(:,:,1) = Sch
       do jl = 1, nc
          Qx(:,jl) = Q(:,pGb(jl))
          Px(:,jl) = P(pGk(jl),:)
       enddo
       Qr = Qx
       Pr = Px
       dN = detN*detN
       dR = dN
       deallocate(N, Ninv, MVG, MGV, Sch, Q, P)
    else
       call cv_tlf_pieces(nmo, nc, nv, nG, M, kt, nkp, K, Sv, diagS, Qx, Px, Qr, Pr, dN, dR)
    endif

    ! contraction over all state pairs
    allocate(Ub(nv,nc,nb), tmp(nv,nc), qxs(nb), qrs(nb), pxs(nk), prs(nk), R(nb,nk), dSv(nv,3))
    do ip = 1, nkp
       do b = 1, nv
          dSv(b,ip) = Sv(b,b,ip)
       enddo
    enddo
    Ub = 0.0_dp
    do v = 1, nb
       do ip = 1, nkp
          call dgemm('N','N', nv, nc, nc, 1.0_dp, Yb(1,1,v), nv, K(1,1,ip), nc, 0.0_dp, tmp, nv)
          if (diagS(ip)) then
             do jl = 1, nc
                Ub(:,jl,v) = Ub(:,jl,v) + dSv(:,ip)*tmp(:,jl)
             enddo
          else
             call dgemm('T','N', nv, nc, nv, 1.0_dp, Sv(1,1,ip), nv, tmp, nv, 1.0_dp, Ub(1,1,v), nv)
          endif
       enddo
       qxs(v) = sum(Yb(:,:,v)*Qx)
       qrs(v) = sum(Yb(:,:,v)*Qr)
    enddo
    do v = 1, nk
       pxs(v) = sum(Yk(:,:,v)*Px)
       prs(v) = sum(Yk(:,:,v)*Pr)
    enddo
    call dgemm('T','N', nb, nk, nv*nc, 1.0_dp, Ub, nv*nc, Yk, nv*nc, 0.0_dp, R, nb)
    fac = sGb*sGk
    do v = 1, nk
       do i = 1, nb
          S(i,v) = S(i,v) + fac*(dN*(R(i,v) + qxs(i)*pxs(v)) + t*dR*qrs(i)*prs(v))
       enddo
    enddo
    deallocate(Ub, tmp, qxs, qrs, pxs, prs, R, dSv, K, Sv, Qx, Px, Qr, Pr, Gsb, Gsk, pGb, pGk, cab, srb, cak, srk)

  end subroutine cv_contract

!######################################################################
! cv_tlf_pieces: TLF(k) closed forms of the CV-CV blocks (class-ordered
! MOs: G = 1..nG, virtuals nc+3..nmo). With d(p) = M(p,p), D_G = prod d,
! E = M_VG D^-1 M_GV, E_VG = M_VG D^-1 M_GG, E_GV = M_GG D^-1 M_VG,
! G2 = M_CG D^-1 M_GC, F = sum_{x<y} M(x,y)M(y,x)/(d(x)d(y)) over G:
!  Xaa/D_G = K1 (x) S1 + K2 (x) diag(d_V) + K3 (x) E + Qx (x) Px
! (the sign factors (-1)^(pG+pG') of the Jacobi identity absorbed), with
!  K1 = delta/d(j) - [k>=2] (1-delta) M(j',j)/(d(j)d(j')),   S1 = M_VV' ([k=0]: diag)
!  K2 = [k=1] -(1-delta) M(j',j)/(d(j)d(j'))
!       [k=2] -delta (F - F_j)/d(j) + (1-delta)(G2(j',j) - 2 M(j',j))/(d(j)d(j'))
!  K3 = -[k>=2] delta/d(j);  Qx = [k>=2] M(b,j)/d(j), Px = M(j',b')/d(j')
! and the Cramer-type minors Xr = (-1)^(nG-1-pG) D_G Qr, Xc likewise:
!  Qr = ([k>=1] M(b,j) - [k>=2] (E_VG(b,j) - M(b,j)))/d(j), Pr analogous;
! detN(TLF) = D_G (k<=1), D_G (1 - F) (k=2); dN = detN D_G, dR = D_G^2
!######################################################################
  subroutine cv_tlf_pieces(nmo, nc, nv, nG, M, kt, nkp, K, Sv, diagS, Qx, Px, Qr, Pr, dN, dR)

    integer(is), intent(in)  :: nmo, nc, nv, nG, kt
    real(dp), intent(in)     :: M(nmo,nmo)
    integer(is), intent(out) :: nkp
    real(dp), intent(out)    :: K(nc,nc,3), Sv(nv,nv,3), Qx(nv,nc), Px(nv,nc), Qr(nv,nc), Pr(nv,nc), dN, dR
    logical, intent(out)     :: diagS(3)

    real(dp), allocatable :: d(:), MVGd(:,:), MGV(:,:), MGGd(:,:), MGC(:,:), E(:,:), EVG(:,:), EGV(:,:), G2(:,:), Fj(:)
    real(dp) :: DG, F, detN
    integer(is) :: p, q, j, jl, b, bl, kk

    allocate(d(nG), Fj(nG))
    do p = 1, nG
       d(p) = M(p,p)
    enddo
    DG = product(d)
    F = 0.0_dp
    Fj = 0.0_dp
    do p = 1, nG
       do q = p + 1, nG
          F = F + M(p,q)*M(q,p)/(d(p)*d(q))
       enddo
       do q = 1, nG
          if (q /= p) Fj(p) = Fj(p) + M(q,p)*M(p,q)/(d(q)*d(p))
       enddo
    enddo
    ! scaled products over G
    allocate(MVGd(nv,nG), MGV(nG,nv), MGGd(nc,nG), MGC(nG,nc), E(nv,nv), EVG(nv,nc), EGV(nc,nv), G2(nc,nc))
    do kk = 1, nG
       do b = 1, nv
          MVGd(b,kk) = M(nc+2+b, kk)/d(kk)
          MGV(kk,b) = M(kk, nc+2+b)
       enddo
       do j = 1, nc
          MGGd(j,kk) = M(j,kk)/d(kk)
          MGC(kk,j) = M(kk,j)
       enddo
    enddo
    call dgemm('N','N', nv, nv, nG, 1.0_dp, MVGd, nv, MGV, nG, 0.0_dp, E, nv)
    call dgemm('N','N', nv, nc, nG, 1.0_dp, MVGd, nv, MGC, nG, 0.0_dp, EVG, nv)
    call dgemm('N','N', nc, nv, nG, 1.0_dp, MGGd, nc, MGV, nG, 0.0_dp, EGV, nc)
    call dgemm('N','N', nc, nc, nG, 1.0_dp, MGGd, nc, MGC, nG, 0.0_dp, G2, nc)

    K = 0.0_dp; Sv = 0.0_dp; Qx = 0.0_dp; Px = 0.0_dp; Qr = 0.0_dp; Pr = 0.0_dp
    diagS = .false.
    ! K1 (x) S1
    nkp = 1
    do jl = 1, nc
       K(jl,jl,1) = 1.0_dp/d(jl)
       if (kt >= 2) then
          do j = 1, nc
             if (j /= jl) K(jl,j,1) = -M(j,jl)/(d(jl)*d(j))
          enddo
       endif
    enddo
    if (kt >= 1) then
       do bl = 1, nv
          do b = 1, nv
             Sv(b,bl,1) = M(nc+2+b, nc+2+bl)
          enddo
       enddo
    else
       diagS(1) = .true.
       do b = 1, nv
          Sv(b,b,1) = M(nc+2+b, nc+2+b)
       enddo
    endif
    ! K2 (x) diag(d_V)
    if (kt >= 1) then
       nkp = 2
       diagS(2) = .true.
       do b = 1, nv
          Sv(b,b,2) = M(nc+2+b, nc+2+b)
       enddo
       do jl = 1, nc
          do j = 1, nc
             if (j == jl) then
                if (kt >= 2) K(jl,j,2) = -(F - Fj(jl))/d(jl)
             else
                if (kt == 1) then
                   K(jl,j,2) = -M(j,jl)/(d(jl)*d(j))
                else
                   K(jl,j,2) = (G2(j,jl) - 2.0_dp*M(j,jl))/(d(jl)*d(j))
                endif
             endif
          enddo
       enddo
    endif
    ! K3 (x) E and the Qx (x) Px piece
    if (kt >= 2) then
       nkp = 3
       Sv(:,:,3) = E
       do jl = 1, nc
          K(jl,jl,3) = -1.0_dp/d(jl)
          do b = 1, nv
             Qx(b,jl) = M(nc+2+b, jl)/d(jl)
             Px(b,jl) = M(jl, nc+2+b)/d(jl)
          enddo
       enddo
    endif
    ! Cramer-type minors
    if (kt >= 1) then
       do jl = 1, nc
          do b = 1, nv
             Qr(b,jl) = M(nc+2+b, jl)/d(jl)
             Pr(b,jl) = M(jl, nc+2+b)/d(jl)
             if (kt >= 2) then
                Qr(b,jl) = Qr(b,jl) - (EVG(b,jl) - M(nc+2+b, jl))/d(jl)
                Pr(b,jl) = Pr(b,jl) - (EGV(jl,b) - M(jl, nc+2+b))/d(jl)
             endif
          enddo
       enddo
    endif
    detN = DG
    if (kt >= 2) detN = DG*(1.0_dp - F)
    dN = detN*DG
    dR = DG*DG
    deallocate(d, Fj, MVGd, MGV, MGGd, MGC, E, EVG, EGV, G2)

  end subroutine cv_tlf_pieces

!######################################################################
! cv_fallback: CV-CV contributions with a singular N = M[G,G']: the
! bordered minors per (j,j') from the LU of N without row j / column j'
! (explicit determinants when that is singular too), the Cramer-type
! minors by explicit determinants
!######################################################################
  subroutine cv_fallback(nmo, nc, nocca, nvirb, nv, nG, hb, pb, hk, pk, nb, nk, M, Yb, Yk, Gsb, Gsk, pGb, pGk, &
       cab, srb, cak, srk, t, S)

    integer(is), intent(in) :: nmo, nc, nocca, nvirb, nv, nG, hb(nocca), pb(nvirb), hk(nocca), pk(nvirb), nb, nk
    integer(is), intent(in) :: Gsb(nG), Gsk(nG), pGb(nc), pGk(nc)
    real(dp), intent(in)    :: M(nmo,nmo), Yb(nv,nc,nb), Yk(nv,nc,nk), cab(nv,nc), srb(nv,nc), cak(nv,nc), srk(nv,nc), t
    real(dp), intent(inout) :: S(nb,nk)

    integer(is), allocatable :: rl(:), cl(:), rows(:), cols(:)
    real(dp), allocatable :: Ytb(:,:,:), Ytk(:,:,:), rI(:), cJ(:), Njj(:,:), Ninv(:,:), Rb(:,:), Cb(:,:)
    real(dp), allocatable :: X(:,:), tmp(:,:), W(:,:)
    real(dp) :: detN, dJ, xr, xc
    integer(is) :: jl, jl2, b, bl, v, i, j, info
    real(dp), parameter :: sing_tol = 1.0e-12_dp

    detN = det_sub(nmo, M, nG, Gsb, Gsk)
    allocate(Ytb(nv,nc,nb), Ytk(nv,nc,nk), rI(nb), cJ(nk), rows(nG), cols(nG))
    do v = 1, nb
       Ytb(:,:,v) = Yb(:,:,v)*cab
    enddo
    do v = 1, nk
       Ytk(:,:,v) = Yk(:,:,v)*cak
    enddo
    ! Cramer-type minors: Xr(j,b) = det M[(G\j)+b, G'], Xc(j',b') = det M[G, (G'\j')+b']
    rI = 0.0_dp; cJ = 0.0_dp
    do jl = 1, nc
       do b = 1, nv
          call remove_one(nG, Gsb, pGb(jl), rows(1:nG-1)); rows(nG) = pb(2+b); call sort_list(nG, rows)
          xr = det_sub(nmo, M, nG, rows, Gsk)
          rI = rI + xr*Ytb(b,jl,:)
          call remove_one(nG, Gsk, pGk(jl), cols(1:nG-1)); cols(nG) = pk(2+b); call sort_list(nG, cols)
          xc = det_sub(nmo, M, nG, Gsb, cols)
          cJ = cJ + xc*Ytk(b,jl,:)
       enddo
    enddo
    do v = 1, nk
       S(:,v) = S(:,v) + 2.0_dp*t*rI*cJ(v)
    enddo
    ! bordered minors per (j,j')
    allocate(rl(nG-1), cl(nG-1), Njj(nG-1,nG-1), Ninv(nG-1,nG-1), Rb(nv,nG-1), Cb(nG-1,nv), X(nv,nv), tmp(nv,nG-1), W(nv,nb))
    do jl2 = 1, nc
       call remove_one(nG, Gsk, pGk(jl2), cl)
       do jl = 1, nc
          call remove_one(nG, Gsb, pGb(jl), rl)
          do j = 1, nG - 1
             do i = 1, nG - 1
                Njj(i,j) = M(rl(i), cl(j))
             enddo
          enddo
          call lu_inverse(nG-1, Njj, dJ, Ninv, info)
          if (info == 0 .and. abs(dJ) >= sing_tol) then
             do i = 1, nG - 1
                do b = 1, nv
                   Rb(b,i) = M(pb(2+b), cl(i))
                   Cb(i,b) = M(rl(i), pk(2+b))
                enddo
             enddo
             do bl = 1, nv
                do b = 1, nv
                   X(b,bl) = M(pb(2+b), pk(2+bl))
                enddo
             enddo
             call dgemm('N','N', nv, nG-1, nG-1, 1.0_dp, Rb, nv, Ninv, nG-1, 0.0_dp, tmp(1,1), nv)
             call dgemm('N','N', nv, nv, nG-1, -1.0_dp, tmp(1,1), nv, Cb, nG-1, 1.0_dp, X, nv)
             do bl = 1, nv
                do b = 1, nv
                   X(b,bl) = srb(b,jl)*srk(bl,jl2)*dJ*X(b,bl)
                enddo
             enddo
          else
             do bl = 1, nv
                cols(1:nG-1) = cl; cols(nG) = pk(2+bl); call sort_list(nG, cols)
                do b = 1, nv
                   rows(1:nG-1) = rl; rows(nG) = pb(2+b); call sort_list(nG, rows)
                   X(b,bl) = det_sub(nmo, M, nG, rows, cols)
                enddo
             enddo
          endif
          ! S += 2 detN Ytb(:,jl,:)^T X Ytk(:,jl2,:)
          call dgemm('T','N', nv, nb, nv, 1.0_dp, X, nv, Ytb(1,jl,1), nv*nc, 0.0_dp, W, nv)
          call dgemm('T','N', nb, nk, nv, 2.0_dp*detN, W, nv, Ytk(1,jl2,1), nv*nc, 1.0_dp, S, nb)
       enddo
    enddo
    deallocate(Ytb, Ytk, rI, cJ, rows, cols, rl, cl, Njj, Ninv, Rb, Cb, X, tmp, W)

  end subroutine cv_fallback

!######################################################################
! cross_contract: contributions <MRSF side | CV side> for the rows of
! Mx (MRSF side: hole/particle maps hm/pm, families Ap/Am) against the
! columns (CV side: hc/pc, amplitudes Yc, G phase sGc):
!   R(I,J) = B_I . g1_J + A_I . g4_J,
!   B_I(i) = sum_a (Ap + t Am)_I(a,i) X2(a),  A_I(a) = sum_i (Am + t Ap)_I(a,i) X3(i),
!   g1_J(i) = sum_k' X1(i;k') ca(k') y'_J(k'),  g4_J(a) = sum_k' X4(a;k') ca(k') y'_J(k')
! (X2 = det Mx[C+a, G'], X3 = det Mx[H\i, G'] supplied by the caller from
! the standard blocks). Exact: X1 by Laplace expansion along the
! particle column with the Jacobi second minors of the hole block, X4
! by Cramer's rule on the bordered core block; TLF: cross_tlf
!######################################################################
  subroutine cross_contract(Mx, nmo, nc, nocca, nvirb, hm, pm, hc, pc, nstm, nstc, Ap, Am, X2, X3, Yc, sGc, t, &
       method, kt, R, ierr)

    integer(is), intent(in)    :: nmo, nc, nocca, nvirb, hm(nocca), pm(nvirb), hc(nocca), pc(nvirb), nstm, nstc
    integer(is), intent(in)    :: method, kt
    real(dp), intent(in)       :: Mx(nmo,nmo), Ap(nvirb,nocca,nstm), Am(nvirb,nocca,nstm), X2(nvirb), X3(nocca)
    real(dp), intent(in)       :: Yc(nvirb-2,nc,nstc), sGc, t
    real(dp), intent(out)      :: R(nstm,nstc)
    integer(is), intent(inout) :: ierr

    integer(is), allocatable :: Hsm(:), Hsc(:), Csm(:), Csc(:), Gsc(:), posHm(:), qCm(:), pHc(:), pCc(:), pGc(:), L(:)
    real(dp), allocatable :: ca(:,:), sc1(:,:), sc4(:,:), sr(:), g1(:,:), g4(:,:), B(:,:), Aa(:,:)
    real(dp), allocatable :: A(:,:), Ainv(:,:), Cc(:,:), Cinv(:,:), MCP(:,:), MPC(:,:), tS(:,:), X1c(:,:), MHV(:,:)
    real(dp), allocatable :: W1(:,:,:), Yt(:,:), Z(:,:), vJ(:), g1s(:)
    integer(is) :: nv, nG, i, j, il, jl, al, bl, kl, v, info, pO2c, pi, pk, rk, i0, i1, nbelow, pj
    real(dp) :: detA, detC, d2, aJ
    real(dp), external :: ddot
    real(dp), parameter :: sing_tol = 1.0e-12_dp

    nv = nvirb - 2
    nG = nc + 1
    R = 0.0_dp
    if (nc == 0 .or. nv == 0) return
    allocate(Hsm(nocca), Hsc(nocca), Csm(nc), Csc(nc), Gsc(nG), posHm(nocca), qCm(nvirb), pHc(nc), pCc(nc), pGc(nc))
    Hsm = hm; call sort_list(nocca, Hsm)
    Hsc = hc; call sort_list(nocca, Hsc)
    Csm = hm(1:nc); call sort_list(nc, Csm)
    Csc = hc(1:nc); call sort_list(nc, Csc)
    Gsc(1:nc) = hc(1:nc); Gsc(nG) = hc(nc+1); call sort_list(nG, Gsc)
    call positions(nc, nocca, nvirb, hm, pm, posHm, qCm)
    do jl = 1, nc
       pHc(jl) = 0; pCc(jl) = 0; pGc(jl) = 0
       do i = 1, nocca
          if (Hsc(i) == hc(jl)) pHc(jl) = i
       enddo
       do i = 1, nc
          if (Csc(i) == hc(jl)) pCc(jl) = i
       enddo
       do i = 1, nG
          if (Gsc(i) == hc(jl)) pGc(jl) = i
       enddo
    enddo
    pO2c = 0
    do i = 1, nocca
       if (Hsc(i) == hc(nc+2)) pO2c = i
    enddo

    ! coefficient and sign tables of the CV side
    allocate(ca(nv,nc), sc1(nv,nc), sc4(nv,nc), sr(nvirb), L(nc+1))
    do jl = 1, nc
       do bl = 1, nv
          nbelow = 0
          do i = 1, nG
             if (Gsc(i) /= hc(jl) .and. Gsc(i) < pc(2+bl)) nbelow = nbelow + 1
          enddo
          ca(bl,jl) = sGc*real((-1)**(pGc(jl) - 1 + nbelow), dp)*isqrt2
          sc1(bl,jl) = real((-1)**(nG - 1 - nbelow), dp)
          L(1:nc) = Csc; L(pCc(jl)) = pc(2+bl); L(nc+1) = hc(nc+1)
          sc4(bl,jl) = real(perm_sign_list(nc+1, L), dp)
       enddo
    enddo
    do al = 1, nvirb
       sr(al) = real((-1)**(nc - qCm(al)), dp)
    enddo

    allocate(g1(nocca,nstc), g4(nvirb,nstc), B(nocca,nstm), Aa(nvirb,nstm))
    g1 = 0.0_dp; g4 = 0.0_dp

    if (method == 0) then
       allocate(A(nocca,nocca), Ainv(nocca,nocca), Cc(nc,nc), Cinv(nc,nc))
       do j = 1, nocca
          do i = 1, nocca
             A(i,j) = Mx(Hsm(i), Hsc(j))
          enddo
       enddo
       do j = 1, nc
          do i = 1, nc
             Cc(i,j) = Mx(Csm(i), Csc(j))
          enddo
       enddo
       call lu_inverse(nocca, A, detA, Ainv, info)
       if (info == 0) call lu_inverse(nc, Cc, detC, Cinv, info)
       if (info /= 0 .or. abs(detA) < sing_tol .or. abs(detC) < sing_tol) then
          if (ierr == 5) then
             ierr = 7
          else if (ierr /= 7) then
             ierr = 2
          endif
          call cross_fallback(Mx, nmo, nc, nocca, nvirb, nv, nG, hm, pm, hc, pc, nstc, Hsm, Csm, Gsc, pGc, ca, Yc, g1, g4)
       else
          ! X4 by Cramer's rule on the bordered core block:
          ! X4(a;j',b') = sr(a) sc4(j',b') detC [tS(a,O1') X1c(j',b') - X1c(j',O1') tS(a,b')]
          allocate(MCP(nc,nvirb), MPC(nvirb,nc), tS(nvirb,nvirb), X1c(nc,nvirb), MHV(nocca,nv))
          do bl = 1, nvirb
             do kl = 1, nc
                MCP(kl,bl) = Mx(Csm(kl), pc(bl))
             enddo
             do al = 1, nvirb
                tS(al,bl) = Mx(pm(al), pc(bl))
             enddo
          enddo
          do kl = 1, nc
             do al = 1, nvirb
                MPC(al,kl) = Mx(pm(al), Csc(kl))
             enddo
          enddo
          call dgemm('N','N', nc, nvirb, nc, 1.0_dp, Cinv, nc, MCP, nc, 0.0_dp, X1c, nc)
          call dgemm('N','N', nvirb, nvirb, nc, -1.0_dp, MPC, nvirb, X1c, nc, 1.0_dp, tS, nvirb)
          allocate(Yt(nv,nc), vJ(nv), Z(nocca,nc), g1s(nocca))
          do v = 1, nstc
             do jl = 1, nc
                Yt(:,pCc(jl)) = sc4(:,jl)*ca(:,jl)*Yc(:,jl,v)
             enddo
             aJ = 0.0_dp
             do pj = 1, nc
                aJ = aJ + ddot(nv, X1c(pj,3), nc, Yt(1,pj), 1)
             enddo
             call dgemv('N', nv, nc, 1.0_dp, Yt, nv, X1c(1,1), 1, 0.0_dp, vJ, 1)
             do al = 1, nvirb
                g4(al,v) = tS(al,1)*aJ
             enddo
             call dgemv('N', nvirb, nv, -1.0_dp, tS(1,3), nvirb, vJ, 1, 1.0_dp, g4(1,v), 1)
             do al = 1, nvirb
                g4(al,v) = sr(al)*detC*g4(al,v)
             enddo
          enddo
          ! X1 by Laplace expansion along the particle column: X1(i;j',b') = sc1(j',b') sum_k W1(i,k;j') Mx(k,b'),
          ! W1 = (-1)^(rk + nG - 1) detA (-1)^(pi+pk+pj'+pO2') det Ainv[{pj',pO2'},{pi,pk}]
          allocate(W1(nocca,nocca,nc))
          W1 = 0.0_dp
          do jl = 1, nc
             pj = pHc(jl)
             do pi = 1, nocca
                do pk = 1, nocca
                   if (pk == pi) cycle
                   rk = pk - 1
                   if (pk > pi) rk = rk - 1
                   i0 = min(pi, pk); i1 = max(pi, pk)
                   d2 = Ainv(pj, i0)*Ainv(pO2c, i1) - Ainv(pj, i1)*Ainv(pO2c, i0)
                   if (pj > pO2c) d2 = -d2
                   W1(pi,pk,jl) = real((-1)**(rk + nG - 1 + pi + pk + pj + pO2c), dp)*detA*d2
                enddo
             enddo
          enddo
          do bl = 1, nv
             do pk = 1, nocca
                MHV(pk,bl) = Mx(Hsm(pk), pc(2+bl))
             enddo
          enddo
          do v = 1, nstc
             do jl = 1, nc
                Yt(:,jl) = sc1(:,jl)*ca(:,jl)*Yc(:,jl,v)
             enddo
             call dgemm('N','N', nocca, nc, nv, 1.0_dp, MHV, nocca, Yt, nv, 0.0_dp, Z, nocca)
             g1s = 0.0_dp
             do jl = 1, nc
                call dgemv('N', nocca, nocca, 1.0_dp, W1(1,1,jl), nocca, Z(1,jl), 1, 1.0_dp, g1s, 1)
             enddo
             do il = 1, nocca
                g1(il,v) = g1s(posHm(il)+1)
             enddo
          enddo
          deallocate(MCP, MPC, tS, X1c, MHV, Yt, vJ, Z, g1s, W1)
       endif
       deallocate(A, Ainv, Cc, Cinv)
    else
       call cross_tlf(Mx, nmo, nc, nocca, nvirb, nv, nG, nstc, Yc, sGc, kt, g1, g4)
    endif

    ! bra-side vectors and the contraction over all pairs
    do v = 1, nstm
       do il = 1, nocca
          B(il,v) = sum((Ap(:,il,v) + t*Am(:,il,v))*X2)
       enddo
       do al = 1, nvirb
          Aa(al,v) = ddot(nocca, Am(al,1,v), nvirb, X3, 1) + t*ddot(nocca, Ap(al,1,v), nvirb, X3, 1)
       enddo
    enddo
    call dgemm('T','N', nstm, nstc, nocca, 1.0_dp, B, nocca, g1, nocca, 0.0_dp, R, nstm)
    call dgemm('T','N', nstm, nstc, nvirb, 1.0_dp, Aa, nvirb, g4, nvirb, 1.0_dp, R, nstm)
    deallocate(Hsm, Hsc, Csm, Csc, Gsc, posHm, qCm, pHc, pCc, pGc, L, ca, sc1, sc4, sr, g1, g4, B, Aa)

  end subroutine cross_contract

!######################################################################
! cross_fallback: g1, g4 from explicit determinants (singular hole or
! core block)
!######################################################################
  subroutine cross_fallback(Mx, nmo, nc, nocca, nvirb, nv, nG, hm, pm, hc, pc, nstc, Hsm, Csm, Gsc, pGc, ca, Yc, g1, g4)

    integer(is), intent(in) :: nmo, nc, nocca, nvirb, nv, nG, hm(nocca), pm(nvirb), hc(nocca), pc(nvirb), nstc
    integer(is), intent(in) :: Hsm(nocca), Csm(nc), Gsc(nG), pGc(nc)
    real(dp), intent(in)    :: Mx(nmo,nmo), ca(nv,nc), Yc(nv,nc,nstc)
    real(dp), intent(out)   :: g1(nocca,nstc), g4(nvirb,nstc)

    integer(is), allocatable :: rows(:), cols(:)
    integer(is) :: jl, bl, il, al, v, pil
    real(dp) :: x1, x4

    allocate(rows(nG), cols(nG))
    g1 = 0.0_dp; g4 = 0.0_dp
    do jl = 1, nc
       do bl = 1, nv
          call remove_one(nG, Gsc, pGc(jl), cols(1:nG-1)); cols(nG) = pc(2+bl); call sort_list(nG, cols)
          do il = 1, nocca
             pil = 0
             do v = 1, nocca
                if (Hsm(v) == hm(il)) pil = v
             enddo
             call remove_one(nocca, Hsm, pil, rows)
             x1 = det_sub(nmo, Mx, nG, rows, cols)
             do v = 1, nstc
                g1(il,v) = g1(il,v) + x1*ca(bl,jl)*Yc(bl,jl,v)
             enddo
          enddo
          do al = 1, nvirb
             rows(1:nc) = Csm; rows(nG) = pm(al); call sort_list(nG, rows)
             x4 = det_sub(nmo, Mx, nG, rows, cols)
             do v = 1, nstc
                g4(al,v) = g4(al,v) + x4*ca(bl,jl)*Yc(bl,jl,v)
             enddo
          enddo
       enddo
    enddo
    deallocate(rows, cols)

  end subroutine cross_fallback

!######################################################################
! cross_tlf: g1, g4 from the TLF(k) closed forms of X1 and X4
! (class-ordered MOs: core 1..nc, O1 = nc+1, O2 = nc+2, virtuals above;
! kappa = sG (-1)^nc / sqrt2 carries the CSF coefficients):
!  X1(i;j',b'): i = O2: -(-1)^(j'+nc) D_G Pr(j',b')
!   i in G, i /= j' [k=2]: -(-1)^(i+j') D_G/(d(i)d(j')) [M(j',i)M(O2,b') - M(j',b')M(O2,i)]
!   i = j' in C: D_G/d(i) {[k>=1] M(O2,b') - [k=2](E_O2V(b') - M(O2,i)M(i,b')/d(i))}
!  X4(a;j',b'): a = O1: as X1 for i = O2
!   a in O2 + V: (-1)^(j'+nc) D_C/d(j') {[k=2][M(j',O1)M(a,b') - M(j',b')M(a,O1)]
!                 + delta_ab' d(a) ([k>=1] M(j',O1) - [k=2] G2c(j'))},
!   G2c(j') = sum_{c in C} M(j',c) M(c,O1)/d(c)
!######################################################################
  subroutine cross_tlf(Mx, nmo, nc, nocca, nvirb, nv, nG, nstc, Yc, sGc, kt, g1, g4)

    integer(is), intent(in) :: nmo, nc, nocca, nvirb, nv, nG, nstc, kt
    real(dp), intent(in)    :: Mx(nmo,nmo), Yc(nv,nc,nstc), sGc
    real(dp), intent(out)   :: g1(nocca,nstc), g4(nvirb,nstc)

    real(dp), allocatable :: d(:), EO2V(:), mO2V(:), MCVt(:,:), PrT(:,:), G2c(:), c4(:), mo1d(:)
    real(dp), allocatable :: u(:), w(:), e(:), z(:,:), MPV(:,:), yc4(:), zc(:)
    real(dp) :: DG, DC, kap, pJ, acc, wsum
    real(dp), external :: ddot
    integer(is) :: p, kk, j, bl, al, i, v, O1, O2

    O1 = nc + 1
    O2 = nc + 2
    allocate(d(nG), EO2V(nv), mO2V(nv), MCVt(nv,nc), PrT(nv,nc), G2c(nc), c4(nc), mo1d(nc))
    allocate(u(nc), w(nc), e(nc), z(nvirb,nc), MPV(nvirb,nv), yc4(nv), zc(nvirb))
    do p = 1, nG
       d(p) = Mx(p,p)
    enddo
    DG = product(d)
    DC = product(d(1:nc))
    kap = sGc*real((-1)**nc, dp)*isqrt2
    ! state-independent tables, all in the (virtual, core) layout of Yc
    do bl = 1, nv
       mO2V(bl) = Mx(O2, nc+2+bl)
       EO2V(bl) = 0.0_dp
       do kk = 1, nG
          EO2V(bl) = EO2V(bl) + Mx(O2,kk)*Mx(kk,nc+2+bl)/d(kk)
       enddo
    enddo
    do j = 1, nc
       do bl = 1, nv
          MCVt(bl,j) = Mx(j, nc+2+bl)
          PrT(bl,j) = 0.0_dp
          if (kt >= 1) PrT(bl,j) = Mx(j,nc+2+bl)/d(j)
          if (kt >= 2) then
             acc = 0.0_dp
             do kk = 1, nG
                acc = acc + Mx(j,kk)*Mx(kk,nc+2+bl)/d(kk)
             enddo
             PrT(bl,j) = PrT(bl,j) - (acc - Mx(j,nc+2+bl))/d(j)
          endif
       enddo
       G2c(j) = 0.0_dp
       do kk = 1, nc
          G2c(j) = G2c(j) + Mx(j,kk)*Mx(kk,O1)/d(kk)
       enddo
       c4(j) = 0.0_dp
       if (kt >= 1) c4(j) = -Mx(j,O1)/d(j)
       if (kt >= 2) c4(j) = c4(j) + G2c(j)/d(j)
       mo1d(j) = Mx(j,O1)/d(j)
    enddo
    do bl = 1, nv
       do al = 1, nvirb
          MPV(al,bl) = Mx(nc+al, nc+2+bl)
       enddo
    enddo

    g1 = 0.0_dp; g4 = 0.0_dp
    do v = 1, nstc
       ! u(j) = sum_b M(O2,b) Y(b,j), e(j) = sum_b E_O2V(b) Y(b,j), w(j) = sum_b M(j,b) Y(b,j)
       call dgemv('T', nv, nc, 1.0_dp, Yc(1,1,v), nv, mO2V, 1, 0.0_dp, u, 1)
       call dgemv('T', nv, nc, 1.0_dp, Yc(1,1,v), nv, EO2V, 1, 0.0_dp, e, 1)
       do j = 1, nc
          w(j) = ddot(nv, MCVt(1,j), 1, Yc(1,j,v), 1)
       enddo
       pJ = ddot(nv*nc, PrT, 1, Yc(1,1,v), 1)
       ! z(a,j) = sum_b M(a,b) Y(b,j); zc(a) = sum_j z(a,j) M(j,O1)/d(j); yc4(b) = sum_j Y(b,j) c4(j)
       call dgemm('N','N', nvirb, nc, nv, 1.0_dp, MPV, nvirb, Yc(1,1,v), nv, 0.0_dp, z, nvirb)
       call dgemv('N', nvirb, nc, 1.0_dp, z, nvirb, mo1d, 1, 0.0_dp, zc, 1)
       call dgemv('N', nv, nc, 1.0_dp, Yc(1,1,v), nv, c4, 1, 0.0_dp, yc4, 1)
       wsum = 0.0_dp
       do j = 1, nc
          wsum = wsum + w(j)/d(j)
       enddo
       ! g1
       g1(O2,v) = kap*real((-1)**nc, dp)*DG*pJ
       do i = 1, nG
          acc = 0.0_dp
          if (kt >= 2) then
             do j = 1, nc
                if (j == i) cycle
                acc = acc + (Mx(j,i)*u(j) - Mx(O2,i)*w(j))/d(j)
             enddo
             acc = real((-1)**i, dp)*kap*DG/d(i)*acc
          endif
          if (i <= nc) then
             if (kt >= 1) acc = acc + DG/d(i)*kap*real((-1)**(i-1), dp)*u(i)
             if (kt >= 2) acc = acc - DG/d(i)*kap*real((-1)**(i-1), dp)*(e(i) - Mx(O2,i)*w(i)/d(i))
          endif
          g1(i,v) = acc
       enddo
       ! g4
       g4(1,v) = kap*real((-1)**nc, dp)*DG*pJ
       do al = 2, nvirb
          acc = 0.0_dp
          if (kt >= 2) acc = Mx(nc+al,O1)*wsum - zc(al)
          if (al >= 3) acc = acc + Mx(nc+al, nc+al)*yc4(al-2)
          g4(al,v) = kap*real((-1)**nc, dp)*DC*acc
       enddo
    enddo
    deallocate(d, EO2V, mO2V, MCVt, PrT, G2c, c4, mo1d, u, w, e, z, MPV, yc4, zc)

  end subroutine cross_tlf

!######################################################################
! align_ket_ext: alignment of the ket MOs for extended objects: core and
! virtual classes by orthogonal Procrustes rotations as in align_ket,
! the SOMO pair only by sign flips diag(+-1) (an O1-O2 rotation or swap
! does not preserve the CV space built on G); offdiag = largest
! off-diagonal SOMO overlap (a measure of the rotation that was not
! aligned); sgn_o = det(U_O)
!######################################################################
  subroutine align_ket_ext(nmo, nc, nocca, nvirb, hb, pb, hk, pk, M, Uh, Up, sgn_o, offdiag, info)

    integer(is), intent(in)    :: nmo, nc, nocca, nvirb, hb(nocca), pb(nvirb), hk(nocca), pk(nvirb)
    real(dp), intent(inout)    :: M(nmo,nmo)
    real(dp), intent(out)      :: Uh(nocca,nocca), Up(nvirb,nvirb), sgn_o, offdiag
    integer(is), intent(out)   :: info

    real(dp), allocatable :: Uc(:,:), Uv(:,:)
    real(dp) :: s1, s2
    integer(is) :: nv

    nv = nvirb - 2
    allocate(Uc(nc,nc), Uv(nv,nv))
    info = 0
    call procrustes(nmo, nc, hb(1:nc), hk(1:nc), M, Uc, info)
    if (info == 0) call procrustes(nmo, nv, pb(3:nvirb), pk(3:nvirb), M, Uv, info)
    offdiag = max(abs(M(hb(nc+1), hk(nc+2))), abs(M(hb(nc+2), hk(nc+1))))
    s1 = sign(1.0_dp, M(hb(nc+1), hk(nc+1)))
    s2 = sign(1.0_dp, M(hb(nc+2), hk(nc+2)))
    M(:, hk(nc+1)) = s1*M(:, hk(nc+1))
    M(:, hk(nc+2)) = s2*M(:, hk(nc+2))
    Uh = 0.0_dp; Up = 0.0_dp
    Uh(1:nc,1:nc) = Uc
    Uh(nc+1,nc+1) = s1
    Uh(nc+2,nc+2) = s2
    Up(1,1) = s1
    Up(2,2) = s2
    Up(3:nvirb,3:nvirb) = Uv
    sgn_o = s1*s2
    deallocate(Uc, Uv)

  end subroutine align_ket_ext

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

  subroutine mrsf_state_overlap_c(nao_b, nao_k, nmo, occ_b, occ_k, mult, method, align, nb, nk, ncolb, ncolk, &
       Cb, Ck, Sao, xb, xk, S, ierr) bind(c, name='mrsf_state_overlap')
    integer(is), intent(in)  :: nao_b, nao_k, nmo, mult, method, align, nb, nk, ncolb, ncolk
    real(dp), intent(in)     :: occ_b(*), occ_k(*), Cb(*), Ck(*), Sao(*), xb(*), xk(*)
    real(dp), intent(out)    :: S(*)
    integer(is), intent(out) :: ierr
    call state_overlap(nao_b, nao_k, nmo, occ_b, occ_k, mult, method, align, nb, nk, ncolb, ncolk, Cb, Ck, Sao, &
         xb, xk, S, ierr)
  end subroutine mrsf_state_overlap_c

end module mrsf_overlap_interface
