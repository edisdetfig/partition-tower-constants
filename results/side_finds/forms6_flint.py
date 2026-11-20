# Modular forms of level 6 on FLINT alone (no PARI), for the partition tower's families. 12 Nov 2026.
#
# THE BASIS. For even k ≥ 4, M_k(Γ0(6)) has dimension k + 1 and the basis H_m = G_k·t^(k−m), m = 0..k, built from two
# eta-quotients (orders at the four cusps checked by the standard formula):
#     G_k = (η(z)·η(6z)⁶ / (η(2z)²·η(3z)³))^k = q^k·φ   weight k, all k zeros at the cusp ∞, none elsewhere;
#     t   = η(z)⁵·η(3z) / (η(2z)·η(6z)⁵)      = u/q   weight 0, a simple pole at ∞, a simple zero at the cusp 0.
# H_m = q^m·φ·u^(k−m) starts with q^m, so a form is fixed by its first k + 1 coefficients (Sturm's bound for Γ0(6) is k):
# the coordinates of a form are its coefficients at q⁰ … q^k.
# THE HECKE OPERATORS, on q-expansions modulo m = ℓ^N: T_p f(n) = a(pn) + p^(k−1)·a(n/p) for p ∤ 6, U_p f(n) = a(pn) for
# p | 6, as (k+1)×(k+1) matrices in those coordinates.
# THE TOWER'S PIECE. A Hecke-stable piece is found over the field with ℓ elements (FLINT's exact linear algebra mod a
# prime) as the kernel of g(A) for a combination A of Hecke operators and the polynomial g of the piece's numbers, then
# lifted to modulo ℓ^N by Newton's method (Hensel's lemma for an invariant subspace, with the Jacobian inverted mod ℓ).
import time
from flint import nmod_poly, nmod_mat, fmpz_mod_poly_ctx, fmpz_mod_ctx, fmpz_mod_mat

NMOD_MAX = 2 ** 62                       # below this the fast word-size types; above it FLINT's multi-precision ones
_ctx = {}


def poly_mod(c, m):
    if m < NMOD_MAX: return nmod_poly(c, m)
    if ('p', m) not in _ctx: _ctx['p', m] = fmpz_mod_poly_ctx(m)
    return _ctx['p', m](c)


def mat_mod(r, c, entries, m):
    if m < NMOD_MAX: return nmod_mat(r, c, entries, m)
    if ('m', m) not in _ctx: _ctx['m', m] = fmpz_mod_ctx(m)
    return fmpz_mod_mat(r, c, entries, _ctx['m', m])


def mod_of(A):
    return int(A.modulus())


def eta_series(d, L, m):
    """∏_{n≥1} (1 − q^{dn}) modulo m, to length L (Euler's pentagonal series)."""
    c = [0] * L; c[0] = 1; k = 1
    while d * k * (3 * k - 1) // 2 < L:
        s = m - 1 if k % 2 else 1
        e = d * k * (3 * k - 1) // 2; c[e] = (c[e] + s) % m
        e = d * k * (3 * k + 1) // 2
        if e < L: c[e] = (c[e] + s) % m
        k += 1
    return poly_mod(c, m)


def identity(n, m):
    return mat_mod(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)], m)


def hecke_matrices(k, m, ops):
    """ops: list of ('T', p) with p ∤ 6 or ('U', p) with p | 6. Returns {op: (k+1)×(k+1) nmod_mat mod m}."""
    n = k + 1; L = max(p for kind, p in ops) * k + 1
    E1, E2, E3, E6 = (eta_series(d, L, m) for d in (1, 2, 3, 6))
    mul = lambda a, b: a.mul_low(b, L)
    u = mul(mul(E1.pow_trunc(5, L), E3), mul(E2, E6.pow_trunc(5, L)).inverse_series_trunc(L))
    g1 = mul(mul(E1, E6.pow_trunc(6, L)), mul(E2.pow_trunc(2, L), E3.pow_trunc(3, L)).inverse_series_trunc(L))
    psi = g1.pow_trunc(k, L)                                  # φ; then φ·u^i, and H_m = q^m·φ·u^(k−m)
    B = [0] * (n * n); C = {op: [0] * (n * n) for op in ops}
    pk = {p: pow(p, k - 1, m) for kind, p in ops}
    for i in range(n):
        mm = k - i; lc = psi.length()
        get = lambda j: int(psi[j]) if 0 <= j < lc else 0      # read only the coefficients needed
        for nn in range(mm, n): B[nn * n + mm] = get(nn - mm)
        for (kind, p) in ops:
            Cp = C[(kind, p)]
            for nn in range(n):
                v = get(p * nn - mm)
                if kind == 'T' and nn % p == 0:
                    v2 = get(nn // p - mm)
                    if v2: v = (v + pk[p] * v2) % m
                Cp[nn * n + mm] = v
        if i < n - 1: psi = mul(psi, u)
    Bm = mat_mod(n, n, B, m)
    # B is lower unitriangular: B⁻¹ = (I + Y)(I + Y²)(I + Y⁴)… with Y = I − B nilpotent
    I = identity(n, m); Y = I - Bm; Binv = I + Y; covered = 2
    while covered < n:
        Y = Y * Y; Binv = Binv * (I + Y); covered *= 2
    return {op: mat_mod(n, n, C[op], m) * Binv for op in ops}


def class_basis_exact(k, s2, s3):
    """The newforms of M_k(Γ0(6)) with U₂ = s2·2^((k−2)/2) and U₃ = s3·3^((k−2)/2) EXACTLY (over the integers; no form of
    lower level or Eisenstein series has these eigenvalues, by Deligne's bound): an integer basis (fmpz_mat, (k+1)×d) in
    coefficient coordinates. Exact series to 3k + 1 terms; 1/∏(1 − q^{dn}) is Σ p(j)·q^{dj}."""
    from flint import fmpz_poly, fmpz_mat, fmpz
    n = k + 1; Lx = 3 * k + 1
    def eta(d):
        c = [0] * Lx; c[0] = 1; j = 1
        while d * j * (3 * j - 1) // 2 < Lx:
            s = -1 if j % 2 else 1
            e = d * j * (3 * j - 1) // 2; c[e] += s
            e = d * j * (3 * j + 1) // 2
            if e < Lx: c[e] += s
            j += 1
        return fmpz_poly(c)
    def peta(d):                                              # 1/∏(1 − q^{dn})
        c = [0] * Lx
        for j in range(0, (Lx - 1) // d + 1): c[d * j] = int(fmpz(j).partitions_p())
        return fmpz_poly(c)
    mul = lambda a, b: a.mul_low(b, Lx)
    E1, E3, E6 = eta(1), eta(3), eta(6); P2, P3, P6 = peta(2), peta(3), peta(6)
    u = mul(mul(mul(E1.pow_trunc(5, Lx), E3), P2), P6.pow_trunc(5, Lx))
    g1 = mul(mul(mul(E1, E6.pow_trunc(6, Lx)), P2.pow_trunc(2, Lx)), P3.pow_trunc(3, Lx))
    psi = g1.pow_trunc(k, Lx)
    B = [[0] * n for _ in range(n)]; C2 = [[0] * n for _ in range(n)]; C3 = [[0] * n for _ in range(n)]
    for i in range(n):
        mm = k - i; lc = psi.length(); get = lambda j: int(psi[j]) if 0 <= j < lc else 0
        for nn in range(mm, n): B[nn][mm] = get(nn - mm)
        for nn in range(n): C2[nn][mm] = get(2 * nn - mm); C3[nn][mm] = get(3 * nn - mm)
        if i < n - 1: psi = mul(psi, u)
    Bm = fmpz_mat(B); I = fmpz_mat([[1 if a == b else 0 for b in range(n)] for a in range(n)])
    Y = I - Bm; Binv = I + Y; covered = 2
    while covered < n:
        Y = Y * Y; Binv = Binv * (I + Y); covered *= 2
    U2 = fmpz_mat(C2) * Binv; U3 = fmpz_mat(C3) * Binv
    c2 = s2 * 2 ** ((k - 2) // 2); c3 = s3 * 3 ** ((k - 2) // 2)
    A2, A3 = U2 - I * c2, U3 - I * c3
    S = fmpz_mat([[int(A2[a, b]) for b in range(n)] for a in range(n)] + [[int(A3[a, b]) for b in range(n)] for a in range(n)])
    K, nul = S.nullspace()
    return fmpz_mat([[int(K[a, b]) for b in range(nul)] for a in range(n)]), nul


def reduce_mat(A, ell):
    return nmod_mat(A.nrows(), A.ncols(), [int(x) % ell for x in A.entries()], ell)


def eval_poly(g, A):
    """g(A) for a polynomial g (nmod_poly) and a square nmod_mat A of the same modulus."""
    co = [int(x) for x in g.coeffs()]; n = A.nrows(); I = identity(n, mod_of(A)); R = I * co[-1]
    for c in reversed(co[:-1]): R = R * A + I * c
    return R


def pivot_rows(X, ell):
    """Rows of the n×d matrix X (nmod_mat over F_ℓ, rank d) whose d×d minor is invertible."""
    Rr, rank = X.transpose().rref(); d = X.ncols(); piv = []
    for i in range(rank):
        for j in range(X.nrows()):
            if int(Rr[i, j]) != 0: piv.append(j); break
    return piv[:d]


def restrict(T, X, R, d):
    """The d×d matrix of T on span(X), when X (n×d) has the identity at the rows R and T preserves span(X)."""
    TX = T * X
    return mat_mod(d, d, [int(TX[R[a], b]) for a in range(d) for b in range(d)], mod_of(T))


def normalized_basis(K, cols, ell):
    """A basis (n×cols) of the column span of K's first `cols` columns, with the identity at its pivot rows."""
    n = K.nrows()
    X = nmod_mat(n, cols, [int(K[i, j]) for i in range(n) for j in range(cols)], ell)
    R = pivot_rows(X, ell)
    XR = nmod_mat(cols, cols, [int(X[R[a], b]) for a in range(cols) for b in range(cols)], ell).inv()
    return X * XR, R


def generalized_kernel(B):
    """Kernel of B^e for e large enough (B over a prime field): returns (basis matrix, dimension)."""
    P = B; K, nul = P.nullspace()
    while True:
        P2 = P * P; K2, nul2 = P2.nullspace()
        if nul2 == nul: return K, nul
        P, K, nul = P2, K2, nul2


def lift_invariant(A, gbar, ell, N):
    """A: (n×n) nmod_mat mod ℓ^N; gbar: polynomial over F_ℓ of degree d. Find X (n×d) with the identity at rows R and
    M (d×d) such that A·X ≡ X·M (mod ℓ^N), lifting the kernel of gbar(A mod ℓ). Returns (X, M, R) or (None, info)."""
    n = A.nrows(); m = ell ** N; Ab = reduce_mat(A, ell)
    G = eval_poly(gbar, Ab); K, d = generalized_kernel(G)  # the whole generalized eigenspace (a level-raising congruence
    if d == 0: return None, ('empty',)                     # can add old forms with the same numbers on lane ℓ)
    Xb, R = normalized_basis(K, d, ell)
    Mb = restrict(Ab, Xb, R, d)
    Rset = set(R); free = [i for i in range(n) if i not in Rset]; pos = {i: t for t, i in enumerate(free)}
    nd = n * d; J = [0] * (nd * nd); off = (n - d) * d
    Al = [int(x) for x in Ab.entries()]; Ml = [int(x) for x in Mb.entries()]; Xl = [int(x) for x in Xb.entries()]
    for i in free:                                           # unknown ΔX[i][c]
        for c in range(d):
            col = pos[i] * d + c
            for i2 in range(n):                              # (A ΔX)[i2][c] += A[i2][i]·ΔX[i][c]
                v = Al[i2 * n + i]
                if v: J[(i2 * d + c) * nd + col] = (J[(i2 * d + c) * nd + col] + v) % ell
            for c2 in range(d):                              # −(ΔX M)[i][c2] −= ΔX[i][c]·M[c][c2]
                v = Ml[c * d + c2]
                if v: J[(i * d + c2) * nd + col] = (J[(i * d + c2) * nd + col] - v) % ell
    for a in range(d):                                       # unknown ΔM[a][b]: −(X ΔM)[i][b] −= X[i][a]·ΔM[a][b]
        for b in range(d):
            col = off + a * d + b
            for i in range(n):
                v = Xl[i * d + a]
                if v: J[(i * d + b) * nd + col] = (J[(i * d + b) * nd + col] - v) % ell
    Jinv = nmod_mat(nd, nd, J, ell).inv()
    X = [Xl[i * d + c] for i in range(n) for c in range(d)]; M = list(Ml)
    for j in range(1, N):
        Xm = mat_mod(n, d, X, m); Mm = mat_mod(d, d, M, m)
        res = [int(v) for v in (A * Xm - Xm * Mm).entries()]
        pj = ell ** j
        if any(v % pj for v in res): return None, ('not divisible at step', j)
        b = nmod_mat(nd, 1, [(v // pj) % ell for v in res], ell)
        delta = [(-int(v)) % ell for v in (Jinv * b).entries()]
        for i in free:
            for c in range(d): X[i * d + c] = (X[i * d + c] + pj * delta[pos[i] * d + c]) % m
        for a in range(d):
            for c in range(d): M[a * d + c] = (M[a * d + c] + pj * delta[off + a * d + c]) % m
    Xm = mat_mod(n, d, X, m); Mm = mat_mod(d, d, M, m)
    if any(int(v) for v in (A * Xm - Xm * Mm).entries()): return None, ('residual not 0 at the end',)
    return (Xm, Mm, R), None


def padic_kernel(B, rows, cols, ell, N):
    """Kernel of a small integer matrix B (list of rows, entries mod ℓ^N) over the ℓ-adic integers: pivoting on the entry
    of least ℓ-adic valuation, column operations tracked in V. Returns (list of kernel columns, each a list of ints mod ℓ^N,
    with B·x ≡ 0 mod ℓ^N and x not ≡ 0 mod ℓ), and the precision N − (largest pivot valuation) to which they are the
    true kernel."""
    m = ell ** N; A = [[x % m for x in row] for row in B]; V = [[1 if i == j else 0 for j in range(cols)] for i in range(cols)]
    def val(x):
        if x % m == 0: return N
        v = 0
        while x % ell == 0: x //= ell; v += 1
        return v
    pivv = []; s = 0
    while s < min(rows, cols):
        best = None
        for i in range(s, rows):
            for j in range(s, cols):
                v = val(A[i][j])
                if best is None or v < best[0]: best = (v, i, j)
        if best[0] >= N: break
        v, i, j = best
        A[s], A[i] = A[i], A[s]
        for row in A: row[s], row[j] = row[j], row[s]
        for row in V: row[s], row[j] = row[j], row[s]
        p = A[s][s]; uinv = pow(p // ell ** v, -1, m)
        for j2 in range(cols):
            if j2 != s and A[s][j2]:
                f = (A[s][j2] // ell ** v) * uinv % m
                for row in A: row[j2] = (row[j2] - f * row[s]) % m
                for row in V: row[j2] = (row[j2] - f * row[s]) % m
        for i2 in range(rows):
            if i2 != s and A[i2][s]:
                f = (A[i2][s] // ell ** v) * uinv % m
                A[i2] = [(a - f * b) % m for a, b in zip(A[i2], A[s])]
        pivv.append(v); s += 1
    return [[V[i][c] for i in range(cols)] for c in range(s, cols)], N - max(pivv, default=0)


def mat_poly(co, T):
    """Σ co[i]·T^i (co lowest first, ints) by Paterson–Stockmeyer: about 2√deg matrix products."""
    n = T.nrows(); md = mod_of(T); I = identity(n, md); deg = len(co) - 1
    if deg <= 0: return I * (co[0] % md if co else 0)
    s = max(1, int(deg ** 0.5)); P = [I, T]
    for _ in range(2, s + 1): P.append(P[-1] * T)
    Ts = P[s]; blocks = [co[i:i + s] for i in range(0, deg + 1, s)]
    R = None
    for blk in reversed(blocks):
        B = I * 0
        for j, c in enumerate(blk):
            if c % md: B = B + P[j] * (c % md)
        R = B if R is None else R * Ts + B
    return R


def inverse_mod(Mt, ell, N):
    """The inverse of a square nmod_mat modulo ℓ^N (its reduction mod ℓ invertible), by Newton–Schulz from the mod-ℓ inverse."""
    n = Mt.nrows(); m = ell ** N
    Z = mat_mod(n, n, [int(x) for x in reduce_mat(Mt, ell).inv().entries()], m); I2 = identity(n, m) * 2
    prec = 1
    while prec < N:
        Z = Z * (I2 - Mt * Z); prec *= 2
    return Z


def ordinary_part(Tl, ell, N):
    """The ordinary part of the Hecke operator T_ℓ (its unit eigenvalues): the idempotent e (a polynomial in T_ℓ, lifted
    from mod ℓ by Newton's e ← 3e² − 2e³), and a basis Y (n×r) of its image with the identity at the rows R."""
    n = Tl.nrows(); m = ell ** N; Tb = reduce_mat(Tl, ell); chi = Tb.charpoly()
    co = [int(x) for x in chi.coeffs()]; a = next(i for i, c in enumerate(co) if c)          # chi = x^a·h
    h = nmod_poly(co[a:], ell); xa = nmod_poly([0] * a + [1], ell)
    gg, s, t = xa.xgcd(h)                                    # s·x^a + t·h = gg (a unit)
    b = s * pow(int(gg.coeffs()[0]), -1, ell)                 # b·x^a ≡ 1 (mod h)
    Ta = identity(n, m); base = Tl; e_ = a
    while e_:
        if e_ & 1: Ta = Ta * base
        base = base * base; e_ >>= 1
    e = mat_poly([int(x) for x in b.coeffs()], Tl) * Ta
    while True:
        e2 = e * e; en = e2 * 3 - e2 * e * 2
        if en == e: break
        e = en
    r = n - a
    if r == 0: return None, None, 0
    eb = reduce_mat(e, ell); Rr, rank = eb.rref()
    cols = []
    for i in range(rank):
        for j in range(n):
            if int(Rr[i, j]) != 0: cols.append(j); break
    Y = mat_mod(n, r, [int(e[i, j]) for i in range(n) for j in cols[:r]], m)
    R = pivot_rows(reduce_mat(Y, ell), ell)
    YR = mat_mod(r, r, [int(Y[R[a_], b_]) for a_ in range(r) for b_ in range(r)], m)
    return Y * inverse_mod(YR, ell, N), R, r


def image_basis(P, r, ell, N):
    """A basis (n×r) of the image of the projector P (rank r) modulo ℓ^N, with the identity at the rows R."""
    n = P.nrows(); m = ell ** N; Rr, rank = reduce_mat(P, ell).rref(); cols = []
    for i in range(rank):
        for j in range(n):
            if int(Rr[i, j]) != 0: cols.append(j); break
    Y = mat_mod(n, r, [int(P[i, j]) for i in range(n) for j in cols[:r]], m)
    R = pivot_rows(reduce_mat(Y, ell), ell)
    YR = mat_mod(r, r, [int(Y[R[a], b]) for a in range(r) for b in range(r)], m)
    return Y * inverse_mod(YR, ell, N), R


def kernel_part(B, ell, N):
    """The generalized kernel of B over the ℓ-adic integers, modulo ℓ^N (the part where B is topologically nilpotent;
    Fitting's splitting): I − e, with e = b(B)·B^a the idempotent onto the part where B is invertible, lifted by
    e ← 3e² − 2e³ exactly as in ordinary_part. Returns (Y n×a with the identity at the rows R, R, a)."""
    n = B.nrows(); m = ell ** N; co = [int(x) for x in reduce_mat(B, ell).charpoly().coeffs()]
    a = next(i for i, c in enumerate(co) if c); I = identity(n, m)
    if a == 0: return None, None, 0
    if a == n: return I, list(range(n)), n
    h = nmod_poly(co[a:], ell); xa = nmod_poly([0] * a + [1], ell)
    gg, s, t = xa.xgcd(h); b = s * pow(int(gg.coeffs()[0]), -1, ell)
    Ba = I; base = B; e_ = a
    while e_:
        if e_ & 1: Ba = Ba * base
        base = base * base; e_ >>= 1
    e = mat_poly([int(x) for x in b.coeffs()], B) * Ba
    while True:
        e2 = e * e; en = e2 * 3 - e2 * e * 2
        if en == e: break
        e = en
    Y, R = image_basis(I - e, a, ell, N)
    return Y, R, a


def class_padic(H, k, ell, N, s2, s3):
    """The class at weight k over the ℓ-adic integers: the forms with U₂ = s2·2^((k−2)/2) and U₃ = s3·3^((k−2)/2) EXACTLY
    (new at 2 and 3; by Deligne's bound no form of lower level and no Eisenstein series has these values), separated from
    the forms that share them only on lane ℓ (a level-raising congruence). H: Hecke matrices modulo ℓ^N with ('U', 2) and
    ('U', 3). Returns (W n×r: the joint generalized eigenspace, its rows RW, r; Wn n×d: the class, its rows Rn, d; N′: the
    class's precision), or None when the joint space is empty."""
    n = k + 1; m = ell ** N
    c2 = s2 * 2 ** ((k - 2) // 2) % m; c3 = s3 * 3 ** ((k - 2) // 2) % m
    Y1, R1, r1 = kernel_part(H[('U', 2)] - identity(n, m) * c2, ell, N)
    if r1 == 0: return None
    Z, Rz, r = kernel_part(restrict(H[('U', 3)], Y1, R1, r1) - identity(r1, m) * c3, ell, N)
    if r == 0: return None
    W = Y1 * Z; RW = [R1[j] for j in Rz]
    U2W = restrict(H[('U', 2)], W, RW, r); U3W = restrict(H[('U', 3)], W, RW, r)
    Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)] + \
         [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(r)] for i in range(r)]
    ker, Np = padic_kernel(Bm, 2 * r, r, ell, N)
    d = len(ker)
    if d == 0 or Np < 1: return W, RW, r, None, None, d, Np
    mp = ell ** Np
    Wn = mat_mod(n, r, [int(v) % mp for v in W.entries()], mp) * mat_mod(r, d, [ker[c][i] % mp for i in range(r) for c in range(d)], mp)
    Rn = pivot_rows(reduce_mat(Wn, ell), ell)
    Wn = Wn * inverse_mod(mat_mod(d, d, [int(Wn[Rn[a], b]) for a in range(d) for b in range(d)], mp), ell, Np)
    return W, RW, r, Wn, Rn, d, Np


def berkowitz(A, m):
    """The characteristic polynomial of the square matrix A (list of rows of ints) modulo any m, without division
    (Berkowitz). Returns its coefficients lowest first, ending in 1."""
    n = len(A); p = [1]                                       # highest degree first while building
    for k in range(n):
        r = [A[k][i] % m for i in range(k)]; v = [A[i][k] % m for i in range(k)]; vals = []
        for _ in range(k):
            vals.append(sum(a * b for a, b in zip(r, v)) % m)
            v = [sum(A[i][l] * v[l] for l in range(k)) % m for i in range(k)]
        T = [1, (-A[k][k]) % m] + [(-x) % m for x in vals]
        p = [sum(T[i - j] * p[j] for j in range(len(p)) if 0 <= i - j < len(T)) % m for i in range(k + 2)]
    return p[::-1]


def slope_root(chi, ell, s, N, pick=None):
    """chi: a monic polynomial (lowest first) known modulo ℓ^N. Its Newton polygon, and the root of valuation exactly s when
    the polygon has one segment of slope s of length 1: Hensel's lemma on h(y) = chi(ℓ^s·y)/ℓ^c. Returns (root mod ℓ^prec,
    prec, the polygon's slopes with their lengths) or (None, 0, slopes). With several roots of slope s, pick (a residue)
    chooses the one ≡ pick·… on lane ℓ (y = root/ℓ^s ≡ pick), when it is a simple root of the segment's polynomial."""
    d = len(chi) - 1; m = ell ** N
    def val(x):
        x %= m
        if x == 0: return N
        v = 0
        while x % ell == 0: x //= ell; v += 1
        return v
    vs = [val(c) for c in chi]
    hull = [0]                                                 # lower convex hull of the points (j, v_j)
    for j in range(1, d + 1):
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            if (vs[b] - vs[a]) * (j - a) >= (vs[j] - vs[a]) * (b - a): hull.pop()
            else: break
        hull.append(j)
    from fractions import Fraction
    slopes = [(Fraction(vs[a] - vs[b], b - a), b - a) for a, b in zip(hull, hull[1:])]
    seg = [(a, b) for a, b in zip(hull, hull[1:]) if Fraction(vs[a] - vs[b], b - a) == s]
    if len(seg) != 1 or (seg[0][1] - seg[0][0] != 1 and pick is None): return None, 0, slopes
    j1, j2 = seg[0]; c = vs[j1] + s * j1; prec = N - c
    if prec <= 0: return None, 0, slopes
    mp = ell ** prec
    h = []
    for j in range(d + 1):
        e = s * j - c
        h.append((chi[j] * ell ** e) % mp if e >= 0 else (chi[j] % m) // ell ** (-e) % mp)
    if j2 - j1 == 1: y = (-(h[j1] % ell) * pow(h[j1 + 1] % ell, -1, ell)) % ell
    else:                                                      # several roots of slope s: the one ≡ pick on lane ℓ
        res = nmod_poly([h[j] % ell for j in range(j1, j2 + 1)], ell)
        rts = {int(rt): e for rt, e in res.roots()}
        if pick % ell not in rts or rts[pick % ell] != 1: return None, 0, slopes
        y = pick % ell
    for _ in range(prec + 2):
        hy = sum(h[j] * pow(y, j, mp) for j in range(d + 1)) % mp
        dh = sum(j * h[j] * pow(y, j - 1, mp) for j in range(1, d + 1)) % mp
        y = (y - hy * pow(dh, -1, mp)) % mp
    return ell ** s * y % (ell ** (prec + s)), prec + s, slopes


def slope_factor(chi, ell, s, N):
    """chi monic (lowest first) known modulo ℓ^N, whose LAST Newton-polygon segment has slope s (all roots of valuation ≥ s,
    the smallest exactly s). Returns (G, prec): G monic (lowest first) whose roots are root/ℓ^s for the roots of valuation
    exactly s, modulo ℓ^prec (Hensel's lemma on h(y) = chi(ℓ^s·y)/ℓ^c ≡ y^(d−e)·g(y) on lane ℓ), or (None, 0)."""
    d = len(chi) - 1; m = ell ** N
    def val(x):
        x %= m
        if x == 0: return N
        v = 0
        while x % ell == 0: x //= ell; v += 1
        return v
    vs = [val(c) for c in chi]
    if not all(vs[d - j] >= s * j for j in range(1, d + 1)): return None, 0
    on = [e for e in range(1, d + 1) if vs[d - e] == s * e]          # the points on the line of slope s through (d, 0)
    if not on: return None, 0
    e = max(on)
    j1 = d - e; c = vs[j1] + s * j1; prec = N - c
    if prec <= 0: return None, 0
    mp = ell ** prec; h = []
    for j in range(d + 1):
        ee = s * j - c
        h.append((chi[j] * ell ** ee) % mp if ee >= 0 else (chi[j] % m) // ell ** (-ee) % mp)
    g = [h[j] % ell for j in range(j1, d + 1)]
    if g[0] == 0: return None, 0
    cof = [0] * j1 + [1]
    G, Hc = hensel_factor(h, g, cof, ell, prec)
    return [x % mp for x in G], prec


def charpoly_small(M, d, m):
    """Characteristic polynomial (lowest coefficient first, ending in 1) of a small d×d matrix mod m (Faddeev–LeVerrier;
    needs 1..d invertible mod m)."""
    A = [[int(M[i, j]) for j in range(d)] for i in range(d)]
    c = [0] * (d + 1); c[d] = 1; Mk = [[0] * d for _ in range(d)]
    for t in range(1, d + 1):
        Mk = [[(sum(A[i][l] * Mk[l][j] for l in range(d)) + (c[d - t + 1] if i == j else 0)) % m for j in range(d)] for i in range(d)]
        AM = [[sum(A[i][l] * Mk[l][j] for l in range(d)) % m for j in range(d)] for i in range(d)]
        c[d - t] = (-pow(t, -1, m) * sum(AM[i][i] for i in range(d))) % m
    return c


def unit_root_matrix(Ml, d, ell, N, k):
    """U with U² − Ml·U + ℓ^(k−1) = 0 and U ≡ Ml (mod ℓ): the U_ℓ-action on the ordinary ℓ-stabilizations."""
    m = ell ** N; Mm = mat_mod(d, d, [int(x) for x in Ml.entries()], m)
    if k - 1 >= N: return Mm
    U = Mm; pk = pow(ell, k - 1, m)
    for _ in range(N + 2):
        c = charpoly_small(U, d, m)                          # U⁻¹ = −(U^(d−1) + c_{d−1}U^(d−2) + … + c₁)/c₀
        I = identity(d, m); P = I * 0; Up = I
        acc = I * c[1]
        for t in range(2, d + 1):
            Up = Up * U; acc = acc + Up * c[t]
        Uinv = acc * ((-pow(c[0], -1, m)) % m)
        U = Mm - Uinv * pk
    return U


def newton_at_minus_one(vals, ell, prec):
    """Newton's series through the points x = 0..X (values mod ℓ^prec) evaluated at x = −1, and the differences'
    ℓ-adic valuations."""
    M = ell ** prec; X = len(vals) - 1
    diffs = []
    for nn in range(X + 1):
        s = 0; b = 1
        for i in range(nn + 1):                              # Δⁿ = Σ (−1)^(n−i) C(n, i) v_i
            s += (-1) ** (nn - i) * b * vals[i]; b = b * (nn - i) // (i + 1)
        diffs.append(s % M)
    val = []
    for dd in diffs:
        v = 0
        if dd == 0: val.append(prec); continue
        while dd % ell == 0: dd //= ell; v += 1
        val.append(v)
    return sum((-1) ** nn * diffs[nn] for nn in range(X + 1)) % M, val, diffs


def hensel_factor(f, g, h, ell, N):
    """f, g, h: coefficient lists (lowest first) with f ≡ g·h (mod ℓ), g monic, gcd(g, h) = 1 mod ℓ. Returns G, H with
    f ≡ G·H (mod ℓ^N), G monic of the same degree as g."""
    m = ell ** N
    P = lambda c: nmod_poly([x % ell for x in c], ell)
    one, s, t = P(g).xgcd(P(h))                              # s·g + t·h = gcd = 1 (made monic)
    inv1 = pow(int(one.coeffs()[0]), -1, ell); s, t = s * inv1, t * inv1
    G, H = [x % m for x in g], [x % m for x in h]
    for j in range(1, N):
        pj = ell ** j
        GH = [0] * (len(G) + len(H) - 1)
        for a, x in enumerate(G):
            for b, y in enumerate(H): GH[a + b] += x * y
        diff = [((f[i] if i < len(f) else 0) - (GH[i] if i < len(GH) else 0)) % m for i in range(max(len(f), len(GH)))]
        if any(x % pj for x in diff): raise ValueError('not a factorization mod ℓ^%d' % j)
        e = P([x // pj for x in diff])
        te = t * e; qg, dg = divmod(te, P(g))
        dh = s * e + qg * P(h)
        dgl = [int(x) for x in dg.coeffs()]; dhl = [int(x) for x in dh.coeffs()]
        G = [(G[i] + pj * (dgl[i] if i < len(dgl) else 0)) % m for i in range(len(G))]
        H = [(H[i] + pj * (dhl[i] if i < len(dhl) else 0)) % m for i in range(max(len(H), len(dhl)))]
        while len(H) > len(h) and H[-1] == 0: H.pop()
    return G, H


def padic_sqrt(a, ell, N):
    """A square root of the unit a (a square mod ℓ) modulo ℓ^N."""
    m = ell ** N; r = int(nmod_poly([(-a) % ell, 0, 1], ell).roots()[0][0])
    for _ in range(N + 1): r = (r + a * pow(r, -1, m)) * pow(2, -1, m) % m
    return r
