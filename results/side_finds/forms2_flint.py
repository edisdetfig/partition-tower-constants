# Modular forms of level 2 on FLINT alone, for the towers of 3-colored partitions. 12 Nov 2026. The level-free tools
# (kernel_part, padic_kernel, restrict, ordinary_part, unit_root_matrix, newton_at_minus_one, …) are forms6_flint's.
#
# THE BASIS. For even k ≥ 4, M_k(Γ0(2)) has dimension ⌊k/4⌋ + 1 (= m + 1) and the basis H_j = G_k·u^j, j = 0..m, from
#     F  = η(2z)^16 / η(z)^8          weight 4, a simple zero at the cusp ∞, none elsewhere;
#     F₂ = 2·E₂(2z) − E₂(z)           weight 2, its one zero (of order 1/2) at the elliptic point;
#     G_k = F^(k/4) (k ≡ 0 mod 4) or F₂·F^((k−2)/4) (k ≡ 2 mod 4): all m zeros at ∞ (and 1/2 at the elliptic point);
#     u  = (η(z)/η(2z))^24 = 1/q·…    weight 0, a simple pole at ∞, a simple zero at the cusp 0.
# (Orders at the cusps by the standard eta-quotient formula; derived.) H_j = q^(m−j)·(1 + …), so a form is fixed by its
# coefficients at q⁰ … q^m (Sturm's bound for Γ0(2) is k/4) and those are its coordinates.
# HECKE OPERATORS modulo m: T_p f(n) = a(pn) + p^(k−1)·a(n/p) for odd p, U₂ f(n) = a(2n).
from flint import nmod_poly, fmpz_poly, fmpz_mat, fmpz
import forms6_flint as F6


def _series(k, L, mod):
    """φ-part: G_k/q^m and u·q as power series to length L, modulo mod (an int) or exactly (mod = None)."""
    def eta(d):                                                 # ∏(1 − q^{dn}) to length L
        c = [0] * L; c[0] = 1; j = 1
        while d * j * (3 * j - 1) // 2 < L:
            s = -1 if j % 2 else 1
            e = d * j * (3 * j - 1) // 2; c[e] += s
            e = d * j * (3 * j + 1) // 2
            if e < L: c[e] += s
            j += 1
        return c
    def sigma1(n):
        s, d = 0, 1
        while d * d <= n:
            if n % d == 0: s += d + (n // d if d * d != n else 0)
            d += 1
        return s
    if mod is None:
        P = fmpz_poly; mk = lambda c: fmpz_poly(c)
        mul = lambda a, b: a.mul_low(b, L)
        inv = lambda a: _inv_exact(a, L)
    else:
        mk = lambda c: F6.poly_mod([x % mod for x in c], mod)
        mul = lambda a, b: a.mul_low(b, L)
        inv = lambda a: a.inverse_series_trunc(L)
    E1, E2 = mk(eta(1)), mk(eta(2))
    Fq = mul(E2.pow_trunc(16, L), inv(E1.pow_trunc(8, L)))     # F / q
    uq = mul(E1.pow_trunc(24, L), inv(E2.pow_trunc(24, L)))    # u·q
    F2 = mk([1] + [24 * sigma1(n) - (48 * sigma1(n // 2) if n % 2 == 0 else 0) for n in range(1, L)])
    if k % 4 == 0: G = Fq.pow_trunc(k // 4, L)
    else: G = mul(F2, Fq.pow_trunc((k - 2) // 4, L))
    return G, uq, mul


def _inv_exact(a, L):
    """The inverse of an integer power series with constant term ±1, to length L."""
    c = [int(x) for x in a.coeffs()] + [0] * L; c0 = c[0]; assert c0 in (1, -1)
    b = [0] * L; b[0] = c0
    for n in range(1, L):
        s = sum(c[i] * b[n - i] for i in range(1, n + 1))
        b[n] = -s * c0
    return fmpz_poly(b)


def hecke_matrices2(k, mod, ops):
    """ops: ('T', p) with p odd, or ('U', 2). Returns {op: (m+1)×(m+1) matrix modulo mod} in coefficient coordinates."""
    m = k // 4; n = m + 1; P = max(p for kind, p in ops); L = P * m + 1
    G, uq, mul = _series(k, L, mod)
    psi = G; B = [0] * (n * n); C = {op: [0] * (n * n) for op in ops}
    pk = {p: pow(p, k - 1, mod) for kind, p in ops}
    for j in range(n):                                          # H_j = q^(m−j)·psi_j, psi_j = G·(u q)^j; column c = m − j,
        c = m - j; lc = psi.length()                            # so B is lower unitriangular
        get = lambda i: int(psi[i]) if 0 <= i < lc else 0
        for nn in range(c, n): B[nn * n + c] = get(nn - c) % mod
        for (kind, p) in ops:
            Cp = C[(kind, p)]
            for nn in range(n):
                v = get(p * nn - c)
                if kind == 'T' and nn % p == 0:
                    v2 = get(nn // p - c)
                    if v2: v += pk[p] * v2
                Cp[nn * n + c] = v % mod
        if j < n - 1: psi = mul(psi, uq)
    Bm = F6.mat_mod(n, n, B, mod)
    I = F6.identity(n, mod); Y = I - Bm; Binv = I + Y; covered = 2
    while covered < n:
        Y = Y * Y; Binv = Binv * (I + Y); covered *= 2
    return {op: F6.mat_mod(n, n, C[op], mod) * Binv for op in ops}


def class_on_lane2(k, ell, s2, N=4):
    """The class (newforms of weight k on Γ0(2) with U₂ = s2·2^((k−2)/2) exactly) on lane ℓ: the generalized eigenspace of
    U₂ for that value, lifted over the ℓ-adic integers (forms6_flint.kernel_part), then the forms whose U₂ is EXACTLY the
    value (an ℓ-adic kernel; a form of level 1 can agree with it on lane ℓ only). Returns (X, R, d): a basis on lane ℓ with
    the identity at the rows R, and the dimension; X is None when the class is empty."""
    from flint import nmod_mat
    mm = ell ** N; n = k // 4 + 1
    H = hecke_matrices2(k, mm, [('U', 2)]); c2 = s2 * 2 ** ((k - 2) // 2) % mm
    Y, R, r = F6.kernel_part(H[('U', 2)] - F6.identity(n, mm) * c2, ell, N)
    if r == 0: return None, None, 0
    U2W = F6.restrict(H[('U', 2)], Y, R, r)
    ker, Np = F6.padic_kernel([[(int(U2W[i, j]) - (c2 if i == j else 0)) % mm for j in range(r)] for i in range(r)], r, r, ell, N)
    d = len(ker)
    if d == 0: return None, None, 0
    Yl = nmod_mat(n, r, [int(v) % ell for v in Y.entries()], ell)
    Kl = nmod_mat(r, d, [ker[c][i] % ell for i in range(r) for c in range(d)], ell)
    X, Rx = F6.normalized_basis(Yl * Kl, d, ell)
    return X, Rx, d


def class_basis_exact2(k, s2):
    """The newforms of M_k(Γ0(2)) with U₂ = s2·2^((k−2)/2) EXACTLY (no form of level 1 and no Eisenstein series has this
    value, by Deligne's bound): an integer basis (fmpz_mat, (m+1)×d) in coefficient coordinates."""
    m = k // 4; n = m + 1; L = 2 * m + 1
    G, uq, mul = _series(k, L, None)
    psi = G; B = [[0] * n for _ in range(n)]; C2 = [[0] * n for _ in range(n)]
    for j in range(n):
        c = m - j; lc = psi.length(); get = lambda i: int(psi[i]) if 0 <= i < lc else 0
        for nn in range(c, n): B[nn][c] = get(nn - c)
        for nn in range(n): C2[nn][c] = get(2 * nn - c)
        if j < n - 1: psi = mul(psi, uq)
    Bm = fmpz_mat(B); I = fmpz_mat([[1 if a == b else 0 for b in range(n)] for a in range(n)])
    Y = I - Bm; Binv = I + Y; covered = 2
    while covered < n:
        Y = Y * Y; Binv = Binv * (I + Y); covered *= 2
    U2 = fmpz_mat(C2) * Binv
    A = U2 - I * (s2 * 2 ** ((k - 2) // 2))
    K, nul = A.nullspace()
    return fmpz_mat([[int(K[a, b]) for b in range(nul)] for a in range(n)]), nul
