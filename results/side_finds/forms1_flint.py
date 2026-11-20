# Modular forms of level 1 on FLINT alone, for towers of weakly holomorphic forms of level 1 (1/Δ = 24-colored partitions,
# j, 1/E₄, …). 12 Nov 2026. The level-free tools (ordinary_part, unit_root_matrix, charpoly_small, …) are forms6_flint's.
#
# THE BASIS. M_k(SL₂(ℤ)) for even k ≥ 0 has dimension d = ⌊k/12⌋ (+1 unless k ≡ 2 mod 12). With w = k − 12(d − 1) one of
# 0, 4, 6, 8, 10, 14 and E_w the product of E₄, E₆ of that weight (E₀ = 1, E₈ = E₄², E₁₀ = E₄E₆, E₁₄ = E₄²E₆):
#     G_k = Δ^(d−1)·E_w = q^(d−1)·(1 + …),   j = E₄³/Δ = 1/q·(1 + …),   H_i = G_k·j^i = q^(d−1−i)·(1 + …), i = 0..d−1.
# A form is fixed by its coefficients at q⁰ … q^(d−1) (Sturm's bound for level 1 is k/12); those are its coordinates.
# HECKE OPERATORS modulo m: T_p f(n) = a(pn) + p^(k−1)·a(n/p).
import forms6_flint as F6


def dim_mk(k):
    if k < 0 or k % 2 or k == 2: return 0
    return k // 12 + (0 if k % 12 == 2 else 1)


def _sigma_list(L, e):
    s = [0] * L
    for d in range(1, L):
        for n in range(d, L, d): s[n] += d ** e
    return s


def hecke_matrices1(k, mod, primes):
    """{p: d×d matrix of T_p on M_k(SL₂(ℤ)) modulo mod} in coefficient coordinates (q⁰ … q^(d−1))."""
    d = dim_mk(k)
    if d == 0: return {}
    P = max(primes); L = P * (d - 1) + 1 if d > 1 else P + 1
    mk = lambda c: F6.poly_mod([x % mod for x in c], mod)
    mul = lambda a, b: a.mul_low(b, L)
    s3, s5 = _sigma_list(L, 3), _sigma_list(L, 5)
    E4 = mk([1] + [240 * s3[n] for n in range(1, L)]); E6 = mk([1] + [-504 * s5[n] for n in range(1, L)])
    E1 = F6.eta_series(1, L, mod); Dq = E1.pow_trunc(24, L)            # Δ/q
    w = k - 12 * (d - 1)
    Ew = {0: mk([1]), 4: E4, 6: E6, 8: mul(E4, E4), 10: mul(E4, E6), 14: mul(mul(E4, E4), E6)}[w]
    G = mul(Dq.pow_trunc(d - 1, L), Ew) if d > 1 else Ew                 # G_k / q^(d−1)
    jq = mul(E4.pow_trunc(3, L), Dq.inverse_series_trunc(L))              # j·q
    n = d; B = [0] * (n * n); C = {p: [0] * (n * n) for p in primes}
    pk = {p: pow(p, k - 1, mod) for p in primes}; psi = G
    for i in range(n):                                                    # H_i = q^(d−1−i)·G·(jq)^i; column c = d−1−i
        c = d - 1 - i; lc = psi.length()
        get = lambda t: int(psi[t]) if 0 <= t < lc else 0
        for nn in range(c, n): B[nn * n + c] = get(nn - c) % mod
        for p in primes:
            Cp = C[p]
            for nn in range(n):
                v = get(p * nn - c)
                if nn % p == 0:
                    v2 = get(nn // p - c)
                    if v2: v += pk[p] * v2
                Cp[nn * n + c] = v % mod
        if i < n - 1: psi = mul(psi, jq)
    Bm = F6.mat_mod(n, n, B, mod); I = F6.identity(n, mod); Y = I - Bm; Binv = I + Y; covered = 2
    while covered < n:
        Y = Y * Y; Binv = Binv * (I + Y); covered *= 2
    return {p: F6.mat_mod(n, n, C[p], mod) * Binv for p in primes}


def newton_at(vals, ell, prec, x0):
    """Newton's series through x = 0..X (values mod ℓ^prec), evaluated at the integer x0 (x0 < 0 allowed), and the
    differences' ℓ-adic valuations."""
    M = ell ** prec; X = len(vals) - 1; diffs = []
    for nn in range(X + 1):
        s = 0; b = 1
        for i in range(nn + 1):
            s += (-1) ** (nn - i) * b * vals[i]; b = b * (nn - i) // (i + 1)
        diffs.append(s % M)
    val = []
    for dd in diffs:
        v = 0
        if dd == 0: val.append(prec); continue
        while dd % ell == 0: dd //= ell; v += 1
        val.append(v)
    tot = 0; c = 1                                                         # C(x0, n) = x0(x0−1)…(x0−n+1)/n!
    for nn in range(X + 1):
        tot += c * diffs[nn]; c = c * (x0 - nn) // (nn + 1)
    return tot % M, val
