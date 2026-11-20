# Forms of level 1 with the multipliers of η⁸ and η¹⁶, for the tower of T_3C = j(3τ)^(1/3) = E₄(3τ)/η(3τ)⁸. 12 Nov 2026.
#
# THE SPACES (derived; checked numerically by the image test in tl_matrix).
#   - η has no zeros in the upper half plane; η⁸ has order 1/3 at ∞ and η¹⁶ order 2/3.
#   - So a form of weight k with η¹⁶'s multiplier is η¹⁶·g, g ∈ M_(k−8)(SL₂(ℤ)) ("class 2"), and with η⁸'s multiplier is
#     η⁸·g, g ∈ M_(k−4) ("class 1").
#   - T_3C has η⁻⁸'s multiplier, which is η¹⁶'s (the multiplier has order 3): class 2.
# SCALED by z → 3z they lie in M_k(Γ0(9)) with trivial character (the eta-quotient theorem: 3¹⁶ and 3⁸ are squares).
# Stored "compressed": c(t) = the coefficient of q^t in E₁⁸·g or E₁¹⁶·g (E₁ = ∏(1 − qⁿ)), the exponent being 3t + r.
# g runs over the basis Δ^(d−1)·E_w·(j·q)^i of forms1_flint, leading index d − 1 − i: triangular coordinates.
# HECKE: (T_ℓF)(n) = A(ℓn) + ℓ^(k−1)·A(n/ℓ) carries class r to class rℓ mod 3.
#   - ℓ ≡ 1 (mod 3): the class is kept, one tower level is T_ℓ.
#   - ℓ ≡ 2 (mod 3): two levels are T_ℓ∘T_ℓ, as for distinct parts.
import forms6_flint as F6
import forms1_flint as F1


def _miller(w, L, mod):
    """The basis of M_w(SL₂(ℤ)) as compressed power series to length L, leading index d − 1 − i, leading coefficient 1."""
    d = F1.dim_mk(w)
    if d == 0: return []
    mk = lambda c: F6.poly_mod([x % mod for x in c], mod)
    mul = lambda a, b: a.mul_low(b, L)
    s3, s5 = F1._sigma_list(L, 3), F1._sigma_list(L, 5)
    E4 = mk([1] + [240 * s3[n] for n in range(1, L)]); E6 = mk([1] + [-504 * s5[n] for n in range(1, L)])
    Dq = F6.eta_series(1, L, mod).pow_trunc(24, L)
    ww = w - 12 * (d - 1)
    Ew = {0: mk([1]), 4: E4, 6: E6, 8: mul(E4, E4), 10: mul(E4, E6), 14: mul(mul(E4, E4), E6)}[ww]
    G = mul(Dq.pow_trunc(d - 1, L), Ew) if d > 1 else Ew
    jq = mul(E4.pow_trunc(3, L), Dq.inverse_series_trunc(L))
    out = []; psi = G
    for i in range(d):
        sh = d - 1 - i
        out.append([0] * sh + [int(psi[t]) % mod for t in range(L - sh)])
        if i < d - 1: psi = mul(psi, jq)
    return out


def basis_compressed(k, r, L, mod):
    """Class r (1: η⁸·M_(k−4); 2: η¹⁶·M_(k−8)) at weight k, compressed, length L."""
    g = _miller(k - 8 if r == 2 else k - 4, L, mod)
    if not g: return []
    Ep = F6.eta_series(1, L, mod).pow_trunc(16 if r == 2 else 8, L)
    out = []
    for c in g:
        s = Ep.mul_low(F6.poly_mod(c, mod), L)
        out.append([int(s[t]) % mod for t in range(L)])
    return out


def tl_matrix(k, r, ell, mod, extra=8):
    """T_ℓ from class r to class r' = rℓ mod 3: a d'×d matrix (columns: images of the class-r basis in the class-r'
    basis), and whether every image lies in the target space on `extra` further coefficients."""
    rp = r * ell % 3
    dt = F1.dim_mk(k - 8 if rp == 2 else k - 4); ds = F1.dim_mk(k - 8 if r == 2 else k - 4)
    if ds == 0 or dt == 0: return [[0] * ds for _ in range(dt)], True
    T = dt + extra; off = (ell * rp - r) // 3
    Lin = ell * (T - 1) + off + 2
    src = basis_compressed(k, r, Lin, mod); tgt = basis_compressed(k, rp, T, mod)
    pk = pow(ell, k - 1, mod); M = [[0] * ds for _ in range(dt)]; ok = True
    for j, c in enumerate(src):
        img = []
        for t in range(T):
            v = c[ell * t + off]; n = 3 * t + rp
            if n % ell == 0:
                q_ = n // ell - r
                if q_ >= 0 and q_ % 3 == 0: v += pk * c[q_ // 3]
            img.append(v % mod)
        res = img[:]; coef = [0] * dt
        for t in range(dt):                                         # the target vector with leading index t is i = dt − 1 − t
            ii = dt - 1 - t; a = res[t] % mod; coef[ii] = a
            if a:
                b = tgt[ii]; res = [(res[x] - a * b[x]) % mod for x in range(T)]
        if any(x % mod for x in res): ok = False
        for i in range(dt): M[i][j] = coef[i]
    return M, ok


def step_operator(k, ell, mod, r=2):
    """One tower level on class r when ℓ ≡ 1 mod 3 (T_ℓ), two levels when ℓ ≡ 2 mod 3 (T_ℓ∘T_ℓ through the other class).
    Returns (matrix as lists, levels per step, images-in-space)."""
    A, ok1 = tl_matrix(k, r, ell, mod)
    if ell % 3 == 1: return A, 1, ok1
    B, ok2 = tl_matrix(k, r * ell % 3, ell, mod)
    d = len(A[0]) if A else 0; dm = len(A)
    C = [[sum(B[i][l] * A[l][j] for l in range(dm)) % mod for j in range(d)] for i in range(d)]
    return C, 2, ok1 and ok2
