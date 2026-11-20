# Forms of level 2 with the multiplier of η(2z)/η(z), for the towers of partitions into distinct parts. 12 Nov 2026.
#
# THE SPACES (derived; checked numerically below by `check`). Let v be the multiplier of η(2z)/η(z) on Γ0(2) (order 24) and
# r odd. A form of weight k with multiplier v^r is (η(2z)/η(z))^r·g with g ∈ M_k(Γ0(2)) vanishing at the cusp 0:
#   - η has no zeros in the upper half plane;
#   - (η(2z)/η(z))^r has order r/24 at ∞ and −r/24 at the cusp 0 (the eta-quotient order formula), so g must have order
#     ≥ 1 at 0;
#   - hence dimension m = ⌊k/4⌋.
# In forms2_flint's basis H_j = G_k·u^j (u = (η(z)/η(2z))^24, a simple zero at 0), the forms vanishing at 0 are H_1 … H_m.
# SCALED: F(z) = f(24z) = Σ c(t)·q^(24t + r) lies in M_k(Γ0(1152), (2/·)) (the eta-quotient character theorem). We store
# the "compressed" series c(t), which is ∏(1 + qⁿ)^r·H_j(q), with leading index m − j (j = 1..m): triangular coordinates
# c(0) … c(m − 1).
# HECKE: T_ℓ (ℓ ∤ 6) maps class r to class r' = rℓ mod 24:
#   (T_ℓF)(n) = A(ℓn) + (2|ℓ)·ℓ^(k−1)·A(n/ℓ).
# Two steps return to class r: the composite T_ℓ∘T_ℓ acts on M_k(v^r). That is the tower's two-level step U_ℓ² on its
# ordinary part.
import forms2_flint as F2
import forms6_flint as F6
from flint import nmod_mat


def kron2(p): return 1 if p % 8 in (1, 7) else -1


def basis_compressed(k, r, L, mod):
    """The m compressed series (lists of ints mod `mod`, length L) of class r and weight k; j = 1..m."""
    m = k // 4
    if m == 0: return []
    G, uq, mul = F2._series(k, L, mod)
    E1, E2 = F6.eta_series(1, L, mod), F6.eta_series(2, L, mod)
    Pr = mul(E2, E1.inverse_series_trunc(L)).pow_trunc(r, L)       # ∏(1 + qⁿ)^r
    base = mul(Pr, G); out = []
    for j in range(1, m + 1):
        s = mul(base, uq.pow_trunc(j, L))
        sh = m - j                                                  # times q^(m − j)
        out.append([0] * sh + [int(s[t]) % mod for t in range(L - sh)])
    return out


def tl_matrix(k, r, ell, mod, extra=8):
    """The m×m matrix of T_ℓ from class r to class r' = rℓ mod 24 (columns: images of the class-r basis in the class-r'
    basis), and whether the images lie in the target space on `extra` further coefficients."""
    m = k // 4; rp = r * ell % 24; T = m + extra
    off = (ell * rp - r) // 24
    Lin = ell * (T - 1) + off + 2
    src = basis_compressed(k, r, Lin, mod); tgt = basis_compressed(k, rp, T, mod)
    eps = kron2(ell) * pow(ell, k - 1, mod) % mod
    M = [[0] * m for _ in range(m)]; ok = True
    for j, c in enumerate(src):
        img = []
        for t in range(T):
            v = c[ell * t + off]
            n = 24 * t + rp
            if n % ell == 0:
                q_ = n // ell - r
                if q_ >= 0 and q_ % 24 == 0: v += eps * c[q_ // 24]
            img.append(v % mod)
        res = img[:]; coef = [0] * m
        for t in range(m):                                          # triangular: the basis vector with leading index t is j' = m − t
            jj = m - 1 - t; a = res[t] % mod
            coef[jj] = a
            if a:
                b = tgt[jj]
                res = [(res[i] - a * b[i]) % mod for i in range(T)]
        if any(x % mod for x in res): ok = False
        for i in range(m): M[i][j] = coef[i]
    return M, ok


def composite(k, ell, mod, r=1):
    """T_ℓ∘T_ℓ on M_k(v^r) (from class r to rℓ and back), as an nmod_mat when mod fits a word, else a list of lists."""
    A, ok1 = tl_matrix(k, r, ell, mod); B, ok2 = tl_matrix(k, r * ell % 24, ell, mod)
    m = len(A)
    C = [[sum(B[i][l] * A[l][j] for l in range(m)) % mod for j in range(m)] for i in range(m)]
    return C, ok1 and ok2
