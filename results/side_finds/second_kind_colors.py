# A dig, 12 Nov 2026: THEOREM A (ORDINARITY) FOR r COLOURS, BY THE SAME THETA CYCLE. For r colours the class has weight
# k = ℓ − r − 2 (level 6 if 3 ∤ r, level 2 if 3 | r) and U_p = ε·(p|ℓ)·p^((k−2)/2) (p = 2, 3; ε = +1 for r = 1, 3 and −1 for
# r = 5; STATEMENT G1). A non-ordinary class form f gives (derived, NOTE.md §3 with c = r + 2): g = θ^(r+3) f̄ of filtration
# w = r + 5, non-ordinary, with U_p(g) ≡ p^(r+3)·a_p(f) ≡ ε·p^((w−2)/2) (mod ℓ), by Euler's criterion. So:
#   exceptions at ℓ can come only from an eigenform G of weight w = r + 5 on the same level whose U₂ (and U₃) are ≡ those
#   values modulo a prime above ℓ and which is non-ordinary at ℓ; Eisenstein series are ordinary.
# Here: every Hecke eigensystem on S_w (exact, over ℚ): its U₂, U₃ minimal polynomials; whether it matches the required
# values exactly (then its non-ordinary primes are the exceptions) or only modulo finitely many primes (the norms of the
# differences). PREDICTED: r = 1 gives only primes dividing 72 or 18 (none ≥ 5); r = 3 and r = 5 likewise only small
# primes, so all three colour classes are ordinary at every ℓ ≥ r + 6. COUNTED (exact linear algebra; sympy factoring).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/second_kind_colors.py [--keep]
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
from flint import fmpz_poly, fmpz_mat, fmpz
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'second_kind_colors.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def exact_series(Lx):
    def eta(d):
        c = [0] * Lx; c[0] = 1; j = 1
        while d * j * (3 * j - 1) // 2 < Lx:
            s = -1 if j % 2 else 1
            c[d * j * (3 * j - 1) // 2] += s
            if d * j * (3 * j + 1) // 2 < Lx: c[d * j * (3 * j + 1) // 2] += s
            j += 1
        return fmpz_poly(c)
    def peta(d):
        c = [0] * Lx
        for j in range(0, (Lx - 1) // d + 1): c[d * j] = int(fmpz(j).partitions_p())
        return fmpz_poly(c)
    return eta, peta


def U_matrices(k, level):
    """Exact U₂ (and U₃ at level 6) on M_k(Γ0(level)) in coefficient coordinates, as sympy matrices.
    Level 6: basis q^m·φ·u^(k−m) (forms6_flint). Level 2: basis q^m·F^?… built the same way from
    G = (η(2z)¹⁶/η(z)⁸)^(k/4) or F₂·(…)^((k−2)/4) and u = (η(z)/η(2z))²⁴ (forms2_flint)."""
    if level == 6:
        n = k + 1; Lx = 3 * k + 1; eta, peta = exact_series(Lx); mul = lambda a, b: a.mul_low(b, Lx)
        E1, E3, E6 = eta(1), eta(3), eta(6); P2, P3, P6 = peta(2), peta(3), peta(6)
        u = mul(mul(mul(E1.pow_trunc(5, Lx), E3), P2), P6.pow_trunc(5, Lx))
        g1 = mul(mul(mul(E1, E6.pow_trunc(6, Lx)), P2.pow_trunc(2, Lx)), P3.pow_trunc(3, Lx))
        psi = g1.pow_trunc(k, Lx); cols = []
        for i in range(n):
            cols.append((k - i, psi)); psi = mul(psi, u) if i < n - 1 else psi
        ops = (2, 3)
    else:                                                   # level 2: forms2_flint's basis H_j = G_k·u^j, column c = m − j
        import forms2_flint as F2
        m = k // 4; n = m + 1; L = 2 * m + 1
        G, uq, mul = F2._series(k, L, None)
        psi = G; B = sp.zeros(n, n); C2 = sp.zeros(n, n)
        for j in range(n):
            c = m - j; lc = psi.length(); get = lambda i: int(psi[i]) if 0 <= i < lc else 0
            for nn in range(c, n): B[nn, c] = get(nn - c)
            for nn in range(n): C2[nn, c] = get(2 * nn - c)
            if j < n - 1: psi = mul(psi, uq)
        return {2: C2 * B.inv()}
    B = sp.zeros(n, n); C = {p: sp.zeros(n, n) for p in ops}
    for mm, ps in cols:
        lc = ps.length(); get = lambda j: int(ps[j]) if 0 <= j < lc else 0
        for nn in range(mm, n): B[nn, mm] = get(nn - mm)
        for p in ops:
            for nn in range(n): C[p][nn, mm] = get(p * nn - mm)
    Binv = B.inv()
    return {p: C[p] * Binv for p in ops}


def eigensystems(U, level):
    """Joint eigensystems of U₂ (and U₃) on the whole space: list of (minpoly of U₂, minpoly of U₃ on that piece)."""
    x = sp.symbols('x'); out = []
    U2 = U[2]; f2 = sp.factor_list(U2.charpoly(x).as_expr())[1]
    for fac, mult in f2:
        P = sp.Poly(fac, x); Kb = (P.as_expr().subs(x, U2) if False else None)
        M = sp.zeros(*U2.shape); co = P.all_coeffs()
        for c in co: M = M * U2 + c * sp.eye(U2.shape[0])
        ker = M.nullspace()
        if level == 6:
            V = sp.Matrix.hstack(*ker); U3V = (V.T * V).inv() * V.T * U[3] * V
            f3 = sp.factor_list(U3V.charpoly(x).as_expr())[1]
            for g3, m3 in f3: out.append((P.as_expr(), sp.Poly(g3, x).as_expr(), len(ker), m3))
        else:
            out.append((P.as_expr(), None, len(ker), mult))
    return out


x = sp.symbols('x')
for r, level, eps in ((1, 6, 1), (5, 6, -1), (3, 2, 1)):
    w = r + 5; c2 = eps * 2 ** ((w - 2) // 2); c3 = eps * 3 ** ((w - 2) // 2) if level == 6 else None
    say(); say('r = %d colours: class weight ℓ − %d on Γ0(%d); partner weight %d; required U₂ ≡ %+d%s' % (
        r, r + 2, level, w, c2, (', U₃ ≡ %+d' % c3) if c3 else ''))
    U = U_matrices(w, level)
    for p2, p3, dimk, m3 in eigensystems(U, level):
        eis = sp.Poly(p2, x).degree() == 1 and int(-sp.Poly(p2, x).all_coeffs()[-1]) in (1, 2 ** (w - 1))
        n2 = int(sp.resultant(sp.Poly(p2, x), sp.Poly(x - c2, x)))
        n3 = int(sp.resultant(sp.Poly(p3, x), sp.Poly(x - c3, x))) if p3 is not None else 0
        if eis and (level == 2 or (sp.Poly(p3, x).degree() == 1 and int(-sp.Poly(p3, x).all_coeffs()[-1]) in (1, 3 ** (w - 1)))):
            say('   U₂: %s, U₃: %s — Eisenstein (ordinary at every ℓ): excluded' % (p2, p3)); continue
        if n2 == 0 and (n3 == 0 or level == 2):
            say('   U₂: %s, U₃: %s — EXACT MATCH: its non-ordinary primes would be exceptions' % (p2, p3))
        else:
            g = sp.gcd(n2, n3) if level == 6 else n2
            pr = sorted(sp.factorint(abs(g)).keys()) if g else 'all (one difference is 0)'
            say('   U₂: %s, U₃: %s — no exact match; congruent only modulo primes dividing %s → %s' % (
                p2, p3, g, pr))
say()
say('seconds: done')
