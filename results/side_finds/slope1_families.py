# A dig, 12 Nov 2026: THE SMALL SLOPES OF U_ℓ AT ANY WEIGHT, NEGATIVE ONES INCLUDED (Lauder's algorithm on overconvergent
# forms of level 1; A. Lauder, LMS J. Comput. Math. 14 (2011), Algorithm 1). Written for the question "what does 13 need?":
# 1/E₄²'s tower at 13 holds four constants of slope 1 (README "Double poles"). One is the pole's ℓ·u⁻¹²; the other three
# are the U₁₃-values at weight −8 of the slope-1 families of level 1, which no classical weight shows directly.
# METHOD (no tower data used):
#   - Katz's expansion: write k = k₀ + j·(ℓ − 1) with 0 ≤ k₀ < ℓ − 1 and E = E_(ℓ−1). Every overconvergent form of weight k₀
#     is Σ c_a·v_a with v_a = Δ^a·E₄^x·E₆^y / E^i(a) (4x + 6y = k₀ + i(a)(ℓ − 1) − 12a), i(a) the least i with
#     a < dim M_(k₀ + i(ℓ−1)): one basis vector per q-order a, a "Katz degree" i(a) for each.
#   - Weight k: U_ℓ at weight k equals h ↦ U(e^j·h) at weight k₀, e = E(q)/E(q^ℓ) (Coleman); j < 0 is allowed.
#   - The matrix C of that operator in the basis v_a (by triangular elimination on q-orders), then A = C scaled by
#     ℓ^(⌊i(c)/(ℓ+1)⌋ − ⌊i(r)/(ℓ+1)⌋). Dwork's bound v(C[r][c]) ≥ (ℓ·i(r) − i(c))/(ℓ + 1) makes A's rows of Katz degree i
#     divisible by ℓ^(about i(ℓ−1)/(ℓ+1)), so det(1 − X·A) on the first n vectors is the characteristic series modulo ℓ^N
#     once n·(ℓ−1)/(ℓ+1) > N + 1. The bound is checked on every matrix, not assumed.
#   - The slope-0 factor ∏(z − λ) and the slope-1 factor ∏(y − λ/ℓ) of the series, by Hensel's lemma on reversed
#     polynomials. The slope-1 step assumes the next slopes are ≥ 2; the Newton polygon is printed to check it.
# CHECKS built in: at a classical weight (k ≥ 2) the same factors from T_ℓ on M_k(SL₂(ℤ)) (forms1_flint); the slope-0
# factor at a negative weight must be the Eisenstein family's x − 1 when forms1_family finds no ordinary cusp family.
# COUNTED (exact integer arithmetic modulo ℓ^N throughout).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/slope1_families.py ℓ k N [--n=…] [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import fmpz_mat, nmod_poly
from sympy import bernoulli
import forms6_flint as F
import forms1_flint as F1
flint.ctx.threads = 4
nums = [a for a in sys.argv[1:] if a.lstrip('-').isdigit()]; L, K, N = (int(a) for a in nums[:3])
NB = next((int(a[4:]) for a in sys.argv[1:] if a.startswith('--n=')), None)
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'slope1_families_%d_w%d_N%d.txt' % (L, K, N)), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, m_exp):
    x %= L ** m_exp
    if x == 0: return m_exp
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def digits(x, k): return [(x // L ** e) % L for e in range(k)]


def eisenstein(w, n, mod):
    """E_w = 1 − (2w/B_w)·Σ σ_(w−1)(m) q^m to n terms modulo mod (the constant needs B_w's numerator prime to ℓ)."""
    B = bernoulli(w); c = -2 * w / B; num, den = int(c.p), int(c.q)
    s = [0] * n
    for d in range(1, n):
        dp = pow(d, w - 1, mod)
        for m_ in range(d, n, d): s[m_] = (s[m_] + dp) % mod
    k_ = num * pow(den, -1, mod) % mod
    return F.poly_mod([1] + [k_ * s[m_] % mod for m_ in range(1, n)], mod)


def factor_slope(c, s, h, prec_in, Tcut):
    """c: the series' coefficients c_t (t ≤ Tcut) known modulo ℓ^prec_in, with c_0 = 1, all its roots X = 1/λ with v(λ) ≥ s
    and exactly h of them with v(λ) = s. Returns the monic ∏(z − λ/ℓ^s) over those h, modulo ℓ^(prec_in − Tcut), using the
    reversed truncated polynomial z^T·c(z^(−1)·ℓ^(−s)) (a monic polynomial since c_0 = 1)."""
    M = prec_in - Tcut
    if M <= 0: return None, 0
    mm = L ** M
    r = []
    for t in range(Tcut + 1):
        x = c[t] % L ** prec_in
        if s * t:
            if x % L ** (s * t): raise ValueError('coefficient %d not divisible by ℓ^%d' % (t, s * t))
            x //= L ** (s * t)
        r.append(x % mm)
    if any(r[t] % L for t in range(h + 1, Tcut + 1)) or r[h] % L == 0: return None, 0
    fstar = r[::-1]                                         # z^T + r_1 z^(T−1) + … + r_T, lowest first
    g = [r[h - e] % L for e in range(h + 1)]                # ḡ(z) = r̄_h + r̄_(h−1) z + … + r̄_0 z^h, monic
    cof = [0] * (Tcut - h) + [1]
    Gf, _ = F.hensel_factor(fstar, g, cof, L, M)
    return [x % mm for x in Gf], M


def unit_roots(G, M):
    """the roots in ℤ_ℓ of the monic G (lowest first) that are simple on lane ℓ, lifted to ℓ^M."""
    mm = L ** M; out = []
    for y0 in range(L):
        fy = sum(c * pow(y0, e, L) for e, c in enumerate(G)) % L
        dfy = sum(e * c * pow(y0, e - 1, L) for e, c in enumerate(G) if e) % L
        if fy == 0 and dfy:
            y = y0
            for _ in range(2 * M + 4):
                fy = sum(c * pow(y, e, mm) for e, c in enumerate(G)) % mm
                dfy = sum(e * c * pow(y, e - 1, mm) for e, c in enumerate(G) if e) % mm
                y = (y - fy * pow(dfy, -1, mm)) % mm
            out.append(y)
    return out


t0 = time.time()
k0 = K % (L - 1); jw = (K - k0) // (L - 1)
dims = []; i = 0
NB = NB or int((N + 3) * (L + 1) / (L - 1)) + 4                    # Katz degree reached must exceed (N + 1)(ℓ + 1)/(ℓ − 1)
while True:                                                          # Katz degrees i(a) of the first NB q-orders
    dims.append(F1.dim_mk(k0 + i * (L - 1)))
    if dims[-1] >= NB and i * (L - 1) / (L + 1) > N + 2: break
    i += 1
deg = []
for a in range(10 ** 6):
    ia = next(ii for ii, d in enumerate(dims) if a < d)
    if ia * (L - 1) / (L + 1) > N + 2: break
    deg.append(ia)
NB = len(deg); top_i = deg[-1]
NP = N + top_i // (L + 1) + 3                                        # working precision: the scaling divides by ≤ ℓ^(top_i/(ℓ+1))
mod = L ** NP; LQ = L * NB + 1
say('ℓ = %d, weight %d = %d + %d·(ℓ − 1): %d basis vectors (Katz degrees 0..%d), working modulo ℓ^%d, the series wanted '
    'modulo ℓ^%d; FLINT %s' % (L, K, k0, jw, NB, top_i, NP, N, flint.__version__))
E4, E6, E = eisenstein(4, LQ, mod), eisenstein(6, LQ, mod), eisenstein(L - 1, LQ, mod)
Dq = F.eta_series(1, LQ, mod).pow_trunc(24, LQ)                     # Δ/q
Einv = E.inverse_series_trunc(LQ)
EV = F.poly_mod([int(E[m_ // L]) if m_ % L == 0 else 0 for m_ in range(LQ)], mod)   # E(q^ℓ)
e = E.mul_low(EV.inverse_series_trunc(LQ), LQ)                       # E(q)/E(q^ℓ)
ej = e.pow_trunc(jw, LQ) if jw >= 0 else e.inverse_series_trunc(LQ).pow_trunc(-jw, LQ)
basis = []
for a in range(NB):
    w2 = k0 + deg[a] * (L - 1) - 12 * a
    x, y = (w2 // 4, 0) if w2 % 4 == 0 else ((w2 - 6) // 4, 1)
    assert w2 >= 0 and 4 * x + 6 * y == w2 and x >= 0
    v = Dq.pow_trunc(a, LQ).mul_low(E4.pow_trunc(x, LQ), LQ)
    if y: v = v.mul_low(E6, LQ)
    if deg[a]: v = v.mul_low(Einv.pow_trunc(deg[a], LQ), LQ)
    basis.append([0] * a + [int(v[m_]) for m_ in range(NB - a)])     # times q^a; coefficients q^0 … q^(NB−1)
C = [[0] * NB for _ in range(NB)]
for c in range(NB):
    v = basis[c]
    f = F.poly_mod([0] * c + [int(x) for x in (Dq.pow_trunc(c, LQ).mul_low(E4.pow_trunc((k0 + deg[c] * (L - 1) - 12 * c) // 4 if (k0 + deg[c] * (L - 1) - 12 * c) % 4 == 0 else (k0 + deg[c] * (L - 1) - 12 * c - 6) // 4, LQ), LQ).mul_low(E6 if (k0 + deg[c] * (L - 1) - 12 * c) % 4 else F.poly_mod([1], mod), LQ).mul_low(Einv.pow_trunc(deg[c], LQ) if deg[c] else F.poly_mod([1], mod), LQ).coeffs())][:LQ - c], mod)
    g = f.mul_low(ej, LQ)
    u = [int(g[L * m_]) if L * m_ < LQ else 0 for m_ in range(NB)]   # U(e^j·v_c), q^0 … q^(NB−1)
    for r in range(NB):                                              # triangular elimination on q-orders
        co = u[r] % mod; C[r][c] = co
        if co:
            br = basis[r]
            for m_ in range(r, NB): u[m_] = (u[m_] - co * br[m_]) % mod
worst = min(((val(C[r][c], NP) if C[r][c] % mod else NP) - (L * deg[r] - deg[c]) / (L + 1)) for r in range(NB) for c in range(NB))
say('   Dwork\'s bound v(C[r][c]) − (ℓ·i(r) − i(c))/(ℓ+1): smallest value %.3f (the scaling below needs ≥ −1)' % worst)
A = [[0] * NB for _ in range(NB)]
for r in range(NB):
    for c in range(NB):
        sh = deg[c] // (L + 1) - deg[r] // (L + 1)
        x = C[r][c] % mod
        if sh >= 0: A[r][c] = x * L ** sh % mod
        else:
            if x % L ** (-sh): raise ValueError('entry (%d, %d) not divisible by ℓ^%d' % (r, c, -sh))
            A[r][c] = x // L ** (-sh)
rowdecay = [min(val(A[r][c], NP) for c in range(NB)) for r in range(NB)]
say('   the scaled matrix: smallest valuation of each row (by Katz degree), the last rows: %s' % rowdecay[-6:])
chi = fmpz_mat(NB, NB, [A[r][c] for r in range(NB) for c in range(NB)]).charpoly().coeffs()   # det(x − A), lowest first
mN = L ** N
cser = [int(chi[NB - t]) * 1 % mN for t in range(NB + 1)]           # det(1 − X·A) = Σ c_t X^t, c_t = coefficient of x^(n−t)
vs = [val(x, N) for x in cser]
say('   det(1 − X·U) modulo ℓ^%d: valuations of c_0 … c_%d: %s' % (N, min(NB, 16), vs[:17]))
hull = [0]
for t in [t for t in range(1, min(NB, 40) + 1) if vs[t] < N]:      # only the coefficients known to be nonzero
    while len(hull) >= 2 and (vs[hull[-1]] - vs[hull[-2]]) * (t - hull[-2]) >= (vs[t] - vs[hull[-2]]) * (hull[-1] - hull[-2]): hull.pop()
    hull.append(t)
segs = [((vs[b] - vs[a]) / (b - a), b - a) for a, b in zip(hull, hull[1:]) if vs[b] < N]
say('   its Newton polygon (slope, length), while below ℓ^%d: %s' % (N, [('%.3f' % s_, ln) for s_, ln in segs]))
h0 = max([t for t in range(NB + 1) if vs[t] == 0])
h1 = 0                                                               # the slope-1 run right after: v(c_(h0+e)) = e, e = 1, 2, …
while h0 + h1 + 1 <= NB and vs[h0 + h1 + 1] == h1 + 1 < N: h1 += 1
say('   slope 0: %d roots;  slope 1: %d roots' % (h0, h1))
# slope 0: the unit eigenvalues λ (roots z = λ of the reversed series)
T0 = min(NB, N + h0 - 2)
fstar = [cser[t] % mN for t in range(T0 + 1)][::-1]
g0 = [cser[h0 - e_] % L for e_ in range(h0 + 1)]
G0, _ = F.hensel_factor(fstar, g0, [0] * (T0 - h0) + [1], L, N)
G0 = [x % mN for x in G0]
say('   the slope-0 factor ∏(z − λ) modulo ℓ^%d, lowest first: %s' % (N, G0))
for lam in unit_roots(G0, N): say('      λ = %s (digits %s)' % (lam, digits(lam, min(N, 12))))
# divide the series by ∏(1 − λX) = Σ G0[s]·X^(h0 − s)
den = [G0[h0 - t] for t in range(h0 + 1)]
q = []
for t in range(NB + 1):
    x = cser[t] - sum(den[s_] * q[t - s_] for s_ in range(1, min(t, h0) + 1))
    q.append(x % mN)
# slope 1 on the rest: r_t = q_t/ℓ^t, known modulo ℓ^(N − t). Truncating at T drops the terms t > T, so the factor is right to
# min(N − T, the least v(r_t) over t > T). v(q_t) is read where q_t ≢ 0 (mod ℓ^N); beyond, it is bounded by the Newton
# polygon's last slope (its slopes only grow) and by N. T is chosen to make that precision largest.
vq = [val(x, N) for x in q]
kn = [t for t in range(len(q)) if vq[t] < N]
tl = kn[-1]; sl = max(((vq[tl] - vq[t]) / (tl - t) for t in kn if t < tl), default=1)   # the slope into the last known point
lb = [vq[t] if vq[t] < N else max(N, vq[tl] + sl * (t - tl)) for t in range(len(q))]
say('   after the slope-0 factor: v(q_t) for t = 0 … %d: %s (the last known at t = %d, slope into it %.2f)' % (
    min(len(q) - 1, 14), vq[:15], tl, sl))
cands = []
for Tc in range(h1 + 1, min(N, len(q))):
    drop = min((lb[t] - t for t in range(Tc + 1, len(q))), default=N)
    cands.append((min(N - Tc, int(drop)), Tc))
Mbest, Tc = max(cands)
G1, M1 = factor_slope(q, 1, h1, N, Tc)
M1 = min(M1, Mbest); G1 = [x % L ** M1 for x in G1]
say('   the slope-1 factor ∏(y − λ/ℓ) modulo ℓ^%d (truncated at t = %d), lowest first: %s' % (M1, Tc, G1))
say('      on lane ℓ: %s; its roots on lane ℓ: %s' % ([x % L for x in G1], [int(r_) for r_, m_ in nmod_poly([x % L for x in G1], L).roots() for _ in range(m_)]))
for y in unit_roots(G1, M1): say('      y = λ/ℓ, digits lowest first: %s' % digits(y, M1))
if K >= 2:                                                          # the classical check: T_ℓ on M_K(SL₂(ℤ))
    d = F1.dim_mk(K); NC = N + d + 6; mC = L ** NC                  # the whole polynomial, no truncation: enough digits for r_d
    Tm = F1.hecke_matrices1(K, mC, [L])[L]
    ch = fmpz_mat(d, d, [int(Tm[r, c]) for r in range(d) for c in range(d)]).charpoly().coeffs()
    cs = [int(ch[d - t]) % mC for t in range(d + 1)]
    say('   CHECK, T_ℓ on M_%d(SL₂(ℤ)) (dimension %d), the whole polynomial modulo ℓ^%d: valuations of det(1 − X·T_ℓ): %s' % (
        K, d, NC, [val(x, NC) for x in cs]))
    fs = cs[::-1]
    G0c, _ = F.hensel_factor(fs, [cs[h0 - e_] % L for e_ in range(h0 + 1)], [0] * (d - h0) + [1], L, NC)
    denc = [G0c[h0 - t] % mC for t in range(h0 + 1)]; qc = []
    for t in range(d + 1): qc.append((cs[t] - sum(denc[s_] * qc[t - s_] for s_ in range(1, min(t, h0) + 1))) % mC)
    G1c, M1c = factor_slope(qc, 1, h1, NC, d)
    agree = min(next((e_ for e_ in range(min(M1, M1c)) if (a_ - b_) % L ** (e_ + 1)), min(M1, M1c)) for a_, b_ in zip(G1, G1c)) if G1c else 0
    say('      classical slope-1 factor modulo ℓ^%d: %s; agrees with the overconvergent one to %d digits' % (M1c, G1c, agree))
say('seconds: %.0f' % (time.time() - t0))
