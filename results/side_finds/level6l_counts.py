# A dig, 14 Nov 2026 (README "Why the waiting starts where it does: the forms of level 6ℓ").
# The forms of level 6ℓ, counted by system beside those of level 6, at the row's weights k = −2 + (ℓ − 1)(x + 1).
#
# THE SPACE M_k(Γ0(6ℓ)) over the ℓ-adic integers.
#   1. A pool of integral forms of level 6ℓ: X(z)·Y(ℓz) with X, Y in forms6_flint's bases of level 6 (weights a + b = k);
#      then the same times one or two factors: E(z) − d·E(dz) (E the weight-2 Eisenstein series, d = ℓ, 2ℓ, 3ℓ, 6ℓ),
#      a(d1·z)·a(d2·z), or the bracket a(d1·z)·θa(d2·z) − θa(d1·z)·a(d2·z) (weight 4), where a = Σ q^(x²+xy+y²) is the
#      weight-1 theta series of level 3 and θ = q·d/dq, d1, d2 among 1, 2, ℓ, 2ℓ. (At weight 4 the Eisenstein products
#      alone fall short: a cusp form meets them only through its central L-value, which is 0 for the forms of sign −1.)
#   2. From it, n = dim M_k(Γ0(6ℓ)) forms independent over Q: their first B + 1 coefficients (B = k(ℓ + 1), Sturm's
#      bound) have rank n modulo a large prime other than ℓ. The script stops if the pool falls short.
#   3. Saturation. A combination of them that is 0 mod ℓ on the first B + 1 coefficients is 0 mod ℓ everywhere (Sturm),
#      so divided by ℓ it is again an integral form, and it replaces one of them. Repeated until the n forms are
#      independent mod ℓ: then they are a basis over Z_ℓ. This is needed: for k ≥ ℓ − 1 the products alone cannot be one,
#      since they reduce mod ℓ into the forms of level 6 and weight kℓ, a space of dimension kℓ + 1 < n. Each round
#      loses one digit. Every coefficient of every combination is checked to be divisible by ℓ.
# THE OPERATORS. Coordinates are the coefficients at n positions ≤ B where the basis is invertible mod ℓ. U_2, U_3 and
#   T_q (q ∤ 6ℓ) are n×n matrices mod ℓ^N. Each is checked on all B + 1 coefficients, and they are checked to commute.
# THE COUNTS, as fading_forms5.py does at level 6. The forms with U_2 = (2|ℓ)^(x+1)·2^((k−2)/2) and
#   U_3 = (3|ℓ)^(x+1)·3^((k−2)/2) exactly (Fitting's splitting, then the exact kernel). They are split by their systems on
#   lane ℓ (T_q, q = 5, 7, 11, 13 other than ℓ). The same at level 6 (forms6_flint.hecke_matrices), with the slopes of T_ℓ
#   on each system. New at ℓ = level 6ℓ − 2·level 6.
# THE CHECKS. The four classes of signs together equal dim S_k(6ℓ) − 2 dim S_k(3ℓ) − 2 dim S_k(2ℓ) + 4 dim S_k(ℓ), and
#   the same at level 6. Every system of level 6 is found at least twice at level 6ℓ.
# COUNTED (the forms), DERIVED (the comparisons).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/level6l_counts.py ℓ X [--from=x]
#              [--N=12] [--M=50] [--keep] [--threads=t]
import os, sys, time, json, itertools
from fractions import Fraction
from math import gcd, prod
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat
import forms6_flint as F


def opt(name, default):
    return next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--' + name + '=')), default)


flint.ctx.threads = int(opt('threads', 4))
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; L, X = nums[0], nums[1]
XS = int(opt('from', 0)); NP = int(opt('N', 12)); MS = int(opt('M', 50)); NP6 = 30
KEEP = '--keep' in sys.argv
PQ = 2 ** 61 - 1                                  # a prime other than ℓ: ranks over Q are read modulo it
QS = [q for q in (5, 7, 11, 13) if q != L]
EDS = (L, 2 * L, 3 * L, 6 * L); TH = (1, 2, L, 2 * L)
FACTORS = [('E', d) for d in EDS] + [('A', d1, d2) for d1 in (1, 2) for d2 in (L, 2 * L)] + \
          [('R', d1, d2) for d1, d2 in itertools.combinations(TH, 2)]
FW = {f: (4 if f[0] == 'R' else 2) for f in FACTORS}
base = os.path.join(HERE, 'level6l_counts_%d_x%d-%d' % (L, XS, X))
OUT, JS = base + '.txt', base + '.jsonl'
if KEEP and (os.path.exists(OUT) or os.path.exists(JS)):
    sys.exit('stopped before computing anything: %s already exists, and results are never overwritten.' % (OUT if os.path.exists(OUT) else JS))
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None
jh = open(JS, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def sign(p, k): return leg(p, L) ** (1 + (k - (L - 3)) // (L - 1))      # the class's sign of U_p (fading_forms5.py)


def vl(x):
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


# ---- dimensions (the standard formula from the genus, the cusps and the elliptic points) ----
def primes_of(N):
    ps, p = [], 2
    while p * p <= N:
        if N % p == 0: ps.append(p)
        while N % p == 0: N //= p
        p += 1
    return ps + ([N] if N > 1 else [])


def dim_S(N, k):
    """(dim S_k(Γ0(N)), number of cusps) for even k ≥ 4; dim M_k = dim S_k + cusps."""
    ps = primes_of(N); mu = Fraction(N)
    for p in ps: mu *= Fraction(p + 1, p)
    k4 = lambda p: 0 if p == 2 else (1 if p % 4 == 1 else -1)               # (−4|p)
    k3 = lambda p: 0 if p == 3 else (1 if p % 3 == 1 else -1)               # (−3|p)
    nu2 = 0 if N % 4 == 0 else prod(1 + k4(p) for p in ps)
    nu3 = 0 if N % 9 == 0 else prod(1 + k3(p) for p in ps)
    c = sum(sum(1 for j in range(1, g + 1) if gcd(j, g) == 1) for d in range(1, N + 1) if N % d == 0 for g in [gcd(d, N // d)])
    g = 1 + mu / 12 - Fraction(nu2, 4) - Fraction(nu3, 3) - Fraction(c, 2)
    S = (k - 1) * (g - 1) + (k // 2 - 1) * c + nu2 * (k // 4) + nu3 * (k // 3)
    assert g.denominator == 1 and S.denominator == 1
    return int(S), c


# ---- q-series ----
def coeffs(f, n):
    c = [int(v) for v in f.coeffs()[:n]]
    return c + [0] * (n - len(c))


class Series:
    """forms6_flint's bases of M_a(Γ0(6)) (q^(a−i)·φ_a·u^i), at z and at ℓz, and E(z) − d·E(dz), mod m to length Lx."""
    def __init__(self, k, Lx, m):
        self.k, self.Lx, self.m = k, Lx, m
        E1, E2, E3, E6 = (F.eta_series(d, Lx, m) for d in (1, 2, 3, 6))
        mul = lambda a, b: a.mul_low(b, Lx)
        u = mul(mul(E1.pow_trunc(5, Lx), E3), mul(E2, E6.pow_trunc(5, Lx)).inverse_series_trunc(Lx))
        g1 = mul(mul(E1, E6.pow_trunc(6, Lx)), mul(E2.pow_trunc(2, Lx), E3.pow_trunc(3, Lx)).inverse_series_trunc(Lx))
        g2 = mul(g1, g1); ph = F.poly_mod([1], m); self.bz = {0: [ph]}
        for a in range(2, k + 1, 2):
            ph = mul(ph, g2); psi = ph; fs = []
            for i in range(a + 1):
                fs.append(mul(psi, F.poly_mod([0] * (a - i) + [1], m)))
                if i < a: psi = mul(psi, u)
            self.bz[a] = fs
        self.bl = {a: [self.at_lz(f) for f in fs] for a, fs in self.bz.items()}
        sg = [0] * Lx; r3 = [0] * Lx
        for d in range(1, Lx):
            ch = 0 if d % 3 == 0 else (1 if d % 3 == 1 else -1)
            for j in range(d, Lx, d): sg[j] += d; r3[j] += ch
        P = lambda c: F.poly_mod([v % m for v in c], m)
        self.fac = {('E', d): P([1 - d] + [-24 * sg[j] + (24 * d * sg[j // d] if j % d == 0 else 0) for j in range(1, Lx)])
                    for d in EDS}
        th0 = [1] + [6 * r3[j] for j in range(1, Lx)]                          # a = Σ q^(x²+xy+y²)
        th = {d: [th0[j // d] if j % d == 0 else 0 for j in range(Lx)] for d in TH}
        dth = {d: [j * v for j, v in enumerate(th[d])] for d in TH}
        for f in FACTORS:
            if f[0] == 'A': self.fac[f] = P(th[f[1]]).mul_low(P(th[f[2]]), Lx)
            if f[0] == 'R': self.fac[f] = P(th[f[1]]).mul_low(P(dth[f[2]]), Lx) - P(dth[f[1]]).mul_low(P(th[f[2]]), Lx)

    def at_lz(self, f):
        c = [0] * self.Lx
        for j, v in enumerate(coeffs(f, (self.Lx - 1) // L + 1)): c[L * j] = v
        return F.poly_mod(c, self.m)

    def form(self, desc):
        fs, a, i, j = desc; b = self.k - sum(FW[f] for f in fs) - a
        f = self.bz[a][i].mul_low(self.bl[b][j], self.Lx)
        for g in fs: f = f.mul_low(self.fac[g], self.Lx)
        return f


def pool(k):
    """The pool's forms, in order: (factors, a, i, j) = factors·X_i(z)·Y_j(ℓz), X of weight a."""
    for e in range(3):
        for fs in itertools.combinations_with_replacement(FACTORS, e):
            w = k - sum(FW[f] for f in fs)
            if w < 0: continue
            for a in range(0, w + 1, 2):
                for i in range(a + 1):
                    for j in range(w - a + 1): yield (fs, a, i, j)


def pivots_of_rref(R, rk, cols):
    piv, c = [], 0
    for i in range(rk):
        while int(R[i, c]) == 0: c += 1
        piv.append(c); c += 1
    return piv


def select_Q(k, n, B):
    """n forms of the pool independent over Q (the first in the pool's order), and how many were tried."""
    S = Series(k, B + 1, PQ); chosen, cols, batch, tried = [], [], [], 0
    def flush():
        allc = cols + [coeffs(S.form(d), B + 1) for d in batch]; descs = chosen + batch
        M = nmod_mat(B + 1, len(allc), [allc[c][t] for t in range(B + 1) for c in range(len(allc))], PQ)
        R, rk = M.rref(); piv = pivots_of_rref(R, rk, len(allc))
        return [descs[c] for c in piv], [allc[c] for c in piv]
    for d in pool(k):
        batch.append(d); tried += 1
        if len(batch) >= 2 * n:
            chosen, cols = flush(); batch = []
            if len(chosen) >= n: break
    if batch and len(chosen) < n: chosen, cols = flush()
    return chosen, tried


def saturate(rows, B, Lx, m):
    """Divide by ℓ every combination that is 0 mod ℓ, until the forms are independent mod ℓ. Returns (rows, rounds,
    how many divisions: the index of the first span in the lattice is ℓ to that power)."""
    rounds = total = 0
    while True:
        n = len(rows)
        Wb = nmod_mat(n, B + 1, [v % L for f in rows for v in coeffs(f, B + 1)], L)
        Kt, nul = Wb.transpose().nullspace()
        if nul == 0: return rows, rounds, total
        Km = nmod_mat(nul, n, [int(Kt[i, c]) for c in range(nul) for i in range(n)], L)
        R, rk = Km.rref(); new = []
        for a in range(rk):
            cv = [int(R[a, i]) for i in range(n)]; piv = next(i for i, v in enumerate(cv) if v); acc = None
            for i, v in enumerate(cv):
                if v: acc = rows[i] * v if acc is None else acc + rows[i] * v
            co = coeffs(acc, Lx)
            bad = sum(1 for v in co if v % L)
            if bad:
                say('   STOP: a combination that is 0 mod ℓ on the first %d coefficients is not 0 mod ℓ at %d of the %d '
                    'coefficients (so the forms are not all of level 6ℓ and weight %d)' % (B + 1, bad, Lx, K_NOW))
                sys.exit(1)
            new.append((piv, F.poly_mod([v // L for v in co], m)))
        for piv, f in new: rows[piv] = f
        rounds += 1; total += rk


def operators(rows, B, k, Np, Lx):
    """U_2, U_3, T_q as n×n matrices mod ℓ^Np in the coordinates (coefficients at the pivot positions); the checks."""
    m = L ** Np; n = len(rows)
    W = [[v % m for v in coeffs(f, Lx)] for f in rows]
    R, rk = nmod_mat(n, B + 1, [W[j][t] % L for j in range(n) for t in range(B + 1)], L).rref()
    assert rk == n
    piv = pivots_of_rref(R, rk, B + 1)
    Pm = F.mat_mod(n, n, [W[j][piv[i]] for i in range(n) for j in range(n)], m); Pinv = F.inverse_mod(Pm, L, Np)
    Wt = F.mat_mod(B + 1, n, [W[j][t] for t in range(B + 1) for j in range(n)], m)
    H, ok = {}, True
    for kind, p in [('U', 2), ('U', 3)] + [('T', q) for q in QS]:
        pk = pow(p, k - 1, m)
        def val(j, t):
            v = W[j][p * t]
            if kind == 'T' and t % p == 0: v = (v + pk * W[j][t // p]) % m
            return v
        Qf = F.mat_mod(B + 1, n, [val(j, t) for t in range(B + 1) for j in range(n)], m)
        Qp = F.mat_mod(n, n, [val(j, piv[i]) for i in range(n) for j in range(n)], m)
        ok = ok and (Wt * (Pinv * Qp) == Qf)
        H[(kind, p)] = Qp * Pinv
    comm = all(H[a] * H[b] == H[b] * H[a] for a, b in itertools.combinations(list(H), 2))
    return H, ok, comm


# ---- the class and its systems (fading_forms5.py's computation) ----
def polygon(chi, Nn):
    """The Newton polygon of chi (lowest first, mod ℓ^Nn): [(slope, length)] smallest first, and whether the smallest
    is exact at this precision."""
    d = len(chi) - 1; m = L ** Nn
    def val(x):
        x %= m
        if x == 0: return Nn
        return vl(x)
    vs = [val(c) for c in chi]; hull = [0]
    for j in range(1, d + 1):
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            if (vs[b] - vs[a]) * (j - a) >= (vs[j] - vs[a]) * (b - a): hull.pop()
            else: break
        hull.append(j)
    sl = [(Fraction(vs[a] - vs[b], b - a), b - a) for a, b in zip(hull, hull[1:])]
    return sl[::-1], vs[hull[-2]] < Nn


def u2_part(H, n, k, s2, Np):
    m = L ** Np; c2 = s2 * 2 ** ((k - 2) // 2) % m
    Y, R, r = F.kernel_part(H[('U', 2)] - F.identity(n, m) * c2, L, Np)
    return Y, R, r, c2


def the_class(H, k, part2, s3, Np, systems=False, slopes=False):
    """(number of forms, precision, {system: (count, slopes of T_ℓ or None)})."""
    Y, R, r, c2 = part2; m = L ** Np; c3 = s3 * 3 ** ((k - 2) // 2) % m
    if r == 0: return 0, Np, {}
    Z, Rz, r2 = F.kernel_part(F.restrict(H[('U', 3)], Y, R, r) - F.identity(r, m) * c3, L, Np)
    if r2 == 0: return 0, Np, {}
    Y, R = Y * Z, [R[j] for j in Rz]; r = r2
    U2W = F.restrict(H[('U', 2)], Y, R, r); U3W = F.restrict(H[('U', 3)], Y, R, r)
    Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)] + \
         [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(r)] for i in range(r)]
    ker, N1 = F.padic_kernel(Bm, 2 * r, r, L, Np); d = len(ker)
    if d == 0 or not systems: return d, N1, {}
    if N1 < 1: return d, N1, None
    m1 = L ** N1
    Xm = F.mat_mod(r, d, [ker[c][i] % m1 for i in range(r) for c in range(d)], m1)
    Rn = F.pivot_rows(nmod_mat(r, d, [ker[c][i] % L for i in range(r) for c in range(d)], L), L)
    XR = F.inverse_mod(F.mat_mod(d, d, [ker[b][Rn[a]] % m1 for a in range(d) for b in range(d)], m1), L, N1)
    def on_new(T):
        TX = F.mat_mod(r, r, [int(v) % m1 for v in T.entries()], m1) * Xm
        return XR * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m1)
    Mq = {q: on_new(F.restrict(H[('T', q)], Y, R, r)) for q in QS}
    Mt = on_new(F.restrict(H[('T', L)], Y, R, r)) if slopes else None
    pieces = [((), F.identity(d, m1), list(range(d)), d)]
    for q in QS:
        newp = []
        for lab, Yp, Rp, rp in pieces:
            Mqr = F.restrict(Mq[q], Yp, Rp, rp)
            for f, ee in F.reduce_mat(Mqr, L).charpoly().factor()[1]:
                co = [int(c) for c in f.coeffs()]
                Bq = Mqr - F.identity(rp, m1) * ((-co[0]) % L) if len(co) == 2 else F.mat_poly(co, Mqr)
                Z, Rz, r2 = F.kernel_part(Bq, L, N1)
                if r2: newp.append((lab + (((-co[0]) % L) if len(co) == 2 else tuple(co),), Yp * Z, [Rp[j] for j in Rz], r2))
        pieces = newp
    out = {}
    for lab, Yp, Rp, rp in pieces:
        sl = None
        if slopes:
            Mtp = F.restrict(Mt, Yp, Rp, rp)
            sl = polygon(F.berkowitz([[int(Mtp[i, j]) for j in range(rp)] for i in range(rp)], m1), N1)
        out[repr(lab)] = (rp, sl)
    return d, N1, out


def fmt_sl(sl):
    if sl is None: return '—'
    segs, exact = sl
    return ', '.join('%s×%d' % (a, b) for a, b in segs) + ('' if exact else ' (smallest NOT exact)')


# ---- the run ----
t00 = time.time(); K_NOW = 0
say('ℓ = %d: the forms of level 6ℓ = %d and of level 6, counted by system, at the weights k = −2 + %d·(x + 1), x = %d..%d'
    % (L, 6 * L, L - 1, XS, X))
say('the class: U_2 = (2|ℓ)^(x+1)·2^((k−2)/2), U_3 = (3|ℓ)^(x+1)·3^((k−2)/2) exactly, with (2|ℓ) = %d, (3|ℓ) = %d; '
    'systems by T_q mod ℓ for q = %s' % (leg(2, L), leg(3, L), QS))
say('precision: the pool mod ℓ^%d before saturation; the counts at level 6ℓ mod ℓ^%d, at level 6 mod ℓ^%d; FLINT %s' % (
    MS, NP, NP6, flint.__version__))
res = {}
for x in range(XS, X + 1):
    k = -2 + (L - 1) * (x + 1); K_NOW = k; t0 = time.time()
    say()
    if k < 4: say('weight %d (x = %d): left out (weight 2 has Eisenstein series with the class\'s numbers)' % (k, x)); continue
    S6l, c6l = dim_S(6 * L, k); n = S6l + c6l; B = k * (L + 1); Lx = max(QS + [2, 3]) * B + 1
    new23 = {'6l': S6l - 2 * dim_S(3 * L, k)[0] - 2 * dim_S(2 * L, k)[0] + 4 * dim_S(L, k)[0],
             '6': dim_S(6, k)[0] - 2 * dim_S(3, k)[0] - 2 * dim_S(2, k)[0] + 4 * dim_S(1, k)[0]}
    chosen, tried = select_Q(k, n, B)
    kinds = {}
    for d in chosen:
        key = '·'.join(sorted(f[0] for f in d[0])) or 'X·Y only'; kinds[key] = kinds.get(key, 0) + 1
    say('weight %d (x = %d): M_k(Γ0(%d)) has dimension %d (formula), Sturm\'s bound %d; %d forms independent over Q '
        'among the first %d of the pool (by their factors: %s)  (%.0f s)' % (
            k, x, 6 * L, n, B, len(chosen), tried, ', '.join('%s %d' % kv for kv in kinds.items()), time.time() - t0))
    if len(chosen) != n:
        say('   STOP: the pool spans %d dimensions over Q, the formula says %d' % (len(chosen), n)); break
    Sr = Series(k, Lx, L ** MS); rows = [Sr.form(d) for d in chosen]; del Sr
    rows, rounds, total = saturate(rows, B, Lx, L ** MS)
    Np = min(NP, MS - rounds)
    if Np < 4:
        say('   STOP: the saturation took %d rounds, leaving ℓ^%d; rerun with a larger --M' % (rounds, MS - rounds)); break
    H, ok, comm = operators(rows, B, k, Np, Lx); del rows
    say('   saturation: %d rounds, %d divisions by ℓ (the products\' span has index ℓ^%d); a basis over Z_ℓ, mod ℓ^%d; '
        'U_2, U_3, T_q right on all %d coefficients: %s; they commute: %s  (%.0f s)' % (
            rounds, total, total, Np, B + 1, 'yes' if ok else 'NO', 'yes' if comm else 'NO', time.time() - t0))
    if not (ok and comm): say('   STOP: the operators fail their checks'); break
    H6 = F.hecke_matrices(k, L ** NP6, [('T', q) for q in QS] + [('T', L), ('U', 2), ('U', 3)])
    cls = (sign(2, k), sign(3, k)); dims = {'6l': {}, '6': {}}; sysd = {'6l': {}, '6': {}}; precs = {}
    for lev, HH, nn, Nn in (('6l', H, n, Np), ('6', H6, k + 1, NP6)):
        for s2 in (1, -1):
            p2 = u2_part(HH, nn, k, s2, Nn)
            for s3 in (1, -1):
                want = (s2, s3) == cls
                d, N1, sy = the_class(HH, k, p2, s3, Nn, systems=want, slopes=(want and lev == '6'))
                dims[lev][(s2, s3)] = d
                if want:
                    precs[lev] = N1
                    if sy is None: say('   STOP: the class at level %s is not exact at this precision' % lev); sys.exit(1)
                    sysd[lev] = sy
    sg = lambda s: '+' if s > 0 else '−'
    for lev, name in (('6l', '6ℓ'), ('6', '6')):
        tot = sum(dims[lev].values())
        say('   forms new at 2 and 3 at level %s, by the signs of U_2, U_3: %s; together %d (formula %d: %s)' % (
            name, ', '.join('(%s,%s) %d' % (sg(a), sg(b), dims[lev][(a, b)]) for a, b in dims[lev]), tot, new23[lev],
            'right' if tot == new23[lev] else 'WRONG'))
    labs = sorted(set(sysd['6l']) | set(sysd['6']), key=lambda s: (-sysd['6'].get(s, (0,))[0], s))
    old_ok = all(sysd['6l'].get(s, (0,))[0] >= 2 * sysd['6'].get(s, (0,))[0] for s in labs)
    say('   the class (%s,%s): %d forms at level 6, %d at level 6ℓ (exact to ℓ^%d and ℓ^%d), in %d systems; each system of '
        'level 6 found at least twice at level 6ℓ: %s' % (sg(cls[0]), sg(cls[1]), dims['6'][cls], dims['6l'][cls],
                                                           precs['6'], precs['6l'], len(labs), 'yes' if old_ok else 'NO'))
    row = {'l': L, 'x': x, 'k': k, 'n': n, 'rounds': rounds, 'divisions': total, 'prec': Np, 'checks': bool(ok and comm),
           'dims6l': {'%d,%d' % a: v for a, v in dims['6l'].items()}, 'dims6': {'%d,%d' % a: v for a, v in dims['6'].items()},
           'new23': new23, 'systems': {}}
    for s in labs:
        d6 = sysd['6'].get(s, (0, None)); d6l = sysd['6l'].get(s, (0, None))
        say('      system %s: level 6: %d, level 6ℓ: %d, new at ℓ: %d;  slopes of T_ℓ at level 6: %s' % (
            s, d6[0], d6l[0], d6l[0] - 2 * d6[0], fmt_sl(d6[1])))
        row['systems'][s] = {'d6': d6[0], 'd6l': d6l[0], 'new': d6l[0] - 2 * d6[0],
                             'slopes6': None if d6[1] is None else [[str(a), b] for a, b in d6[1][0]],
                             'smallest_exact': None if d6[1] is None else d6[1][1]}
    res[x] = row
    if jh: jh.write(json.dumps(row) + '\n'); jh.flush(); os.fsync(jh.fileno())
    say('   (%.0f s for this weight)' % (time.time() - t0))

# ---- the prediction (README, written before the run), then the same reading for every system ----
PRED = {5: [('(2, 2, 1)', "5's law form", [2, 3], [1])],
        7: [('(4, 5, 1)', "7's law form", [1], [0]), ('(5, 6, 0)', "7's quiet form", [2, 3, 4, 5], [0, 1])],
        11: [('(7, 3, 7)', "11's law form", [0, 1], []), ('(0, 4, 1)', "11's quiet form", [1, 2, 3], [0])],
        13: [('(11, 8, 0)', "13's quiet form", [0, 1], [])]}
get = lambda x, s: res[x]['systems'].get(s, {'d6': 0, 'd6l': 0, 'new': 0})
say(); say('THE PREDICTION (README, written before the run): at a waiting step, 0 forms at level 6 and at least 2 new at ℓ; '
           'at an earlier step, at most 1 new at ℓ')
for s, name, wait, early in PRED.get(L, []):
    for x in sorted(wait + early):
        if x not in res: say('   %s %s, x = %d: not computed' % (name, s, x)); continue
        g = get(x, s); w = x in wait
        okp = g['d6'] == 0 and (g['new'] >= 2 if w else g['new'] <= 1)
        say('   %s %s, x = %d (%s): level 6: %d, new at ℓ: %d — %s' % (name, s, x, 'waiting step' if w else 'earlier step',
                                                                      g['d6'], g['new'], 'as predicted' if okp else 'NOT as predicted'))
say(); say('THE SAME READING FOR EVERY SYSTEM of the class at level 6 (the ghost picture for its first form): the steps '
           'before it first appears with at least 2 new at ℓ are its zeros z; its smallest slope at step x is predicted to '
           'be the sum over z of 1 + v_ℓ(x − z)%s' % ('' if XS == 0 else ' (steps before x = %d not computed)' % XS))
agree = tried_n = 0
for s in sorted({s for r in res.values() for s, v in r['systems'].items() if v['d6'] > 0}):
    xs = sorted(res); a = next(x for x in xs if get(x, s)['d6'] > 0)
    zs = [x for x in xs if x < a and get(x, s)['new'] >= 2]
    parts = []
    for x in [x for x in xs if x >= a]:
        v = res[x]['systems'][s]; pred = sum(1 + vl(x - z) for z in zs)
        if v['slopes6'] is None or not v['smallest_exact']: parts.append('x = %d: — (not exact)' % x); continue
        got = Fraction(v['slopes6'][0][0]); tried_n += 1; agree += (got == pred)
        parts.append('x = %d: %s/%s%s' % (x, pred, got, '' if got == pred else ' DIFFERENT'))
    say('   system %s: first at x = %d with %d form%s; its zeros z = %s; predicted/counted smallest slope: %s' % (
        s, a, get(a, s)['d6'], '' if get(a, s)['d6'] == 1 else 's', zs, '; '.join(parts)))
say('   agreeing: %d of %d (system, weight) pairs' % (agree, tried_n))
say(); say('seconds: %.0f' % (time.time() - t00))
