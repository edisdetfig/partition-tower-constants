# A dig, 11 Nov 2026: WHERE THE FRACTAL CONSTANT LEADS. partitions_fractal.txt found that with the laws divided out, two levels
# of the partition numbers differ by a constant on lane ℓ (−1 for 5, −2 for 7, −2 for 11), where Folsom, Kent and Ono's
# theorem (2012) has only 0 for these three primes ("dimension 0", their remark 4). Here:
#   1  The constant to many digits: does it hold beyond the first digit, and does it settle?
#   2  Inside each law: is p(N) divided by the law a Hecke eigenform, like Atkin's class (partitions_cells.txt)?
#   3  Where its multiplier is 0: one more factor of ℓ than the law forces.
#   4  For 5, Atkin's class as level 0 of the same tower.
# Level b: the N with ℓᵇ dividing 24N − 1; Hᵦ(n) = p(ℓᵇn + δ) divided by the power of ℓ that divides every term. p(n) for every
# n up to 2,400,000 by Euler's pentagonal rule, whole numbers only, kept on the lanes 5²¹, 7¹⁸, 11¹⁴, 13¹³. Shares, members,
# failures and constants are COUNTED; formulas FITTED as said; the sources' statements QUOTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_tower.py [--keep]
import os, sys, time, math, itertools, collections
import numpy as np
import sympy as sp
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_tower.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); N = 2_400_000
G, S, k = [], [], 1
while True:
    a, b = k * (3 * k - 1) // 2, k * (3 * k + 1) // 2
    if a > N: break
    s = 1 if k % 2 else -1; G.append(a); S.append(s)
    if b <= N: G.append(b); S.append(s)
    k += 1
o = np.argsort(G); G = np.array(G, dtype=np.int64)[o]; S = np.array(S, dtype=np.int64)[o]
P = {}
for l, e in ((5, 21), (7, 18), (11, 14), (13, 13)):
    L = l ** e; p = np.zeros(N + 1, dtype=np.int64); p[0] = 1; c = 0
    for n in range(1, N + 1):
        while c < len(G) and G[c] <= n: c += 1
        p[n] = int(np.dot(S[:c], p[n - G[:c]])) % L
    P[l] = (p, e)
say('p(n) for every n up to %d on the lanes 5²¹, 7¹⁸, 11¹⁴, 13¹³ (%.0f s); p(100) right on each: %s' % (N, time.time() - t0, all(int(P[l][0][100]) == 190569292 % l ** P[l][1] for l in P)))


def level(l, b):
    p, e = P[l]; A = l ** b; r0 = pow(24, -1, A); Ns = np.arange(r0, N + 1, A, dtype=np.int64); seq = p[Ns]; v = 0
    while v < e and np.all(seq % l ** (v + 1) == 0): v += 1
    return Ns, seq, v, r0


say(); say('1  THE CONSTANT TO MANY DIGITS: Hᵦ₊₂ ≡ λ·Hᵦ (mod ℓʲ) for every n, for the largest j that one λ allows (digits lowest first)')
for l in (5, 7, 11, 13):
    e = P[l][1]; H = {}
    for b in range(1, 12):
        if l ** b > N: break
        Ns, seq, v, r0 = level(l, b)
        if len(seq) < 3: break
        H[b] = (seq // l ** v, v, r0)
    rows = []
    for b in sorted(H):
        if b + 2 not in H: continue
        (a1, v1, r1), (a2, v2, r2) = H[b], H[b + 2]; m = len(a2); a1 = a1[:m]; prec = e - max(v1, v2)
        n0 = int(np.nonzero(a1 % l != 0)[0][0]); best, lam = 0, None
        for j in range(1, prec + 1):
            mod = l ** j; cand = (int(a2[n0]) * pow(int(a1[n0]) % mod, -1, mod)) % mod
            if np.all((a2 - cand * a1) % mod == 0): best, lam = j, cand
            else: break
        ds, t = [], lam
        for _ in range(best): ds.append(t % l); t //= l
        rows.append((b, best, ds, m))
        say('   ℓ = %2d: p(%d^%d·n + %d)/%d^%d against p(%d^%d·n + %d)/%d^%d: holds to %d digit%s (%d terms): %s' % (
            l, l, b + 2, r2, l, v2, l, b, r1, l, v1, best, '' if best == 1 else 's', m, ', '.join(map(str, ds))))
    if l == 5:
        lam12 = sum(d * 5 ** i for i, d in enumerate(rows[-2][2][:12])); m12 = 5 ** 12; B = math.isqrt(m12 // 2)
        r0_, r1_, s0_, s1_ = m12, lam12, 0, 1
        while r1_ > B: q = r0_ // r1_; r0_, r1_, s0_, s1_ = r1_, r0_ - q * r1_, s1_, s0_ - q * s1_
        say('      the smallest fraction that matches the 12 settled digits of 5 is %d/%d: not a simple fraction' % (r1_ * (1 if s1_ > 0 else -1), abs(s1_)))
say('   For 13 (no laws to divide), Folsom, Kent and Ono state A₁₃(2, 4, 2) = 45, i.e. 6 + 3·13: the second row above, digits 6, 3.')


def leg(x, q):
    tab = np.array([0] + [1 if pow(r, (q - 1) // 2, q) == 1 else -1 for r in range(1, q)], dtype=np.int64); return tab[np.asarray(x) % q]


def kron(a, Q): return int(sp.jacobi_symbol(a % Q, Q))


def hecke(l, b, qmax=2000, minmem=40):
    """for level b of ℓ: per prime Q, the multiplier μ with c(N) ≡ μ·c(N*) where 24N − 1 = Q²·(24N* − 1) (Q not dividing the
    second), by the sign of (24N − 1)/Q² at Q; and the residue classes of M' = (24N − 1)/ℓᵇ on lane ℓ that vanish (the square rule)"""
    p, e = P[l]; Ns, seq, v, r0 = level(l, b); A = l ** b; inv = {r: pow(r, -1, l) for r in range(1, l)}
    c = np.zeros(N + 1, dtype=np.int64); c[Ns] = (seq // l ** v) % l; Mp = (24 * Ns - 1) // A; out = {}; bad = []; vanish = set()
    for Q in sp.primerange(5, qmax):
        if Q == l: continue
        s = (Mp % (Q * Q) == 0) & ((Mp // (Q * Q)) % Q != 0)
        if s.sum() < minmem: break
        NN = Ns[s]; m = Mp[s] // (Q * Q); Nst = (A * m + 1) // 24; x = c[NN]; y = c[Nst]; e2 = leg(m, l); EM = leg(m, Q) * kron(A, Q)
        mu = {}; ok = True
        for E2 in (1, -1):
            for sg in (1, -1):
                cell = (e2 == E2) & (EM == sg)
                if cell.sum() < 6: continue
                xx, yy = x[cell], y[cell]; nz = yy != 0
                if not nz.any() and np.all(xx == 0): vanish.add(E2); continue
                rs = set(((xx[nz] * np.array([inv[int(t)] for t in yy[nz]])) % l).tolist())
                if len(rs) == 1 and np.all(xx[~nz] == 0): mu.setdefault(sg, set()).add(rs.pop())
                else: ok = False
        if not ok or any(len(vs) > 1 for vs in mu.values()): bad.append(Q); continue
        out[Q] = {sg: vs.pop() for sg, vs in mu.items()}
    return out, bad, vanish, v


say(); say('2  INSIDE EACH LAW: c(N) = p(N)/ℓᵛ on lane ℓ. Is c(N) ≡ μ·c(N*) when 24N − 1 = Q²·(24N* − 1), one μ per sign of (24N − 1)/Q² at Q?')
CH = {'1': 1, '(−1|Q)': -1, '(2|Q)': 2, '(−2|Q)': -2, '(3|Q)': 3, '(−3|Q)': -3, '(6|Q)': 6, '(−6|Q)': -6}
K = 1200


def eta_quotient(parts, mod):
    shift = sum(d * x for d, x in parts.items()) // 24; f = np.zeros(K + 1, dtype=np.int64); f[0] = 1
    for d, x in parts.items():
        for kk in range(1, K // d + 1):
            st = d * kk
            for _ in range(abs(x)):
                if x > 0: f[st:] = (f[st:] - f[:K + 1 - st]) % mod
                else:
                    for i in range(st, K + 1): f[i] = (f[i] + f[i - st]) % mod
    g = np.zeros(K + 1, dtype=np.int64); g[shift:] = f[:K + 1 - shift]; return g


FORMS = {'η(z)²η(11z)² (curve 11a)': {1: 2, 11: 2}, 'η(z)η(2z)η(7z)η(14z) (14a)': {1: 1, 2: 1, 7: 1, 14: 1}, 'η(z)η(3z)η(5z)η(15z) (15a)': {1: 1, 3: 1, 5: 1, 15: 1},
         'η(2z)²η(10z)² (20a)': {2: 2, 10: 2}, 'η(2z)η(4z)η(6z)η(12z) (24a)': {2: 1, 4: 1, 6: 1, 12: 1}, 'η(3z)²η(9z)² (27a)': {3: 2, 9: 2}, 'η(4z)²η(8z)² (32a)': {4: 2, 8: 2},
         'η(6z)⁴ (36a)': {6: 4}, 'η(8z)⁸/η(4z)²η(16z)² (64a)': {8: 8, 4: -2, 16: -2}, 'η(z)⁴η(5z)⁴': {1: 4, 5: 4}, 'η(2z)⁴η(4z)⁴': {2: 4, 4: 4},
         'η(z)²η(2z)²η(3z)²η(6z)²': {1: 2, 2: 2, 3: 2, 6: 2}, 'η(3z)⁸': {3: 8}, 'η(2z)¹²': {2: 12}, 'η(z)⁶η(3z)⁶': {1: 6, 3: 6}, 'η(z)⁸η(2z)⁸': {1: 8, 2: 8}, 'Δ (Ramanujan\'s τ)': {1: 24}}
ZEROS = {}
for l, levels in ((5, (1, 2, 3)), (7, (1, 2)), (11, (1, 2))):
    for b in levels:
        mu, bad, vanish, v = hecke(l, b); i2 = pow(2, -1, l)
        full = {Q: r for Q, r in mu.items() if 1 in r and -1 in r}
        lam = {Q: ((r[1] + r[-1]) * i2) % l for Q, r in full.items()}; w = {Q: ((r[-1] - r[1]) * i2) % l for Q, r in full.items()}
        wfit = all((kron(-3, Q) * pow(Q, (l - 3) % (l - 1), l) - w[Q]) % l == 0 for Q in w)
        eis = [(cn, a_, cc) for cn, a_, cc in itertools.product(CH, range(l - 1), range(l - 1)) if all((kron(CH[cn], Q) * pow(Q, a_, l) * (1 + pow(Q, cc, l)) - lam[Q]) % l == 0 for Q in lam)]
        zs = sorted((Q, sg) for Q, r in mu.items() for sg, val in r.items() if val == 0); ZEROS[(l, b)] = zs
        say('   ℓ = %2d, level %d (÷%d^%d): one μ for every one of %d primes Q (exceptions: %s);  classes of M\' on lane %d that vanish (square rule): %s' % (
            l, b, l, v, len(mu), bad if bad else 'none', l, sorted(vanish) if vanish else 'none'))
        say('      the half-difference of the two μ is (−3|Q)·Q⁻² for all %d: %s.  The half-sum λ = χ(Q)·Qᵃ·(1 + Q^c)?  %s' % (
            len(w), wfit, ('yes: χ = %s, Q^%d·(1 + Q^%d), as in Atkin\'s class' % eis[0]) if eis else 'no such formula'))
        if b == 1 and l in (7, 11):
            g_hits = []
            for name, parts in FORMS.items():
                g = eta_quotient(parts, l)
                for cn, a_ in itertools.product(CH, range(l - 1)):
                    if all((kron(CH[cn], Q) * pow(Q, a_, l) * int(g[Q]) - lam[Q]) % l == 0 for Q in lam if Q <= K): g_hits.append((name, cn, a_))
            say('      λ = χ(Q)·Qᵃ·(Q-th coefficient of a classic form), tried on %d forms (curves 11a to 64a, Δ, and others): %s' % (len(FORMS), g_hits if g_hits else 'none fits'))
        say('      μ = 0 at: %s' % (', '.join('Q = %d (sign %+d, Q ≡ %d mod %d)' % (Q, sg, Q % l, l) for Q, sg in zs[:9]) or 'none in range'))

say(); say('3  WHERE μ = 0: ONE MORE FACTOR OF ℓ THAN THE LAWS FORCE (progressions; every member up to %d checked)' % N)


def progression(l, b, Q, sign):
    A = l ** b * Q ** 3; out = []
    for B in range(A):
        Mv = 24 * B - 1
        if Mv % l ** b or Mv % (Q * Q) or (Mv // (Q * Q)) % Q == 0: continue
        if kron(Mv // (Q * Q), Q) == sign: out.append(B)
    return A, out


for l, b, Q, sign, need in ((5, 1, 13, -1, 2), (5, 1, 23, 1, 2), (5, 2, 13, -1, 3), (7, 1, 13, 1, 2), (7, 1, 23, -1, 2), (7, 2, 13, 1, 3), (11, 1, 5, -1, 2)):
    assert (Q, sign) in ZEROS[(l, b)]
    p = P[l][0]; A, Bs = progression(l, b, Q, sign); allN = np.concatenate([np.arange(B, N + 1, A) for B in Bs]); bad = int((p[allN] % l ** need != 0).sum())
    say('   p(%d·n + B) ≡ 0 (mod %d^%d), where the laws force %d^%d: %d values of B (%s%s);  %d members, %d not multiples' % (
        A, l, need, l, need - 1, len(Bs), ', '.join(map(str, Bs[:6])), ', …' if len(Bs) > 6 else '', len(allN), bad))
say('   Status: for 5 these are Lovejoy and Ono\'s (Crelle 542, 2002; their example: p(10985N + 1514) ≡ 0 (mod 25), our first row).')
say('   For 7 and 11 the analogous results are in Marie Jameson\'s "A refinement of Ramanujan\'s congruences modulo powers of 7')
say('   and 11", which describes the multipliers with two elliptic curves (her paper would not download, so the curves are not named here).')

say(); say('4  FOR 5, ATKIN\'S CLASS IS LEVEL 0 OF THE TOWER: p(N\')/ℓ^v against p(N), where 24N\' − 1 = ℓ²·(24N − 1) and ℓ does not divide 24N − 1')
for l in (5, 7, 11, 13):
    p, e = P[l]; Ns, seq, v, r0 = level(l, 2); Mp = (24 * Ns - 1) // l ** 2; keep = Mp % l != 0; Nh = Ns[keep]; Mp = Mp[keep]; Nl = (Mp + 1) // 24
    x = (p[Nh] // l ** v) % l; y = p[Nl] % l; inv = {r: pow(r, -1, l) for r in range(1, l)}; cls = leg(-Mp, l); parts = []
    for cc, name in ((-1, 'Atkin\'s class'), (1, 'the other class')):
        s = cls == cc; nz = s & (y != 0); rs = collections.Counter(((x[nz] * np.array([inv[int(t)] for t in y[nz]])) % l).tolist())
        if np.all(x[s] == 0): parts.append('%s: the deeper value is always 0, one more factor of %d (the square rule; %d pairs)' % (name, l, int(s.sum())))
        elif len(rs) == 1 and np.all(x[s & (y == 0)] == 0): parts.append('%s: the deeper value is always %d times (%d pairs)' % (name, next(iter(rs)), int(s.sum())))
        else: parts.append('%s: no single ratio (%d ratios)' % (name, len(rs)))
    say('   ℓ = %2d: %s' % (l, ';  '.join(parts)))
say(); say('seconds: %.0f' % (time.time() - t0))
