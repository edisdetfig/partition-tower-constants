# A dig, 11 Nov 2026: THE RULE OF EACH PRIME, AS A MODULAR FORM, AND WHICH ONES ARE CURVES. Ahlgren, Allen and Tang
# ("Congruences like Atkin's for the partition function", arXiv 2112.09481, Prop. 3.3 and around it; QUOTED) work with the
# newforms of weight ℓ − 3 and level 6: S₂(6) = {0}; S₄^new(6) is spanned by q − 2q² − 3q³ + 4q⁴ + 6q⁵ + …; S₁₀^new(6) by
# q − 16q² + 81q³ + 256q⁴ + 2694q⁵ + …; S₁₆^new(6) holds three newforms. Here, for every prime ℓ from 7 to 31:
#   1  those newforms, computed from scratch (products of E₂(z) − d·E₂(dz), d = 2, 3, 6, which give all forms of weight k on
#      Γ0(6); newforms picked out by their U₂ and U₃ eigenvalues ±2^((k−2)/2), ±3^((k−2)/2) and split by T₅);
#   2  which one is the partition tower's rule: λ(Q) ≡ ±(d|Q)·Q^b·a_Q(f) (mod ℓ) at every prime Q where our data give the step
#      rule's eigenvalue λ (levels 1 and 2, as in partitions_curves.txt);
#   3  which elliptic curves have the same numbers on lane ℓ: a_Q(f) ≡ ±(d|Q)·Q^b·a_Q(E) (mod ℓ) at every prime Q up to 97 that
#      does not divide 6ℓ·N(E), over ALL isogeny classes of Cremona's tables with conductor below 70,000 (a curve may have
#      extra bad primes); then every survivor rechecked at all primes up to 389.
# Coefficients and point counts COUNTED; the matches are counts. Cremona's tables: the folder given as the argument.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_newforms.py <folder> [--keep]
import os, re, sys, time, collections
import numpy as np
import sympy as sp
from sympy import QQ
from sympy.polys.matrices import DomainMatrix
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv; SRC = [a for a in sys.argv[1:] if a != '--keep'][0]
fh = open(os.path.join(HERE, 'partitions_newforms.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); LEN = 400; PRIMES = (7, 11, 13, 17, 19, 23, 29, 31)
leg = lambda a, q: int(sp.jacobi_symbol(a % q, q))

# ---------- 1. the newforms of weight ℓ − 3 and level 6 ----------
E2 = [1] + [-24 * int(sp.divisor_sigma(n)) for n in range(1, LEN)]
def scale(f, d): g = [0] * LEN; g[::d] = f[:(LEN - 1) // d + 1]; return g
GENS = [[a - d * b for a, b in zip(E2, scale(E2, d))] for d in (2, 3, 6)]
def mul(f, g): return list(np.convolve(np.array(f, dtype=object), np.array(g, dtype=object))[:LEN])


def nullspace(rows):
    """basis of the null space, each vector scaled to whole numbers"""
    ns = DomainMatrix([[QQ(int(x)) for x in r] for r in rows], (len(rows), len(rows[0])), QQ).nullspace().to_Matrix(); out = []
    for i in range(ns.rows):
        v = [sp.Rational(x) for x in ns.row(i)]; den = sp.ilcm(*[x.q for x in v]) if len(v) > 1 else v[0].q
        out.append([int(x * den) for x in v])
    return out


def newforms(k):
    cur = {(0, 0, 0): [1] + [0] * (LEN - 1)}
    for _ in range(k // 2):
        nxt = {}
        for key, f in cur.items():
            for i in range(3):
                nk = list(key); nk[i] += 1; nk = tuple(nk)
                if nk not in nxt: nxt[nk] = mul(f, GENS[i])
        cur = nxt
    P = 2 ** 61 - 1; chosen, red = [], []                      # independent monomials, picked by elimination mod P
    for f in cur.values():
        v = [x % P for x in f[:3 * (k + 1) + 20]]
        for piv, r in red:
            if v[piv]: c = v[piv] * pow(r[piv], -1, P) % P; v = [(a - c * b) % P for a, b in zip(v, r)]
        nz = next((i for i, x in enumerate(v) if x), None)
        if nz is not None: red.append((nz, v)); chosen.append(f)
        if len(chosen) == k + 1: break
    h = (k - 2) // 2; out = []
    for e2 in (1, -1):
        for e3 in (1, -1):
            ns = nullspace([[b[2 * n] - e2 * 2 ** h * b[n] for b in chosen] for n in range(min(LEN // 2, 4 * k + 20))])
            if not ns: continue
            B2 = [[sum(v[i] * chosen[i][n] for i in range(len(chosen))) for n in range(LEN)] for v in ns]
            ns3 = nullspace([[b[3 * n] - e3 * 3 ** h * b[n] for b in B2] for n in range(min(LEN // 3, 4 * k + 20))])
            if not ns3: continue
            forms = [[sum(v[i] * B2[i][n] for i in range(len(B2))) for n in range(LEN)] for v in ns3]
            if len(forms) == 1:
                f = forms[0]; c = f[1]; out.append((e2, e3, [sp.Rational(x, 1) / c for x in f]))
            else:
                ns5 = None                                     # split by T5 when the eigenvalue is a whole number
                T5 = [[(b[5 * n] if 5 * n < LEN else None) for n in range(LEN // 5)] for b in forms]
                out.append((e2, e3, '%d-dimensional' % len(forms)))
    return out


say('1  THE NEWFORMS OF WEIGHT ℓ − 3 AND LEVEL 6 (coefficients up to q^%d; a₂ = −ε₂·2^((k−2)/2), a₃ = −ε₃·3^((k−2)/2))' % (LEN - 1))
NFW = {}
for l in PRIMES:
    k = l - 3; t = time.time(); res = newforms(k); NFW[k] = []
    for e2, e3, f in res:
        if isinstance(f, str): say('   ℓ = %2d, weight %2d, a₂ = %+d·2^%d, a₃ = %+d·3^%d: %s (not split here)' % (l, k, -e2, (k - 2) // 2, -e3, (k - 2) // 2, f)); continue
        assert all(x.q == 1 for x in f[:LEN]), 'not integral'
        a = [int(x) for x in f]; NFW[k].append(a)
        say('   ℓ = %2d, weight %2d: q %+d·q² %+d·q³ %+d·q⁴ %+d·q⁵ %+d·q⁷ …  (%.0f s)' % (l, k, a[2], a[3], a[4], a[5], a[7], time.time() - t))

# ---------- 2. which newform is the tower's rule ----------
NMAX = 2_400_000
G, S, kk = [], [], 1
while True:
    x1, x2 = kk * (3 * kk - 1) // 2, kk * (3 * kk + 1) // 2
    if x1 > NMAX: break
    s = 1 if kk % 2 else -1; G.append(x1); S.append(s)
    if x2 <= NMAX: G.append(x2); S.append(s)
    kk += 1
o = np.argsort(G); G = np.array(G, dtype=np.int64)[o]; S = np.array(S, dtype=np.int64)[o]; PP = []
for L in (7 ** 3 * 11 ** 3 * 13, 17 * 19 * 23 * 29 * 31):
    p = np.zeros(NMAX + 1, dtype=np.int64); p[0] = 1; c = 0
    for n in range(1, NMAX + 1):
        while c < len(G) and G[c] <= n: c += 1
        p[n] = int(np.dot(S[:c], p[n - G[:c]])) % L
    PP.append(p)
LAW = lambda l, b: b if l == 11 else (b // 2 + 1 if l == 7 else 0)


def vl(x, l):
    v = 0
    while x % l == 0: x //= l; v += 1
    return v


def eigen(l, b):
    p = PP[0] if l <= 13 else PP[1]; A = l ** b; r0 = pow(24, -1, A); lam = {}
    for Q in sp.primerange(5, 114):
        if Q == l: continue
        g = collections.defaultdict(set); cnt = collections.Counter()
        for Ns in range(r0, (NMAX + (Q * Q - 1) // 24) // (Q * Q) + 1, A):
            M = 24 * Ns - 1
            if vl(M, l) != b or M % Q == 0 or (l == 7 and b % 2 == 1 and leg(M // A, 7) == 1): continue
            x = (int(p[Ns]) // l ** LAW(l, b)) % l
            if x == 0: continue
            y = (int(p[Q * Q * Ns - (Q * Q - 1) // 24]) // l ** LAW(l, b)) % l; e = leg(M, Q)
            g[e].add(y * pow(x, -1, l) % l); cnt[e] += 1
        if len(g) == 2 and min(cnt.values()) >= 2 and all(len(v) == 1 for v in g.values()):
            lam[Q] = (next(iter(g[1])) + next(iter(g[-1]))) * pow(2, -1, l) % l
    return lam


DS = (1, -1, 2, -2, 3, -3, 6, -6)
say(); say('2  WHICH NEWFORM IS THE TOWER\'S RULE: λ(Q) ≡ ±(d|Q)·Q^b·a_Q(f) (mod ℓ) at every prime Q with data (levels 1 and 2). On lane ℓ')
say('   weights that differ by a multiple of ℓ − 1 behave alike, so every computed weight k ≡ ℓ − 3 (mod ℓ − 1) is tried')
TOWER = {}
for l in PRIMES:
    cands = [(k, i, a) for k in sorted(NFW) if (k - (l - 3)) % (l - 1) == 0 for i, a in enumerate(NFW[k])]
    for b_lev in (1, 2):
        lam = eigen(l, b_lev); res = []
        for k, i, a in cands:
            am = [x % l for x in a]
            fits = [(d, sg, b) for d in DS for sg in (1, -1) for b in range(l - 1) if all(lam[Q] == sg * leg(d, Q) * pow(Q, b, l) * am[Q] % l for Q in lam)]
            if fits: res.append((k, i, len(fits), fits[0]))
            if fits and b_lev == 1: TOWER.setdefault(l, (k, i))
        say('   ℓ = %2d, level %d (%d primes): %s' % (l, b_lev, len(lam), '; '.join('weight %d, a₂ = %+d, a₃ = %+d: %d way%s, e.g. d = %d, sign %+d, Q^%d' % (
            k, NFW[k][i][2], NFW[k][i][3], nw, '' if nw == 1 else 's', *ft) for k, i, nw, ft in res) or 'no newform fits'))

# ---------- 3. which elliptic curves have the same numbers on lane ℓ ----------
say(); say('3  CURVES WITH THE SAME NUMBERS: a_Q(f) ≡ ±(d|Q)·Q^b·a_Q(E) (mod ℓ), every class of Cremona\'s tables (conductor < 70,000)')
LINE = re.compile(r'^(\d+)\s+([a-z]+)\s+1\s+\[(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+)\]')
labs, conds, AI = [], [], []
for fn in sorted(os.listdir(SRC)):
    if fn.startswith('allcurves.'):
        for line in open(os.path.join(SRC, fn), encoding='utf-8'):
            m = LINE.match(line)
            if m: labs.append(m.group(1) + m.group(2)); conds.append(int(m.group(1))); AI.append([int(m.group(i)) for i in range(3, 8)])
conds = np.array(conds, dtype=np.int64); AI = np.array(AI, dtype=object); NC = len(labs)
QF = [Q for Q in sp.primerange(5, 98)]
say('   %d isogeny classes; filter primes Q = 5 to 97; every survivor then checked at all primes up to %d' % (NC, LEN - 11))
t = time.time(); AQ = np.zeros((NC, len(QF)), dtype=np.int64)
for qi, Q in enumerate(QF):
    xs = np.arange(Q, dtype=np.int64); sq = np.zeros(Q, dtype=np.int64); sq[(xs * xs) % Q] = 1
    legt = 2 * sq - 1; legt[0] = 0
    r = np.array([[int(v) % Q for v in row] for row in AI], dtype=np.int64)
    a1, a2, a3, a4, a6 = r.T
    for c1 in (0, 1):
        for c3 in (0, 1):
            for c2 in range(Q):
                sel = np.nonzero((a1 % 2 == c1) & (a3 % 2 == c3) & (a2 == c2))[0]
                if len(sel) == 0: continue
                # D(x) = (a1 x + a3)² + 4(x³ + a2 x² + a4 x + a6), for the given residues
                base = ((c1 * xs + c3) ** 2 + 4 * (xs ** 3 + c2 * xs * xs)) % Q
                D = (base[None, :] + 4 * (a4[sel][:, None] * xs[None, :] + a6[sel][:, None])) % Q
                AQ[sel, qi] = -legt[D].sum(axis=1)
say('   point counts for all classes at the filter primes: %.0f s' % (time.time() - t))


def full_trace(curve, Q):
    a1, a2, a3, a4, a6 = [int(v) for v in curve]; xs = np.arange(Q, dtype=np.int64); x2 = xs * xs % Q; x3 = x2 * xs % Q
    D = ((a1 * xs + a3) ** 2 + 4 * (x3 + (a2 % Q) * x2 + (a4 % Q) * xs + (a6 % Q))) % Q; sq = np.zeros(Q, dtype=bool); sq[x2] = True
    return -int(np.where(D == 0, 0, np.where(sq[D], 1, -1)).sum())


for l in PRIMES:
    if l not in TOWER: say('   ℓ = %2d: (no newform matched the data)' % l); continue
    kw, iw = TOWER[l]; a = [x % l for x in NFW[kw][iw]]; qcols = [qi for qi, Q in enumerate(QF) if Q != l]
    bad = np.zeros((NC, len(QF)), dtype=bool)
    for qi, Q in enumerate(QF): bad[:, qi] = (conds % Q == 0)
    hits = []
    for d in DS:
        for sg in (1, -1):
            for b in range(l - 1):
                ok = np.ones(NC, dtype=bool)
                for qi in qcols:
                    Q = QF[qi]; coef = sg * leg(d, Q) * pow(Q, b, l) % l
                    ok &= bad[:, qi] | ((coef * AQ[:, qi] - a[Q]) % l == 0)
                    if not ok.any(): break
                for ci in np.nonzero(ok)[0]:
                    if (~bad[ci, qcols]).sum() < 15: continue
                    full = all((sg * leg(d, Q) * pow(Q, b, l) * full_trace(AI[ci], Q) - a[Q]) % l == 0
                               for Q in sp.primerange(5, LEN - 10) if Q != l and conds[ci] % Q)
                    if full: hits.append((labs[ci], int(conds[ci]), [int(v) for v in AI[ci]], d, sg, b))
    byclass = collections.OrderedDict()
    for lab, N, cur, d, sg, b in sorted(hits, key=lambda h: h[1]): byclass.setdefault((lab, N, tuple(cur)), []).append((d, sg, b))
    say('   ℓ = %2d (the tower\'s newform, weight %d): %s' % (l, kw, 'no curve' if not byclass else '%d classes, smallest conductors:' % len(byclass)))
    for (lab, N, cur), ways in list(byclass.items())[:4]:
        say('      %s = %s, conductor %d = %s: d = %d, sign %+d, Q^%d' % (lab, list(cur), N, sp.factorint(N), *ways[0]))
say(); say('seconds: %.0f' % (time.time() - t0))
