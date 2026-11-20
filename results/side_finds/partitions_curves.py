# A dig, 11 Nov 2026: A CURVE FOR EVERY PRIME? Jameson proved that inside the laws of 7 and 11 the step rule's eigenvalue
# λ(Q) on lane ℓ is (3|Q)·Q^a·(1 + Q − #E(F_Q)) for one elliptic curve (E7 = 42a1, E11 = 726d1 in Cremona's labels; a = 1, 3);
# for 5, Lovejoy and Ono have 1 + Q in place of the count. Every prime from 5 to 31 follows a step rule at levels 1 and 2
# (section 1). Which rule? Here every elliptic curve over Q whose conductor is below 70,000 and has no prime factor other
# than 2, 3 and ℓ is tested against the eigenvalues, with any twist χ = (d|Q), d ∈ {±1, ±2, ±3, ±6}, a sign and a power Q^a:
#     λ(Q) ≡ ±(d|Q)·Q^a·a_Q(E)  (mod ℓ),   a_Q(E) = 1 + Q − #E(F_Q), counted,
# at every prime Q from 5 to 113 where the data give the eigenvalue (about 26 primes, so a fit by chance has odds near ℓ^−26).
# Curves: J. Cremona's tables (github.com/JohnCremona/ecdata, files allcurves.00000-09999 to 60000-69999, downloaded to the
# scratchpad; the folder is the argument), one curve per isogeny class. Eigenvalues and point counts COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_curves.py <folder> [--keep]
import os, re, sys, time, collections
import numpy as np
import sympy as sp
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv; SRC = [a for a in sys.argv[1:] if a != '--keep'][0]
fh = open(os.path.join(HERE, 'partitions_curves.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); NMAX = 2_400_000
G, S, k = [], [], 1
while True:
    a, b = k * (3 * k - 1) // 2, k * (3 * k + 1) // 2
    if a > NMAX: break
    s = 1 if k % 2 else -1; G.append(a); S.append(s)
    if b <= NMAX: G.append(b); S.append(s)
    k += 1
o = np.argsort(G); G = np.array(G, dtype=np.int64)[o]; S = np.array(S, dtype=np.int64)[o]
LANES = (5 ** 3 * 7 ** 3 * 11 ** 3 * 13, 17 * 19 * 23 * 29 * 31); PP = []
for L in LANES:
    p = np.zeros(NMAX + 1, dtype=np.int64); p[0] = 1; c = 0
    for n in range(1, NMAX + 1):
        while c < len(G) and G[c] <= n: c += 1
        p[n] = int(np.dot(S[:c], p[n - G[:c]])) % L
    PP.append(p)
say('p(n) for every n up to %d on the lanes 5³·7³·11³·13 and 17·19·23·29·31 (%.0f s)' % (NMAX, time.time() - t0))
PRIMES = (5, 7, 11, 13, 17, 19, 23, 29, 31)
LAW = lambda l, b: b if l in (5, 11) else (b // 2 + 1 if l == 7 else 0)
leg = lambda a, q: int(sp.jacobi_symbol(a % q, q))


def vl(x, l):
    v = 0
    while x % l == 0: x //= l; v += 1
    return v


def eigen(l, b, qmax=113):
    p = PP[0] if l <= 13 else PP[1]; A = l ** b; r0 = pow(24, -1, A); lam = {}; bad = []
    for Q in sp.primerange(5, qmax + 1):
        if Q == l: continue
        g = collections.defaultdict(set); cnt = collections.Counter()
        for Ns in range(r0, (NMAX + (Q * Q - 1) // 24) // (Q * Q) + 1, A):
            M = 24 * Ns - 1
            if vl(M, l) != b or M % Q == 0: continue
            if (l == 5 and b % 2 == 0 or l == 7 and b % 2 == 1) and leg(M // A, l) == 1: continue      # square-rule classes
            x = (int(p[Ns]) // l ** LAW(l, b)) % l
            if x == 0: continue
            y = (int(p[Q * Q * Ns - (Q * Q - 1) // 24]) // l ** LAW(l, b)) % l; e = leg(M, Q)
            g[e].add(y * pow(x, -1, l) % l); cnt[e] += 1
        if len(g) < 2 or min(cnt.values()) < 2: continue
        if any(len(v) > 1 for v in g.values()): bad.append(Q); continue
        mp, mm = next(iter(g[1])), next(iter(g[-1])); i2 = pow(2, -1, l)
        if (mp - mm) * i2 % l != (-leg(-3, Q) * pow(Q, -2, l)) % l: bad.append(Q); continue
        lam[Q] = (mp + mm) * i2 % l
    return lam, bad


say(); say('1  THE STEP RULE OF EVERY PRIME: one multiplier per class, with Atkin\'s middle term −(−3|Q)·Q⁻², at levels 1 and 2')
LAM = {}
for l in PRIMES:
    for b in (1, 2):
        lam, bad = eigen(l, b); LAM[(l, b)] = lam
        say('   ℓ = %2d, level %d: %d primes Q fit (Q from 5 to %d), %d do not%s' % (l, b, len(lam), max(lam), len(bad), '' if not bad else ' (%s)' % bad))

say(); say('2  IS IT THE SIMPLE 1 + Q RULE?  λ(Q) = ±(d|Q)·Q^a·(1 + Q^c) on lane ℓ, any d, a, c:')
DS = (1, -1, 2, -2, 3, -3, 6, -6)
for l in PRIMES:
    res = []
    for b in (1, 2):
        lam = LAM[(l, b)]; fits = [(d, sg, a_, c_) for d in DS for sg in (1, -1) for a_ in range(l - 1) for c_ in range(l - 1)
                                  if all(lam[Q] == sg * leg(d, Q) * pow(Q, a_, l) * (1 + pow(Q, c_, l)) % l for Q in lam)]
        res.append('level %d: %s' % (b, 'yes, e.g. d = %d, sign %+d, Q^%d·(1 + Q^%d)' % fits[0] if fits else 'no'))
    say('   ℓ = %2d: %s' % (l, ';  '.join(res)))

say(); say('3  THE CURVES: λ(Q) ≡ ±(d|Q)·Q^a·a_Q(E) (mod ℓ) at every Q, for which curves E (conductor < 70,000, only 2, 3, ℓ)?')
LINE = re.compile(r'^(\d+)\s+([a-z]+)\s+(\d+)\s+\[(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+)\]')
classes = {}
for fn in sorted(os.listdir(SRC)):
    if not fn.startswith('allcurves.'): continue
    for line in open(os.path.join(SRC, fn), encoding='utf-8'):
        m = LINE.match(line)
        if m and m.group(3) == '1': classes['%s%s' % (m.group(1), m.group(2))] = (int(m.group(1)), tuple(int(m.group(i)) for i in range(4, 9)))
say('   isogeny classes read: %d (conductors up to %d)' % (len(classes), max(N for N, _ in classes.values())))


def trace(curve, Q):
    a1, a2, a3, a4, a6 = curve; xs = np.arange(Q, dtype=np.int64); x2 = xs * xs % Q; x3 = x2 * xs % Q
    D = ((a1 * xs + a3) ** 2 + 4 * (x3 + a2 * x2 + a4 * xs + a6)) % Q; sq = np.zeros(Q, dtype=bool); sq[x2] = True
    return -int(np.where(D == 0, 0, np.where(sq[D], 1, -1)).sum())


FOUND = {}
for l in PRIMES:
    mine = {lab: v for lab, v in classes.items() if set(sp.primefactors(v[0])) <= {2, 3, l}}
    for b in (1, 2):
        lam = LAM[(l, b)]; Qs = sorted(lam); hits = []
        for lab, (N, cur) in sorted(mine.items(), key=lambda kv: kv[1][0]):
            aq = {Q: trace(cur, Q) for Q in Qs}
            for d in DS:
                for sg in (1, -1):
                    for a_ in range(l - 1):
                        if all(lam[Q] == sg * leg(d, Q) * pow(Q, a_, l) * aq[Q] % l for Q in Qs): hits.append((lab, N, cur, d, sg, a_))
        FOUND[(l, b)] = hits
        rule = [a_ for a_ in range(l - 1) if (2 + 2 * a_) % (l - 1) == (l - 1 - 2) % (l - 1)]
        say('   ℓ = %2d, level %d (%d classes tried, %d primes): %s' % (l, b, len(mine), len(Qs), 'no curve fits' if not hits else ''))
        shown = collections.OrderedDict()
        for lab, N, cur, d, sg, a_ in hits: shown.setdefault((lab, N, cur), []).append((d, sg, a_))
        for (lab, N, cur), ways in list(shown.items())[:8]:
            say('      %s = %s (conductor %d = %s): %s;  a ≡ −2 (mod (ℓ − 1)/2) for a in %s' % (lab, list(cur), N, sp.factorint(N),
                ', '.join('d = %d, sign %+d, Q^%d' % w for w in ways[:4]), rule))
        if len(shown) > 8: say('      … and %d more classes' % (len(shown) - 8))
say(); say('seconds: %.0f' % (time.time() - t0))
