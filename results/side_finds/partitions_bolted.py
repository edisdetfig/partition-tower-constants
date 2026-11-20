# A dig, 11 Nov 2026: WHAT THE LEVEL-0 LINK CARRIES. For 5, p(25N − 1) ≡ 25·N(N − 3)·p(N) (mod 125) (partitions_tower_link.txt):
# the tower hangs from p(N) itself. For 7 to 31 nothing of that kind exists; their towers have a constant only from level 1
# down. Here, for every prime from 5 to 31: (1) how much knowing p(N) on lane ℓ tells about the value two levels down
# (24N′ − 1 = ℓ²·(24N − 1), the law divided out), as conditional shares and as mutual information in bits, against the same
# measure with the deeper values shuffled; (2) the step rule's multipliers at level 0 (Atkin's class, where p(N) has
# one, partitions_cells.txt) and two levels down, side by side: a link can only join levels that follow the same rule.
# p(n) for every n up to 2,400,000 by Euler's pentagonal rule. Shares, bits and multipliers COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_bolted.py [--keep]
import os, sys, time, math, random, collections
import numpy as np
import sympy as sp
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_bolted.txt'), 'x', encoding='utf-8') if KEEP else None


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
LANES = (5 ** 4 * 7 ** 4 * 11 ** 4 * 13, 17 * 19 * 23 * 29 * 31)
PP = []
for L in LANES:
    p = np.zeros(NMAX + 1, dtype=np.int64); p[0] = 1; c = 0
    for n in range(1, NMAX + 1):
        while c < len(G) and G[c] <= n: c += 1
        p[n] = int(np.dot(S[:c], p[n - G[:c]])) % L
    PP.append(p)
say('p(n) for every n up to %d on the lanes 5⁴·7⁴·11⁴·13 and 17·19·23·29·31 (%.0f s)' % (NMAX, time.time() - t0))
lane_of = lambda l: PP[0] if l <= 13 else PP[1]
LAW = lambda l, b: b if l in (5, 11) else (b // 2 + 1 if (l == 7 and b >= 1) else 0)
leg = lambda a, q: int(sp.jacobi_symbol(a % q, q))


def vl(x, l):
    v = 0
    while x % l == 0: x //= l; v += 1
    return v


def mutual(pairs):
    """mutual information (bits) between x and y within each group g, averaged over groups; pairs = (g, x, y)"""
    by = collections.defaultdict(list)
    for g, x, y in pairs: by[g].append((x, y))
    tot = 0.0; n_all = len(pairs)
    for g, xy in by.items():
        n = len(xy); cx = collections.Counter(x for x, _ in xy); cy = collections.Counter(y for _, y in xy); cxy = collections.Counter(xy)
        tot += sum(c_ / n_all * math.log2(c_ * n / (cx[x] * cy[y])) for (x, y), c_ in cxy.items())
    return tot


say(); say('1  HOW MUCH p(N) TELLS ABOUT THE VALUE TWO LEVELS DOWN (on lane ℓ; N with ℓ not dividing 24N − 1; grouped by N mod ℓ)')
say('   share of deeper values that are multiples of ℓ, when p(N) is / is not one;  information in bits (shuffled: the same')
say('   measure after shuffling the deeper values within each group, the noise floor);  the most there is to know: log₂ ℓ')
rng = random.Random(1)
for l in (5, 7, 11, 13, 17, 19, 23, 29, 31):
    p = lane_of(l); sh = (l * l - 1) // 24; pairs = []
    for N in range(1, (NMAX + sh) // (l * l) + 1):
        if (24 * N - 1) % l == 0: continue
        x = int(p[N]) % l; y = (int(p[l * l * N - sh]) // l ** LAW(l, 2)) % l; pairs.append((N % l, x, y))
    z = [y == 0 for g, x, y in pairs if x == 0]; nz = [y == 0 for g, x, y in pairs if x != 0]
    shuf = []; by = collections.defaultdict(list)
    for i, (g, x, y) in enumerate(pairs): by[g].append(i)
    ys = [y for _, _, y in pairs]
    for g, idx in by.items():
        vals = [ys[i] for i in idx]; rng.shuffle(vals)
        for i, v in zip(idx, vals): shuf.append((g, pairs[i][1], v))
    say('   ℓ = %2d: %6d pairs;  p(N) a multiple: %.4f;  not: %.4f;  information %.4f bits (shuffled %.4f; at most %.2f)' % (
        l, len(pairs), sum(z) / len(z), sum(nz) / len(nz), mutual(pairs), mutual(shuf), math.log2(l)))

say(); say('2  THE STEP RULE AT LEVEL 0 (Atkin\'s class) AND TWO LEVELS DOWN: the multipliers (μ₊, μ₋) on lane ℓ, one per class of')
say('   24N* − 1 at Q (levels: N* with ℓ not dividing 24N* − 1, in Atkin\'s class; and N* with ℓ² exactly dividing it, the law divided')
say('   out; for 5 the square-rule class is left out)')


def mults(l, b, Q):
    p = lane_of(l); g = collections.defaultdict(set); A = l ** b; r0 = pow(24, -1, A) if b else 1
    for Ns in range(r0, (NMAX + (Q * Q - 1) // 24) // (Q * Q) + 1, A if b else 1):
        M = 24 * Ns - 1
        if vl(M, l) != b or M % Q == 0: continue
        if b == 0 and leg(-M, l) != -1: continue                       # Atkin's class at level 0
        if b == 2 and l == 5 and leg(M // 25, 5) == 1: continue        # 5's square rule at level 2
        x = (int(p[Ns]) // l ** LAW(l, b)) % l
        if x == 0: continue
        N = Q * Q * Ns - (Q * Q - 1) // 24; y = (int(p[N]) // l ** LAW(l, b)) % l
        g[leg(M, Q)].add(y * pow(x, -1, l) % l)
    return tuple(next(iter(g[e])) if len(g[e]) == 1 else None for e in (1, -1))


for l in (5, 7, 13):
    rows = []; same = 0; tot = 0
    for Q in (q for q in sp.primerange(7, 60) if q != l):
        m0, m2 = mults(l, 0, Q), mults(l, 2, Q)
        if None in m0 or None in m2: continue
        tot += 1; same += m0 == m2; rows.append('%d: %s / %s' % (Q, m0, m2))
    say('   ℓ = %2d: the same at %d of %d primes;  Q: level 0 / level 2:  %s' % (l, same, tot, ';  '.join(rows[:8])))
say('   (for 11, 17, 19, 23, 29 and 31, p(N) has no step rule at level 0 at all: partitions_cells.txt, section 1)')
say(); say('seconds: %.0f' % (time.time() - t0))
