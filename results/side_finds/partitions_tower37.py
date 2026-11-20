# A dig, 11 Nov 2026: IS THE TOWER'S RULE STILL ONE RULE AT 37, 41, 43, 47? Up to 31 the step rule at level 1 had one
# multiplier per class (partitions_curves.txt), i.e. the tower was a single Hecke eigenform mod ℓ. Ahlgren, Allen and Tang note
# that from 37 on the relevant space of forms is no longer one-dimensional (two-dimensional for 37). Here, from the partition
# numbers: at level 1 (ℓ exactly dividing 24N* − 1) and every prime Q from 5 to 43 with at least 4 pairs in each class,
# is p(N) ≡ μ·p(N*) (mod ℓ) with one μ per class, where 24N − 1 = Q²·(24N* − 1)? p(n) for every n up to 2,400,000 by Euler's
# pentagonal rule (lane 37·41·43·47). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_tower37.py [--keep]
import os, sys, time, collections
import numpy as np
import sympy as sp
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_tower37.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); NMAX = 2_400_000; L = 37 * 41 * 43 * 47
G, S, k = [], [], 1
while True:
    a, b = k * (3 * k - 1) // 2, k * (3 * k + 1) // 2
    if a > NMAX: break
    s = 1 if k % 2 else -1; G.append(a); S.append(s)
    if b <= NMAX: G.append(b); S.append(s)
    k += 1
o = np.argsort(G); G = np.array(G, dtype=np.int64)[o]; S = np.array(S, dtype=np.int64)[o]
p = np.zeros(NMAX + 1, dtype=np.int64); p[0] = 1; c = 0
for n in range(1, NMAX + 1):
    while c < len(G) and G[c] <= n: c += 1
    p[n] = int(np.dot(S[:c], p[n - G[:c]])) % L
say('p(n) for every n up to %d on the lane 37·41·43·47 (%.0f s)' % (NMAX, time.time() - t0))
leg = lambda a, q: int(sp.jacobi_symbol(a % q, q))


def vl(x, l):
    v = 0
    while x % l == 0: x //= l; v += 1
    return v


say(); say('LEVEL 1: one multiplier per class (ε = sign of 24N* − 1 at Q)?  For 13 to 31 every prime Q had one (partitions_curves.txt)')
for l in (37, 41, 43, 47):
    r0 = pow(24, -1, l); single = tot = 0; rows = []
    for Q in sp.primerange(5, 44):
        if Q == l: continue
        g = collections.defaultdict(list)
        for Ns in range(r0, (NMAX + (Q * Q - 1) // 24) // (Q * Q) + 1, l):
            M = 24 * Ns - 1
            if vl(M, l) != 1 or M % Q == 0 or int(p[Ns]) % l == 0: continue
            g[leg(M, Q)].append(int(p[Q * Q * Ns - (Q * Q - 1) // 24]) % l * pow(int(p[Ns]) % l, -1, l) % l)
        if len(g) < 2 or min(len(v) for v in g.values()) < 4: continue
        tot += 1; one = all(len(set(v)) == 1 for v in g.values()); single += one
        rows.append('%d: %s' % (Q, 'one each' if one else '/'.join(str(len(set(g[e]))) for e in (1, -1)) + ' values'))
    say('   ℓ = %d: one multiplier per class at %d of %d primes;  %s' % (l, single, tot, ';  '.join(rows[:8])))
say(); say('seconds: %.0f' % (time.time() - t0))
