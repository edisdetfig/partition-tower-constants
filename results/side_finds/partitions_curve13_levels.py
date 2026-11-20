# A dig, 11 Nov 2026: THE 13 THEOREM AT EVERY LEVEL. partitions_curve13_proof.txt proved, for every prime Q ≥ 5 (Q ≠ 13) and
# EVERY n ≥ 0:  p(13Q²n + β) ≡ (3|Q)·Q⁴·[a_Q − ((−24n − 11)|Q)]·p(13n + 6) − Q⁹·p(13n/Q² + β′)  (mod 13), E = 1014b1.
# 13 has no Ramanujan law, so nothing is divided out on the way down: a partition number at level b (13ᵇ divides 24N − 1)
# is also a level-1 number p(13n + 6), with 13^(b−1) dividing 24n + 11. The theorem holds for every n, so it already covers
# every level (DERIVED); written in level-b terms the sign changes by (13|Q)^(b−1), because 24n + 11 = 13^(b−1)·m. No new Sturm
# check is needed. Shown here directly: the statement at n of exact level 2 and 3 (every case up to 2,400,000) and at levels
# 4 to 7 (exact values of p(N), Hardy–Ramanujan–Rademacher formula, sympy). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_curve13_levels.py [--keep]
import os, sys, time
import numpy as np
import sympy as sp
from sympy.functions.combinatorial.numbers import partition
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_curve13_levels.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); E = (1, 0, 1, -10482, 722308); NMAX = 2_400_000
leg = lambda a, q: int(sp.jacobi_symbol(a % q, q))


def vl(x, l):
    v = 0
    while x % l == 0: x //= l; v += 1
    return v


def trace(curve, Q):
    a1, a2, a3, a4, a6 = curve; xs = np.arange(Q, dtype=np.int64); x2 = xs * xs % Q; x3 = x2 * xs % Q
    D = ((a1 * xs + a3) ** 2 + 4 * (x3 + (a2 % Q) * x2 + (a4 % Q) * xs + (a6 % Q))) % Q; sq = np.zeros(Q, dtype=bool); sq[x2] = True
    return -int(np.where(D == 0, 0, np.where(sq[D], 1, -1)).sum())


G, S, kk = [], [], 1
while True:
    x1, x2 = kk * (3 * kk - 1) // 2, kk * (3 * kk + 1) // 2
    if x1 > NMAX: break
    s = 1 if kk % 2 else -1; G.append(x1); S.append(s)
    if x2 <= NMAX: G.append(x2); S.append(s)
    kk += 1
o = np.argsort(G); G = np.array(G, dtype=np.int64)[o]; S = np.array(S, dtype=np.int64)[o]
p = np.zeros(NMAX + 1, dtype=np.int64); p[0] = 1; c = 0
for n in range(1, NMAX + 1):
    while c < len(G) and G[c] <= n: c += 1
    p[n] = int(np.dot(S[:c], p[n - G[:c]])) % 13
say('p(n) mod 13 for every n up to %d (%.0f s)' % (NMAX, time.time() - t0))


def rhs(Q, n, pval):
    m = 24 * n + 11; r = leg(3, Q) * pow(Q, 4, 13) * (trace(E, Q) - leg(-m, Q)) * pval
    if m % (Q * Q) == 0: r -= pow(Q, 9, 13) * pval_of((13 * (m // (Q * Q)) + 1) // 24)
    return r % 13


def pval_of(N): return int(p[N]) if N <= NMAX else int(partition(N)) % 13


say(); say('LEVELS 2 AND 3 (n with 13 or 13² exactly dividing 24n + 11; every case with 13Q²n + β ≤ 2,400,000)')
for b in (2, 3):
    checks = wrong = 0; qs = set()
    for Q in sp.primerange(5, 400):
        if Q == 13: continue
        beta = (143 * Q * Q + 1) // 24
        for n in range(0, (NMAX - beta) // (13 * Q * Q) + 1):
            if vl(24 * n + 11, 13) != b - 1: continue
            checks += 1; qs.add(Q); wrong += (int(p[13 * Q * Q * n + beta]) - rhs(Q, n, int(p[13 * n + 6]))) % 13 != 0
    say('   level %d: %d cases at %d primes Q, wrong: %d' % (b, checks, len(qs), wrong))

say(); say('LEVELS 4 TO 7 (exact values: the first two n of exact level b, primes Q = 5, 7, 11, 17, 19, 23)')
for b in (4, 5, 6, 7):
    A = 13 ** (b - 1); n_list = []; n = 0
    r0 = (-11 * pow(24, -1, A)) % A                                   # 24n + 11 ≡ 0 (mod 13^(b−1))
    n = r0
    while len(n_list) < 2:
        if vl(24 * n + 11, 13) == b - 1: n_list.append(n)
        n += A
    checks = wrong = 0; t = time.time(); big = 0
    for Q in (5, 7, 11, 17, 19, 23):
        beta = (143 * Q * Q + 1) // 24
        for n in n_list:
            N = 13 * Q * Q * n + beta; big = max(big, N); checks += 1
            wrong += (int(partition(N)) % 13 - rhs(Q, n, int(partition(13 * n + 6)) % 13)) % 13 != 0
    say('   level %d: %d cases, wrong: %d  (largest N = %d; %.0f s)' % (b, checks, wrong, big, time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
