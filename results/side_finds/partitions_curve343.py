# A dig, 11 Nov 2026: HOW FAR JAMESON'S CURVE REACHES. partitions_curve49.txt: for 7, from level 2 on, the step rule's
# eigenvalue on lane 49 is (3|Q)·Q¹⁹·(1 + Q − #E7(F_Q)). Here lane 343 at levels 3 to 6: is there still one multiplier per
# class (3 members per class, primes 11 to 59, at level 6 to 29), and does any formula ±(3|Q)·Q^a·(1 + Q − #E7(F_Q)),
# a = 0 to 293, fit the eigenvalue? Exact values of p(N) (Hardy–Ramanujan–Rademacher formula, sympy). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_curve343.py [--keep]
import os, sys, time
import numpy as np
import sympy as sp
from sympy.functions.combinatorial.numbers import partition
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_curve343.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); leg = lambda a, q: int(sp.jacobi_symbol(a % q, q)); LAW7 = lambda b: b // 2 + 1; MOD = 343


def vl(x, l):
    v = 0
    while x % l == 0: x //= l; v += 1
    return v


def trace(Q, curve=(1, 1, 1, -4, 5)):
    a1, a2, a3, a4, a6 = curve; xs = np.arange(Q, dtype=np.int64); x2 = xs * xs % Q; x3 = x2 * xs % Q
    D = ((a1 * xs + a3) ** 2 + 4 * (x3 + a2 * x2 + a4 * xs + a6)) % Q; sq = np.zeros(Q, dtype=bool); sq[x2] = True
    return -int(np.where(D == 0, 0, np.where(sq[D], 1, -1)).sum())


P = {}
def c(N, b):
    if N not in P: P[N] = int(partition(N))
    x = P[N]; assert x % 7 ** LAW7(b) == 0
    return (x // 7 ** LAW7(b)) % MOD


say('level b: N* at exact level b (odd levels outside the square rule), N with 24N − 1 = Q²·(24N* − 1); c = p/7^law on lane 343')
for b in (3, 4, 5, 6):
    A = 7 ** b; r0 = pow(24, -1, A); MU = {}; t = time.time()
    for Q in sp.primerange(11, 60 if b < 6 else 30):
        g = {}
        for eps in (1, -1):
            vals = []; Ns = r0
            while len(vals) < 3:
                M = 24 * Ns - 1
                if vl(M, 7) == b and M % Q and leg(M, Q) == eps and not (b % 2 and leg(M // A, 7) == 1):
                    x = c(Ns, b)
                    if x % 7: vals.append(c(Q * Q * Ns - (Q * Q - 1) // 24, b) * pow(x, -1, MOD) % MOD)
                Ns += A
            g[eps] = vals
        MU[Q] = g
    ok = sum(1 for g in MU.values() if all(len(set(v)) == 1 for v in g.values()))
    fits = [(sg, a) for a in range(294) for sg in (1, -1) if ok == len(MU) and all(
        (MU[Q][1][0] + MU[Q][-1][0]) * pow(2, -1, MOD) % MOD == sg * leg(3, Q) * pow(Q, a, MOD) * trace(Q) % MOD for Q in MU)]
    say('   level %d: one multiplier per class for %d of %d primes;  E7 formula fits on lane 343: %s  (%.0f s)' % (
        b, ok, len(MU), fits if fits else 'none', time.time() - t))
say('   (on lane 49 the same search finds a = 19 at levels 2, 3, 4: partitions_curve49.txt)')
say(); say('seconds: %.0f' % (time.time() - t0))
