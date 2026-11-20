# A dig, 11 Nov 2026: JAMESON'S CURVE ON LANE 49. partitions_multipliers.txt found that for 7, from level 2 on, the step
# rule's eigenvalue on lane 49 is (3|Q)·Q¹⁹·a_Q, where a_Q = 1 + Q − #E7(F_Q) counts the points of Jameson's curve
# E7: y² + xy + y = x³ + x² − 4x + 5. Her theorem gives this on lane 7 only, with Q¹ in place of Q¹⁹. The exponent is the one
# that puts the weight at −2, the weight of p(n)'s own series under the Shimura lift: 2 + 2·1 ≡ −2 (mod 6) and
# 2 + 2·19 ≡ −2 (mod 42) (my reading). Here: the formula at levels 2, 3, 4 for primes far beyond the data; what it predicts
# (laws with TWO extra factors of 7 from a single prime Q, where a_Q ≡ 0 or ±1 on lane 49); and the share of such primes.
# Exact values of p(N) by the Hardy–Ramanujan–Rademacher formula (sympy). Checks COUNTED; predictions DERIVED from the
# formula, then counted; the share DERIVED by counting matrices of GL2(Z/49).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_curve49.py [--keep]
import os, sys, time
import numpy as np
import sympy as sp
from sympy.functions.combinatorial.numbers import partition
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_curve49.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time()
leg = lambda a, q: int(sp.jacobi_symbol(a % q, q))
LAW7 = lambda b: b // 2 + 1


def vl(x, l):
    v = 0
    while x % l == 0: x //= l; v += 1
    return v


def trace(Q, curve=(1, 1, 1, -4, 5)):
    a1, a2, a3, a4, a6 = curve; xs = np.arange(Q, dtype=np.int64); x2 = xs * xs % Q; x3 = x2 * xs % Q
    D = ((a1 * xs + a3) ** 2 + 4 * (x3 + a2 * x2 + a4 * xs + a6)) % Q; sq = np.zeros(Q, dtype=bool); sq[x2] = True
    return -int(np.where(D == 0, 0, np.where(sq[D], 1, -1)).sum())


P = {}
def pv(N):
    if N not in P: P[N] = int(partition(N))
    return P[N]


def c49(N, b):
    """p(N) divided by Watson's power at level b, on lane 49"""
    x = pv(N); assert x % 7 ** LAW7(b) == 0
    return (x // 7 ** LAW7(b)) % 49


say('1  THE FORMULA ON LANE 49: for N* at exact level b (24N* − 1 = 7ᵇ·m, 7 not dividing m; at odd levels outside the square rule)')
say('   and N with 24N − 1 = Q²·(24N* − 1): c(N) ≡ μ_ε·c(N*), with μ_ε = λ + ε·δ, ε = the sign of 24N* − 1 at Q, and the prediction')
say('   λ = (3|Q)·Q¹⁹·a_Q, δ = −(−3|Q)·Q⁴⁰ (on lane 49). One member of each class per prime, exact values throughout.')
for b, qs in ((2, list(sp.primerange(11, 400))), (3, list(sp.primerange(11, 200))), (4, list(sp.primerange(11, 100)))):
    A = 7 ** b; r0 = pow(24, -1, A); ok = tot = 0; bad = []
    for Q in qs:
        lam = leg(3, Q) * pow(Q, 19, 49) * trace(Q) % 49; dl = -leg(-3, Q) * pow(Q, 40, 49) % 49; good = True
        for eps in (1, -1):
            Ns = r0
            while True:
                M = 24 * Ns - 1
                if vl(M, 7) == b and M % Q and leg(M, Q) == eps and not (b % 2 and leg(M // A, 7) == 1) and c49(Ns, b) % 7: break
                Ns += A
            N = Q * Q * Ns - (Q * Q - 1) // 24
            good = good and c49(N, b) == (lam + eps * dl) * c49(Ns, b) % 49
        tot += 1; ok += good
        if not good: bad.append(Q)
    say('   level %d: %d of %d primes from 11 to %d right%s' % (b, ok, tot, qs[-1], '' if not bad else ' (wrong at %s)' % bad[:8]))

say(); say('2  WHAT IT PREDICTS: one prime, two extra factors of 7 at level 2 (Watson\'s law gives 7², Jameson\'s corollaries 7³)')
say('   (a) a_Q = 0 makes λ ≡ 0 on lane 49 (any power, in fact), so p(N) ≡ 0 (mod 7⁴) whenever 24N − 1 = 49·Q³·m with Q and 7 not')
say('       dividing m. The first 12 such N for each Q; Q = 41 (a_Q = −6) as the control:')
for Q in (31, 47, 41):
    vals = []; m = 1
    while len(vals) < 12:
        M = 49 * Q ** 3 * m
        if M % 24 == 23 and m % Q and m % 7: vals.append(vl(pv((M + 1) // 24), 7))
        m += 1
    say('       Q = %d (a_Q = %+d): powers of 7 in p(N): %s' % (Q, trace(Q), vals))
say('   (b) a_Q ≡ ±1 on lane 49 (a_Q = ±48 or ±50, possible only from Q = 577 on) makes μ_ε ≡ 0 for one class ε: 8 members of that')
say('       class at level 2 for each such Q, and 8 of the other class:')
cand = [(Q, a) for Q, a in ((Q, trace(Q)) for Q in sp.primerange(577, 2000)) if a % 49 in (1, 48)]
say('       primes 577 to 2000 with a_Q ≡ ±1 on lane 49: %s' % ', '.join('%d (%+d)' % qa for qa in cand))
for Q, a in cand[:6]:
    lam = leg(3, Q) * pow(Q, 19, 49) * a % 49; dl = -leg(-3, Q) * pow(Q, 40, 49) % 49
    eps = [e for e in (1, -1) if (lam + e * dl) % 49 == 0][0]; vals, ctl = [], []; Ns = 47
    while len(vals) < 8 or len(ctl) < 8:
        M = 24 * Ns - 1
        if vl(M, 7) == 2 and M % Q:
            N = Q * Q * Ns - (Q * Q - 1) // 24
            if leg(M, Q) == eps and len(vals) < 8: vals.append(vl(pv(N), 7))
            elif leg(M, Q) == -eps and len(ctl) < 8: ctl.append(vl(pv(N), 7))
        Ns += 49
    say('       Q = %4d (a_Q = %+d), class %+d: %s;  other class: %s' % (Q, a, eps, vals, ctl))

say(); say('3  HOW MANY PRIMES GIVE SUCH LAWS. Counted: primes 5 ≤ Q ≤ 30,000 (Q ≠ 7) with a_Q ≡ 0, +1, −1 on lane 49. Derived: the share')
say('   of matrices in GL2(Z/49) with those traces (the share of primes, if the curve\'s 7-adic Galois image is all of GL2(Z₇))')
qs = [Q for Q in sp.primerange(5, 30001) if Q != 7]; tr = np.array([trace(Q) % 49 for Q in qs])
x = np.arange(49); a_, d_ = np.meshgrid(x, x, indexing='ij'); tr_ad = (a_ + d_) % 49; ad = a_ * d_ % 49
bc = (x[:, None] * x[None, :]) % 49; cnt_bc = np.bincount(bc.ravel(), minlength=49)
total = 0; by_tr = np.zeros(49, dtype=np.int64)
for t_ in range(49):
    for det_mod in range(49):
        if det_mod % 7 == 0: continue
        need = (ad - det_mod) % 49; n_ = (cnt_bc[need] * (tr_ad == t_)).sum(); by_tr[t_] += n_
total = by_tr.sum()
say('   |GL2(Z/49)| = %d (counted);  traces 0, +1, −1: %d, %d, %d each' % (total, by_tr[0], by_tr[1], by_tr[48]))
for t_, nm in ((0, '0'), (1, '+1'), (48, '−1')):
    say('   a_Q ≡ %s on lane 49: counted %.4f of %d primes;  matrices: %.4f' % (nm, (tr == t_).mean(), len(qs), by_tr[t_] / total))
say('   a_Q is even for %d of the %d primes (counted), so a_Q ≡ ±1 on lane 49 needs |a_Q| ≥ 48, and with |a_Q| ≤ 2√Q those two shares' % (
    int(sum(1 for Q in qs if trace(Q) % 2 == 0)), len(qs)))
say('   come only from primes above 577')
say(); say('seconds: %.0f' % (time.time() - t0))
