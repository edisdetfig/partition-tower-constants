# A dig, 11 Nov 2026: THE CURVE BEHIND 13. partitions_curves.txt found that the step rule of 13 (levels 1 and 2; there are no
# laws to divide out for 13) has the eigenvalue
#     λ(Q) ≡ (3|Q)·Q¹⁰·a_Q (mod 13),   a_Q = 1 + Q − #E13(F_Q),   E13 = 1014f1: y² + xy = x³ − 62x + 324 (conductor 2·3·13²),
# the same shape as Jameson's curves for 7 (42a1, Q¹) and 11 (726d1, Q³), with the power again the one that puts the weight at
# −2 (2 + 2·10 ≡ −2 mod 12). Here: (1) the formula at primes far beyond the data, with exact values of p(N) (Hardy–Ramanujan–
# Rademacher formula, sympy), at levels 1, 2 and 3; (2) what it predicts: laws p(N) ≡ 0 (mod 13) where none is forced,
# from every prime with a_Q ≡ 0 or ±1 on lane 13, checked with exact values; (3) the share of such primes.
# Checks COUNTED; predictions DERIVED from the formula, then counted; the share DERIVED by counting matrices of GL2(F13).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_curve13.py [--keep]
import os, sys, time
import numpy as np
import sympy as sp
from sympy.functions.combinatorial.numbers import partition
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_curve13.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); L = 13; E13 = (1, 0, 0, -62, 324)
leg = lambda a, q: int(sp.jacobi_symbol(a % q, q))


def vl(x, l):
    v = 0
    while x % l == 0: x //= l; v += 1
    return v


def trace(Q, curve=E13):
    a1, a2, a3, a4, a6 = curve; xs = np.arange(Q, dtype=np.int64); x2 = xs * xs % Q; x3 = x2 * xs % Q
    D = ((a1 * xs + a3) ** 2 + 4 * (x3 + a2 * x2 + a4 * xs + a6)) % Q; sq = np.zeros(Q, dtype=bool); sq[x2] = True
    return -int(np.where(D == 0, 0, np.where(sq[D], 1, -1)).sum())


def disc(a1, a2, a3, a4, a6):
    b2, b4, b6 = a1 * a1 + 4 * a2, 2 * a4 + a1 * a3, a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6


P = {}
def pm(N):
    if N not in P: P[N] = int(partition(N)) % L
    return P[N]


lam = lambda Q: leg(3, Q) * pow(Q, 10, L) * trace(Q) % L
dlt = lambda Q: -leg(-3, Q) * pow(Q, -2, L) % L
d = disc(*E13)
say('E13 = 1014f1 = %s: discriminant %d = %s, bad primes %s;  a_Q for Q = 5 to 47: %s' % (
    list(E13), d, sp.factorint(d), sorted(sp.factorint(abs(d))), ', '.join('%d: %+d' % (Q, trace(Q)) for Q in sp.primerange(5, 48) if Q != 13)))

say(); say('1  THE FORMULA FAR BEYOND THE DATA: p(N) ≡ (λ + ε·δ)·p(N*) (mod 13), λ = (3|Q)·Q¹⁰·a_Q, δ = −(−3|Q)·Q⁻², for N* at exact')
say('   level b, 24N − 1 = Q²·(24N* − 1), ε the sign of 24N* − 1 at Q; two members of each class per prime, exact values')
for b, qmax in ((1, 600), (2, 300), (3, 150)):
    A = L ** b; r0 = pow(24, -1, A); ok = tot = 0; bad = []
    for Q in sp.primerange(5, qmax):
        if Q == L: continue
        good = True
        for eps in (1, -1):
            got = 0; Ns = r0
            while got < 2:
                M = 24 * Ns - 1
                if vl(M, L) == b and M % Q and leg(M, Q) == eps and pm(Ns):
                    N = Q * Q * Ns - (Q * Q - 1) // 24
                    good = good and pm(N) == (lam(Q) + eps * dlt(Q)) * pm(Ns) % L; got += 1
                Ns += A
        tot += 1; ok += good
        if not good: bad.append(Q)
    say('   level %d: %d of %d primes from 5 to %d right%s' % (b, ok, tot, qmax, '' if not bad else ' (wrong at %s)' % bad[:8]))

say(); say('2  WHAT IT PREDICTS: laws mod 13 at level 1 (where nothing is forced). a_Q ≡ ±1 makes one class\'s multiplier λ + ε·δ zero:')
say('   p(13·Q²·n + β) ≡ 0 for n in that class (progressions modulo 13·Q³). a_Q ≡ 0 makes λ zero: p(N) ≡ 0 whenever 24N − 1 = 13·Q³·m')
say('   with Q and 13 not dividing m (progressions modulo 13·Q⁴). Ten members of each (n from 1,000 on), and ten of the other class:')
lists = {i: [Q for Q in sp.primerange(5, 400) if Q != 13 and trace(Q) % L == i % L] for i in (0, 1, -1)}
say('   primes below 400 with a_Q ≡ 0: %s;  ≡ +1: %s;  ≡ −1: %s' % tuple(', '.join(map(str, lists[i])) for i in (0, 1, -1)))
beta = lambda Q: (11 * 13 * Q * Q + 1) // 24
for Q in [q for q in lists[1] + lists[-1] if q > 5][:4]:
    eps = [e for e in (1, -1) if (lam(Q) + e * dlt(Q)) % L == 0][0]; law, ctl = [], []; n = 1000
    while len(law) < 10 or len(ctl) < 10:
        M = 13 * (24 * n + 11)
        if (24 * n + 11) % 13 and M % Q:
            N = 13 * Q * Q * n + beta(Q); (law if leg(M, Q) == eps else ctl).append(N) if len(law if leg(M, Q) == eps else ctl) < 10 else None
        n += 1
    hit = sum(1 for N in law if pm(N) == 0); hc = sum(1 for N in ctl if pm(N) == 0)
    say('   Q = %3d (a_Q = %+d): class %+d: %d of 10 multiples of 13 (largest N = %d);  other class: %d of 10' % (Q, trace(Q), eps, hit, max(law), hc))
for Q in lists[0][:3]:
    law, m = [], 1
    while len(law) < 10:
        M = 13 * Q ** 3 * m
        if M % 24 == 23 and m % Q and m % 13: law.append((M + 1) // 24)
        m += 1
    say('   Q = %3d (a_Q = %+d): 24N − 1 = 13·Q³·m: %d of 10 multiples of 13 (largest N = %d)' % (Q, trace(Q), sum(1 for N in law if pm(N) == 0), max(law)))

say('   The first prime with a_Q = 0 is 59 (#E13(F₅₉) = 60). QUOTED: Ono\'s example (Annals of Math. 151, 2000, eq. (2), credited to')
say('   Atkin and O\'Brien), p(59⁴·13n + 111247) ≡ 0 (mod 13). DERIVED: 24·111247 − 1 = 13·59³, so it is the case Y ≡ 1 (mod 59)')
say('   of the curve\'s law 24N − 1 = 13·59³·Y. Every Y prime to 59 (13 may divide it: levels 2 and 3 follow the same formula) gives')
say('   a progression p(13·59⁴·n + B) ≡ 0 (mod 13): 58 values of B. Two members (n = 0, 1) of each, exact values:')
Bs = []; Y = 1
while len(Bs) < 58:
    M = 13 * 59 ** 3 * Y
    if M % 24 == 23 and Y % 59: Bs.append((M + 1) // 24)
    Y += 1
Bs = sorted(set(B % (13 * 59 ** 4) for B in Bs))
hits = sum(1 for B in Bs for n in (0, 1) if pm(13 * 59 ** 4 * n + B) == 0)
say('   %d values of B (111247 among them: %s): %d of %d members are multiples of 13 (largest N = %d)' % (
    len(Bs), 111247 in Bs, hits, 2 * len(Bs), 13 * 59 ** 4 + max(Bs)))

say(); say('3  HOW MANY PRIMES GIVE SUCH LAWS: counted over the primes 5 ≤ Q ≤ 30,000 (Q ≠ 13); derived from the matrices of GL2(F13)')
qs = [Q for Q in sp.primerange(5, 30001) if Q != 13]; tr = np.array([trace(Q) % L for Q in qs])
g = [(a_, b_, c_, d_) for a_ in range(L) for b_ in range(L) for c_ in range(L) for d_ in range(L) if (a_ * d_ - b_ * c_) % L]
by = {i: sum(1 for m4 in g if (m4[0] + m4[3]) % L == i % L) for i in (0, 1, -1)}
for i, nm in ((0, '0'), (1, '+1'), (-1, '−1')):
    say('   a_Q ≡ %s: counted %.4f of %d primes;  matrices %d of %d = %.4f' % (nm, (tr == i % L).mean(), len(qs), by[i], len(g), by[i] / len(g)))
say('   together: %.1f%% of all primes (matrices: %.1f%%)' % (100 * np.isin(tr, [0, 1, L - 1]).mean(), 100 * (by[0] + by[1] + by[-1]) / len(g)))
say(); say('seconds: %.0f' % (time.time() - t0))
