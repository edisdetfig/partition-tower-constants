# A dig, 11 Nov 2026: MARIE JAMESON'S CURVES. Her paper "A refinement of Ramanujan's congruences modulo powers of 7 and 11"
# (https://web.math.utk.edu/~mjameso2/research/03.pdf, read as text) proves that inside Watson's and Atkin's laws for 7 and
# 11 the step rule at a prime Q (24N − 1 multiplied by Q²) multiplies by a number fixed by counting the points of an
# elliptic curve over the field with Q elements:
#   E7:  y² + xy + y = x³ + x² − 4x + 5  (conductor 42)        E11: y² + xy + y = x³ − 14x + 20  (conductor 726)
# Lovejoy and Ono (Crelle 542, 2002) proved the same for 5, where 1 + Q stands in place of the count. Here: their theorems
# checked on every value we have; the primes that give laws one factor higher (her Corollaries 1.2 and 1.3) set beside the
# ones counted in partitions_tower.txt; the share of all primes that give them; and new laws from primes beyond the
# recurrence, checked with exact values of p(N) (Hardy–Ramanujan–Rademacher formula, sympy).
# Formulas QUOTED from the papers; checks, members and shares COUNTED; the matrix counts DERIVED (counted in GL2).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_jameson.py [--keep]
import os, sys, time
import numpy as np
import sympy as sp
from sympy.functions.combinatorial.numbers import partition
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_jameson.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); NMAX = 2_400_000; L = 5 ** 6 * 7 ** 4 * 11 ** 3
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
say('p(n) for every n up to %d on the lane 5⁶·7⁴·11³ (%.0f s); p(100) right: %s' % (NMAX, time.time() - t0, int(p[100]) == 190569292 % L))

CURVE = {7: (1, 1, 1, -4, 5), 11: (1, 0, 1, -14, 20)}


def trace(curve, Q):
    """1 + Q − #E(F_Q), counted: for each x the quadratic in y has 1 + (D(x)|Q) roots"""
    a1, a2, a3, a4, a6 = curve; xs = np.arange(Q, dtype=np.int64); x2 = xs * xs % Q; x3 = x2 * xs % Q
    D = ((a1 * xs + a3) ** 2 + 4 * (x3 + a2 * x2 + a4 * xs + a6)) % Q
    sq = np.zeros(Q, dtype=bool); sq[x2] = True
    return -int(np.where(D == 0, 0, np.where(sq[D], 1, -1)).sum())


def disc(a1, a2, a3, a4, a6):
    b2, b4, b6 = a1 * a1 + 4 * a2, 2 * a4 + a1 * a3, a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6


say(); say('1  THE TWO CURVES')
for l in (7, 11):
    d = disc(*CURVE[l])
    say('   E%d %s: discriminant %d = %s, so its bad primes are %s;  1 + Q − #E(F_Q) for Q = 5 to 47: %s' % (
        l, list(CURVE[l]), d, sp.factorint(d), sorted(sp.factorint(abs(d))),
        ', '.join('%d: %+d' % (Q, trace(CURVE[l], Q)) for Q in sp.primerange(5, 48) if Q != l)))

leg = lambda a, Q: int(sp.jacobi_symbol(int(a) % Q, Q))


def step(l, j, Q, m_sign):
    """the multiplier of p(N*) in the theorems, given the sign of m = (24N* − 1)/l^j at Q (m_sign may be 0)"""
    if l == 5:
        return sp.jacobi_symbol(15, Q) * (1 + Q - (Q * Q * leg(-1, Q) * m_sign if j % 2 else leg(-1, Q) * m_sign))
    cQ = trace(CURVE[l], Q); e = 1 if l == 7 else 3; w = 3 if l == 7 else 5
    return Q ** e * leg(3, Q) * (cQ - (m_sign if j % 2 else Q ** w * leg(-1, Q) * m_sign))


THIRD = {5: 1, 7: 3, 11: 7}
R = {5: (19, 23), 7: (17, 23), 11: (13, 23)}                       # m ≡ r (mod 24) at odd / even levels
VLAW = {5: lambda j: j, 7: lambda j: j // 2 + 1, 11: lambda j: j}

say(); say('2  THEIR THEOREMS ON EVERY VALUE WE HAVE: p(N) ≡ [multiplier]·p(N*) − Q^e·p(N**) (mod ℓ^(law + 1)), where')
say('   24N − 1 = Q²(24N* − 1), 24N** − 1 = (24N* − 1)/Q², N* runs over the law\'s progression ℓʲn + β(j), and the multiplier')
say('   is Lovejoy–Ono\'s (15|Q)(1 + Q − …) for 5 and Jameson\'s Q^e·(3|Q)·(1 + Q − #E(F_Q) − …) for 7 and 11')
for l, js in ((5, (1, 2, 3, 4, 5)), (7, (1, 2, 3, 4, 5)), (11, (1, 2))):
    for j in js:
        mod = l ** (VLAW[l](j) + 1); r = R[l][(j + 1) % 2]; A = l ** j; nQ = checks = wrong = 0
        beta = (r * A + 1) // 24
        for Q in sp.primerange(7 if l == 5 else 5, 10 ** 6):
            if Q == l: continue
            if Q * Q * beta - (Q * Q - 1) // 24 > NMAX: break
            ns = np.arange(0, (NMAX - (Q * Q * beta - (Q * Q - 1) // 24)) // (A * Q * Q) + 1, dtype=np.int64)
            Ns = A * ns + beta; Nb = Q * Q * Ns - (Q * Q - 1) // 24; m = 24 * ns + r
            sg = np.array([leg(x, Q) for x in (m % Q).tolist()])
            mult = {s_: int(step(l, j, Q, s_)) % mod for s_ in (-1, 0, 1)}
            lhs = p[Nb] % mod
            rhs = (np.array([mult[s_] for s_ in sg.tolist()], dtype=np.int64) * (p[Ns] % mod)) % mod
            div = (m % (Q * Q) == 0)
            if div.any():
                Nss = (A * (m[div] // (Q * Q)) + 1) // 24
                rhs[div] = (rhs[div] - pow(Q, THIRD[l], mod) * (p[Nss] % mod)) % mod
            nQ += 1; checks += len(ns); wrong += int((lhs != rhs).sum())
        say('   ℓ = %2d, j = %d (mod %d^%d): %d primes Q, %d checks, %d wrong' % (l, j, l, VLAW[l](j) + 1, nQ, checks, wrong))

say(); say('3  THE PRIMES THAT GIVE LAWS ONE FACTOR HIGHER: 1 + Q − #E(F_Q) ≡ ±1 (mod ℓ) gives her Corollary 1.3 (progressions')
say('   modulo ℓʲ·Q³: our "μ = 0"), ≡ 0 gives her Corollary 1.2 (progressions modulo ℓʲ·Q⁴)')
HERS = {(7, 0): '31, 47, 79, 97, 113, 127, 191', (7, 1): '23, 41, 71, 89, 103, 131, 173, 199', (7, -1): '13, 53, 61, 149, 151, 163, 167',
        (11, 0): '11, 19, 23, 41, 59, 107, 193', (11, 1): '67, 109, 131, 139, 149, 179', (11, -1): '5, 101, 137, 167, 181'}
for l in (7, 11):
    tr = {Q: trace(CURVE[l], Q) % l for Q in sp.primerange(5, 200) if Q != l}
    for i in (0, 1, -1):
        ours = [Q for Q in tr if tr[Q] == i % l]
        say('   ℓ = %2d, ≡ %+d: counted %s;  her table: %s' % (l, i, ', '.join(map(str, ours)), HERS[(l, i)]))
say('   partitions_tower.txt (section 2) counted μ = 0 at Q = 13, 23, 41, 53, 61, 71, 89 for 7 and at Q = 5, 67 for 11: the ≡ ±1 primes')
say('   below 90 and 70 above. The 1375 law, p(1375n + 149) and p(1375n + 424) ≡ 0 (mod 121), is her Corollary 1.3 with Q = 5.')

say(); say('4  WHAT SHARE OF ALL PRIMES GIVES SUCH LAWS. Counted over the primes 5 ≤ Q ≤ 30,000 (Q ≠ ℓ); her counts for the')
say('   same range QUOTED; the limit DERIVED from her argument by counting the matrices of GL2 with each trace (Chebotarev)')
qs = [Q for Q in sp.primerange(5, 30001)]
for l in (7, 11):
    tr = np.array([trace(CURVE[l], Q) % l for Q in qs if Q != l]); nq = len(tr)
    g = [(a_, b_, c_, d_) for a_ in range(l) for b_ in range(l) for c_ in range(l) for d_ in range(l) if (a_ * d_ - b_ * c_) % l]
    bytr = {i: sum(1 for m4 in g if (m4[0] + m4[3]) % l == i % l) for i in (0, 1, -1)}
    hers = {7: ('0.1433', '0.1458', '0.1529'), 11: ('0.0860', '0.0878', '0.0924')}[l]
    say('   ℓ = %2d: %d primes; trace ≡ 0: %.4f (hers %s; limit %d/%d = %.4f);  ≡ +1: %.4f (hers %s; %d/%d = %.4f);  ≡ −1: %.4f (hers %s; %.4f)' % (
        l, nq, (tr == 0).mean(), hers[0], bytr[0], len(g), bytr[0] / len(g), (tr == 1).mean(), hers[1], bytr[1], len(g), bytr[1] / len(g),
        (tr == l - 1).mean(), hers[2], bytr[-1] / len(g)))
    say('           together, a law one factor higher comes from %.1f%% of all primes (limit %.1f%%)' % (
        100 * ((tr == 0) | (tr == 1) | (tr == l - 1)).mean(), 100 * (bytr[0] + bytr[1] + bytr[-1]) / len(g)))
q5 = np.array([Q % 5 for Q in qs])
say('   ℓ =  5 (Lovejoy–Ono, no curve): Q ≡ 3 (mod 5) zeroes a class (%.4f of the primes), Q ≡ 4 zeroes both (%.4f); limit 1/4 each, 1/2 together' % (
    (q5 == 3).mean(), (q5 == 4).mean()))

say(); say('5  LAWS AT PRIMES WE HAD NOT TESTED, FAR BEYOND THE RECURRENCE (her corollaries with j = 1), checked with exact values of')
say('   p(N). Members: N = ℓ·Q²·n + β(1, Q) for the first 10 n from 1,000 on that satisfy the corollary\'s condition (every N above')
say('   2,400,000); control: the first 10 such n of the other class.')
say('   For 7 only n outside the square rule are used (p(7n + 5) is already a multiple of 49 when 24n + 17 is a square on lane 7)')


def law_n(l, Q, i, n):
    m = 24 * n + R[l][0]
    if l == 5: return leg(-m, Q) == 1                       # Lovejoy–Ono, Q ≡ 3 (mod 5): the class that their multiplier zeroes
    if i == 0: return m % Q == 0 and m % (Q * Q) != 0       # Corollary 1.2
    return leg(m, Q) == i                                   # Corollary 1.3


def ctrl_n(l, Q, i, n):
    m = 24 * n + R[l][0]
    if l == 5: return leg(-m, Q) == -1
    if i == 0: return m % Q != 0
    return leg(m, Q) == -i


for l, Q, i in [(7, 103, 1), (7, 149, -1), (7, 31, 0), (11, 101, -1), (11, 109, 1), (11, 19, 0), (5, 113, None)]:
    beta = (R[l][0] * l * Q * Q + 1) // 24; mod = l ** (VLAW[l](1) + 1); free = lambda n: l != 7 or leg(24 * n + 17, 7) == -1
    ns_law, ns_ctl, n = [], [], 1000
    while len(ns_law) < 10 or len(ns_ctl) < 10:
        if free(n):
            if law_n(l, Q, i, n) and len(ns_law) < 10: ns_law.append(n)
            elif ctrl_n(l, Q, i, n) and len(ns_ctl) < 10: ns_ctl.append(n)
        n += 1
    hit = sum(1 for n in ns_law if int(partition(l * Q * Q * n + beta)) % mod == 0)
    hc = sum(1 for n in ns_ctl if int(partition(l * Q * Q * n + beta)) % mod == 0)
    cond = {None: '((−24n − 19)|%d) = +1' % Q, 0: '%d divides 24n + %d once' % (Q, R[l][0])}.get(i, '((24n + %d)|%d) = %+d' % (R[l][0], Q, i or 0))
    trq = 'Lovejoy–Ono, Q ≡ 3 (mod 5)' if l == 5 else '1 + Q − #E ≡ %+d' % (((trace(CURVE[l], Q) + l // 2) % l) - l // 2)
    say('   ℓ = %2d, Q = %3d (%s): p(%d·%d²·n + %d) with %s ≡ 0 (mod %d^%d): %d of %d (largest N = %d);  control: %d of %d' % (
        l, Q, trq, l, Q, beta, cond, l, VLAW[l](1) + 1, hit, len(ns_law), l * Q * Q * max(ns_law) + beta, hc, len(ns_ctl)))
say(); say('seconds: %.0f' % (time.time() - t0))
