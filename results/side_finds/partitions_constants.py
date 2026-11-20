# A dig, 11 Nov 2026: ONE NUMBER PER PRIME. partitions_tower.txt (section 1) found that, with the laws divided out, partition
# numbers two levels apart differ by a constant that settles digit by digit into one ℓ-adic number per prime. Here:
#   1  what the data up to 2,400,000 predict for levels far beyond them, checked against exact values of p(N)
#      (Hardy–Ramanujan–Rademacher formula, sympy), the largest with about 425,000 digits;
#   2  the constants to many more digits, and how fast they settle;
#   3  the same constants for 17, 19, 23, 29, 31;
#   4  whether they are simple numbers (fractions, roots of small polynomials, unit roots of Hecke polynomials).
# Level b: H_b(n) = p(ℓᵇn + δ_b)/ℓ^v, δ_b = 1/24 on lane ℓᵇ, v = Watson's or Atkin's exponent (0 for ℓ ≥ 13). Digits are
# written lowest first (4, 2, 3 means 4 + 2ℓ + 3ℓ² + …). Digits COUNTED; predictions DERIVED from the data, then compared;
# relation sizes COUNTED (lattice reduction, LLL); first digits from the literature QUOTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_constants.py [--keep]
import os, sys, time, math, random
import numpy as np
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.functions.combinatorial.numbers import partition
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_constants.txt'), 'x', encoding='utf-8') if KEEP else None


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
P = {}
for l, e in ((5, 21), (7, 18), (11, 14), (13, 13)):
    L = l ** e; p = np.zeros(NMAX + 1, dtype=np.int64); p[0] = 1; c = 0
    for n in range(1, NMAX + 1):
        while c < len(G) and G[c] <= n: c += 1
        p[n] = int(np.dot(S[:c], p[n - G[:c]])) % L
    P[l] = (p, e)
say('p(n) for every n up to %d on the lanes 5²¹, 7¹⁸, 11¹⁴, 13¹³ (%.0f s)' % (NMAX, time.time() - t0))
random.seed(1); sample = [NMAX] + random.sample(range(1, NMAX), 200)
say('exact values (sympy) against the recurrence on lane 5²¹ for %d values of n up to %d: %d differ' % (
    len(sample), NMAX, sum(1 for n in sample if int(partition(n)) % 5 ** 21 != int(P[5][0][n]))))

VLAW = lambda l, b: b if l in (5, 11) else (b // 2 + 1 if l == 7 else 0)
dig = lambda x, l, k: [(x // l ** i) % l for i in range(k)]
show = lambda ds: ', '.join(map(str, ds))


def agree(a, b, l, k=80):
    n = 0
    while n < k and (a - b) % l ** (n + 1) == 0: n += 1
    return n


say(); say('1  PREDICTIONS FROM THE DATA ALONE, CHECKED WITH EXACT VALUES')
say('   From the data: for each pair of levels, the constant λ with H_(b+2)(n) ≡ λ·H_b(n) for every n, to as many digits as one λ')
say('   allows; "settled" = the digits on which the last two pairs agree. Prediction: H_b(0) ≡ λ^((b − b₀)/2)·H_b₀(0) to the settled')
say('   digits, from the deepest data level b₀ of the same parity.')
SETTLED, BASE = {}, {}
for l in (5, 7, 11, 13):
    p, e = P[l]; Hd = {}
    for b in range(1, 12):
        A = l ** b; r0 = pow(24, -1, A)
        if r0 > NMAX: break
        Ns = np.arange(r0, NMAX + 1, A, dtype=np.int64); v = VLAW(l, b)
        Hd[b] = (p[Ns] // l ** v, e - v)
    pairs = []
    for b in sorted(Hd):
        if b + 2 not in Hd or len(Hd[b + 2][0]) < 2: continue
        a1, a2 = Hd[b][0][:len(Hd[b + 2][0])], Hd[b + 2][0]; prec = min(Hd[b][1], Hd[b + 2][1])
        n0 = int(np.nonzero(a1 % l != 0)[0][0]); lam, best = None, 0
        for j in range(1, prec + 1):
            mod = l ** j; cand = int(a2[n0]) * pow(int(a1[n0]) % mod, -1, mod) % mod
            if np.all((a2 - cand * a1) % mod == 0): lam, best = cand, j
            else: break
        pairs.append((b, lam, best, len(a2)))
    (b1, l1, j1, _), (b2, l2, j2, _) = pairs[-2], pairs[-1]; K = min(agree(l1, l2, l), j1, j2)
    SETTLED[l] = (l2 % l ** K, K); tops = sorted(Hd)[-2:]; BASE[l] = {b: int(Hd[b][0][0]) for b in tops}
    say('   ℓ = %2d: pairs %s;  settled: %d digits: %s;  deepest data levels with n = 0: %s' % (
        l, ', '.join('(%d,%d) %d digits' % (b, b + 2, j) for b, _, j, _ in pairs), K, show(dig(l2, l, K)), tops))

TOP = {5: 16, 7: 13, 11: 10, 13: 9, 17: 8, 19: 8, 23: 7, 29: 7, 31: 7}
if '--quick' in sys.argv: TOP = {5: 12, 7: 10, 11: 8, 13: 7, 17: 6, 19: 6, 23: 6, 29: 5, 31: 5}     # a fast trial run only
EX = {}


def exact_levels(l):
    out = {}
    for b in range(1, TOP[l] + 1):
        d = pow(24, -1, l ** b); rb = (24 * d - 1) // l ** b; n = 0
        while (24 * n + rb) % l == 0: n += 1                                  # stay at exact level b (only 23 needs n = 1)
        N = l ** b * n + d; t = time.time(); x = int(partition(N)); v = VLAW(l, b)
        assert x % l ** v == 0, (l, b)
        out[b] = (x // l ** v, N, n, len(str(x)), time.time() - t)
    return out


for l in (5, 7, 11, 13):
    EX[l] = exact_levels(l); alpha, K = SETTLED[l]; rows = []
    for b in sorted(EX[l]):
        if b <= max(BASE[l]): continue
        b0 = max(bb for bb in BASE[l] if bb % 2 == b % 2); pred = pow(alpha, (b - b0) // 2, l ** K) * BASE[l][b0] % l ** K
        h, N, n, nd, sec = EX[l][b]; ok = agree(pred, h, l, K)
        rows.append('b = %d (N = %d, %d digits, %.0f s): %d of %d' % (b, N, nd, sec, ok, K))
    say('   ℓ = %2d, the %d settled digits predicted at each deeper level, right:  %s' % (l, K, ';  '.join(rows)))
say('   every exact value is also divisible by the law\'s power (checked)')

say(); say('2  THE CONSTANTS FROM THE EXACT VALUES: r_b = H_(b+2)(0)/H_b(0) as an ℓ-adic number; how many digits r_b shares with r_(b+1)')
ALPHA = {}


def chain(l):
    H = EX[l]; r = {}
    for b in sorted(H):
        if b + 2 in H and H[b][0] % l: r[b] = H[b + 2][0] * pow(H[b][0], -1, l ** 60) % l ** 60
    bs = sorted(r); ag = [agree(r[b], r[b + 1], l) for b in bs if b + 1 in r]
    K = ag[-1]; ALPHA[l] = (r[bs[-1]] % l ** K, K)
    return ag, K, r[bs[-1]]


for l in (5, 7, 11, 13):
    ag, K, last = chain(l)
    say('   ℓ = %2d: shared digits %s;  settled: %d digits: %s' % (l, ag, K, show(dig(last, l, K))))
say('   how fast: digits gained per two levels, from the last shared counts: %s' % ', '.join(
    '%d: %s' % (l, (lambda g: g[-1] - g[-3])(chain(l)[0])) for l in (5, 7, 11, 13)))

say(); say('3  THE SAME FOR 17, 19, 23, 29, 31 (no laws; level b uses n = 0, except 23 at even levels, where n = 0 lies one level deeper)')
for l in (17, 19, 23, 29, 31):
    EX[l] = exact_levels(l); ag, K, last = chain(l)
    say('   ℓ = %2d: levels 1 to %d (largest N = %d); shared digits %s;  settled: %d digits: %s' % (
        l, TOP[l], max(v_[1] for v_ in EX[l].values()), ag, K, show(dig(last, l, K))))
say('   first digits in the literature (QUOTED): 5 for 7 and 4·5 ≡ 9 for 11 (the factors per level in Watson\'s and Atkin\'s')
say('   generating functions, as Jameson\'s Cor. 2.2 and 2.3 state them: 5^⌊j/2⌋ and 4^(3(j−1)/2)); 6 and 45 = 6 + 3·13 for 13 (Folsom,')
say('   Kent and Ono, A₁₃(1,3,1) and A₁₃(2,4,2)); −1 = 4 for 5 is section 2 of partitions_tower_link.txt. Beyond the first digit (and')
say('   45 for 13) I found none of these digits stated anywhere.')

say(); say('4  ARE THEY SIMPLE NUMBERS?  The smallest whole-number relation c₀ + c₁α + … + c_dα^d ≡ 0 (lane ℓ^K) that lattice')
say('   reduction (LLL) finds, against the same search on 50 random ℓ-adic numbers with as many digits (their smallest and median);')
say('   and unit roots of x² − a·x ± ℓ^k with |a| ≤ 2ℓ^(k/2) (Hecke polynomials of classical eigenforms of weight k + 1)')


def lll_relation(alpha, l, K, d):
    # The kept partitions_constants.txt stops here: sympy's DomainMatrix.lll failed an internal assertion on the 33-digit
    # constant of 5. Section 4 was then run on its own, with a plain exact LLL, as partitions_simple.py (.txt). This version
    # uses that same LLL, so a new run would not stop there (with --keep it refuses to start, since the kept file exists).
    from fractions import Fraction
    m = l ** K; B = [[m] + [0] * d] + [[(-pow(alpha, i, m)) % m] + [1 if j == i - 1 else 0 for j in range(d)] for i in range(1, d + 1)]
    n = len(B)

    def gso(B):
        Bs, mu = [], [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = sum(Fraction(a) * b for a, b in zip(B[i], Bs[j])) / sum(b * b for b in Bs[j])
                v = [a - mu[i][j] * b for a, b in zip(v, Bs[j])]
            Bs.append(v)
        return Bs, mu
    Bs, mu = gso(B); k_ = 1
    while k_ < n:
        for j in range(k_ - 1, -1, -1):
            q = round(mu[k_][j])
            if q: B[k_] = [a - q * b for a, b in zip(B[k_], B[j])]; Bs, mu = gso(B)
        if sum(x * x for x in Bs[k_]) >= (Fraction(3, 4) - mu[k_][k_ - 1] ** 2) * sum(x * x for x in Bs[k_ - 1]): k_ += 1
        else: B[k_], B[k_ - 1] = B[k_ - 1], B[k_]; Bs, mu = gso(B); k_ = max(k_ - 1, 1)
    best = min(B, key=lambda r_: max(abs(x) for x in r_))
    return max(abs(x) for x in best), best


for l in (5, 7, 13):
    alpha, K = ALPHA[l]
    parts = []; rng = random.Random(7)
    for d in (1, 2, 3, 4):
        sz, vec = lll_relation(alpha, l, K, d)
        rnd = sorted(lll_relation(rng.randrange(1, l ** K), l, K, d)[0] for _ in range(50))
        parts.append('degree %d: %d (random: smallest %d, median %d)' % (d, sz, rnd[0], rnd[25]))
    hits = 0
    for kk in range(1, K):
        for eps in (1, -1):
            a_ = (alpha + eps * l ** kk * pow(alpha, -1, l ** K)) % l ** K
            a_ = a_ - l ** K if a_ > l ** K // 2 else a_
            hits += abs(a_) <= 2 * l ** (kk / 2)
    say('   ℓ = %2d (%d digits): largest coefficient of the smallest relation: %s;  Hecke unit roots that fit: %d' % (l, K, ';  '.join(parts), hits))
say(); say('seconds: %.0f' % (time.time() - t0))
