# A dig, 11 Nov 2026: ARE THE CONSTANTS SIMPLE NUMBERS? Section 4 of partitions_constants.py, run on its own:
# that run stopped at section 4, because sympy's lattice reduction (DomainMatrix.lll) failed an internal assertion on
# the 33-digit constant of 5; sections 1 to 3 of partitions_constants.txt are complete. This script reads the settled digits
# printed there (the one input; read only) and uses a plain exact LLL (fractions, written out below) instead of sympy's.
# Tests: (1) the smallest whole-number relation c₀ + c₁α + … + c_dα^d ≡ 0 on lane ℓ^K that LLL finds, against the same
# search on 100 random ℓ-adic units with as many digits; (2) unit roots of x² − a·x ± ℓ^k with |a| ≤ 2ℓ^(k/2) (Hecke
# polynomials of classical eigenforms of weight k + 1). Sizes COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_simple.py [--keep]
import os, re, sys, time, random
from fractions import Fraction
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_simple.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def lll(B, delta=Fraction(3, 4)):
    """textbook LLL in exact rational arithmetic (rows of B are the basis)"""
    B = [list(r) for r in B]; n = len(B)

    def gso(B):
        Bs, mu = [], [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = sum(Fraction(a) * b for a, b in zip(B[i], Bs[j])) / sum(b * b for b in Bs[j])
                v = [a - mu[i][j] * b for a, b in zip(v, Bs[j])]
            Bs.append(v)
        return Bs, mu
    Bs, mu = gso(B); k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q: B[k] = [a - q * b for a, b in zip(B[k], B[j])]; Bs, mu = gso(B)
        if sum(x * x for x in Bs[k]) >= (delta - mu[k][k - 1] ** 2) * sum(x * x for x in Bs[k - 1]): k += 1
        else: B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu = gso(B); k = max(k - 1, 1)
    return B


def relation(alpha, l, K, d):
    m = l ** K; rows = [[m] + [0] * d] + [[(-pow(alpha, i, m)) % m] + [1 if j == i - 1 else 0 for j in range(d)] for i in range(1, d + 1)]
    red = lll(rows); best = min(red, key=lambda r_: max(abs(x) for x in r_))
    assert sum(c_ * pow(alpha, i, m) for i, c_ in enumerate(best)) % m == 0
    return max(abs(x) for x in best), best


t0 = time.time()
txt = open(os.path.join(HERE, 'partitions_constants.txt'), encoding='utf-8').read()
CONST = {int(l): [int(x) for x in ds.split(', ')] for l, ds in re.findall(r'ℓ = +(\d+): (?:levels .*?; )?shared digits \[[^\]]*\];  settled: \d+ digits: ([\d, ]+)', txt)}
say('settled digits read from partitions_constants.txt: %s' % ', '.join('%d: %d' % (l, len(v)) for l, v in sorted(CONST.items())))
say(); say('1  THE SMALLEST RELATION LLL FINDS (largest coefficient), against 100 random ℓ-adic units with as many digits:')
say('   the constant\'s value, then the random ones\' smallest / median / largest')
for l in sorted(CONST):
    ds = CONST[l]; K = len(ds)
    if K < 5: continue
    alpha = sum(d_ * l ** i for i, d_ in enumerate(ds)); rng = random.Random(7); parts = []
    for d in (1, 2, 3, 4):
        sz, vec = relation(alpha, l, K, d); rnd = []
        while len(rnd) < 100:
            x = rng.randrange(1, l ** K)
            if x % l: rnd.append(relation(x, l, K, d)[0])
        rnd.sort(); rank = sum(1 for r_ in rnd if r_ < sz)
        parts.append('degree %d: %d (%d / %d / %d; %d of 100 random below it)' % (d, sz, rnd[0], rnd[50], rnd[-1], rank))
    say('   ℓ = %2d (%d digits): %s' % (l, K, ';  '.join(parts)))
say(); say('2  HECKE UNIT ROOTS: for k = 1 to K − 1 and both signs, a = α ± ℓ^k/α on lane ℓ^K; a fit means |a| ≤ 2ℓ^(k/2)')
for l in sorted(CONST):
    ds = CONST[l]; K = len(ds)
    if K < 5: continue
    alpha = sum(d_ * l ** i for i, d_ in enumerate(ds)); m = l ** K; hits = []
    for kk in range(1, K):
        for eps in (1, -1):
            a_ = (alpha + eps * l ** kk * pow(alpha, -1, m)) % m; a_ = a_ - m if a_ > m // 2 else a_
            if abs(a_) <= 2 * l ** (kk / 2): hits.append((kk, eps, a_))
    say('   ℓ = %2d: %d fits%s' % (l, len(hits), '' if not hits else ': ' + str(hits)))
say(); say('seconds: %.0f' % (time.time() - t0))
