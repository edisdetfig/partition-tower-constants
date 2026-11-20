# A dig, 11 Nov 2026: THE LEFTOVER. partitions_cells.txt: beyond every law, among the N whose 24N − 1 has no square factor at
# all (no Hecke rule can reach them), p(N) is a multiple of 5 for 20.64% of Atkin's class instead of 20%, and of 7 for 14.46%
# instead of 14.29%. Here: where exactly it shows (every prime 5 to 23, every class, and inside the laws once divided out),
# what it does not depend on, how deep the extra multiples go, and whether it is the class numbers of quadratic fields.
# "Atkin's class": −(24N − 1) is not a square on lane ℓ (Radu 2013: laws live only there or where ℓ divides 24N − 1).
# p(n) for every n up to 2,400,000 by Euler's pentagonal rule, whole numbers only, on the lane 5⁴·7⁴·11³·13·17·19·23.
# Shares and members COUNTED; z DERIVED from the counts; class numbers COUNTED (reduced forms).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_leftover.py [--keep]
import os, sys, time, math, collections
import numpy as np
import sympy as sp
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_leftover.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); N = 2_400_000; L = 5 ** 4 * 7 ** 4 * 11 ** 3 * 13 * 17 * 19 * 23
G, S, k = [], [], 1
while True:
    a, b = k * (3 * k - 1) // 2, k * (3 * k + 1) // 2
    if a > N: break
    s = 1 if k % 2 else -1; G.append(a); S.append(s)
    if b <= N: G.append(b); S.append(s)
    k += 1
o = np.argsort(G); G = np.array(G, dtype=np.int64)[o]; S = np.array(S, dtype=np.int64)[o]
p = np.zeros(N + 1, dtype=np.int64); p[0] = 1; c = 0
for n in range(1, N + 1):
    while c < len(G) and G[c] <= n: c += 1
    p[n] = int(np.dot(S[:c], p[n - G[:c]])) % L
say('p(n) for every n up to %d on the lane 5⁴·7⁴·11³·13·17·19·23 (%.0f s); p(100) right: %s' % (N, time.time() - t0, int(p[100]) == 190569292 % L))

nn = np.arange(N + 1, dtype=np.int64); M = 24 * nn - 1; QMAX = math.isqrt(24 * N - 1); PR = list(sp.primerange(5, QMAX + 1))
sqf = np.ones(N + 1, dtype=bool); omega = np.zeros(N + 1, dtype=np.int64); rem = M.copy(); has = {l: np.zeros((N + 1, l), dtype=bool) for l in (5, 7)}
sqf_but = {l: np.ones(N + 1, dtype=bool) for l in (5, 7, 11)}          # no square factor other than a power of ℓ
for q in PR:
    r = pow(24, -1, q); idx = np.arange(r, N + 1, q); omega[idx] += 1; rem[idx] //= q
    for l in (5, 7): has[l][idx, q % l] = True
    r2 = pow(24, -1, q * q); i2 = np.arange(r2, N + 1, q * q); sqf[i2] = False
    for l in sqf_but:
        if q != l: sqf_but[l][i2] = False
big = rem > 1; omega[big] += 1
for l in (5, 7): has[l][np.nonzero(big)[0], rem[big] % l] = True
sqf[0] = False


def tab(l): return np.array([0] + [1 if pow(r, (l - 1) // 2, l) == 1 else -1 for r in range(1, l)], dtype=np.int64)


def zline(z, sel, l):
    sh = z[sel].mean(); se = math.sqrt((1 / l) * (1 - 1 / l) / sel.sum()); return sh, (sh - 1 / l) / se


say(); say('1  WHERE IT SHOWS: N whose 24N − 1 has no square factor at all; share of p(N) that are multiples of ℓ (z: standard errors above 1/ℓ)')
for l in (5, 7, 11, 13, 17, 19, 23):
    t = tab(l); z = (p % l == 0); parts = []
    for name, sel in (('Atkin\'s class', t[(-M) % l] == -1), ('the other class', t[(-M) % l] == 1), ('ℓ divides 24N − 1 once', (M % l == 0) & (M % (l * l) != 0))):
        if l in (5, 7, 11) and name.startswith('ℓ'): continue
        sh, zz = zline(z, sel & sqf, l); parts.append('%s %.5f (z = %+.1f)' % (name, sh, zz))
    say('   ℓ = %2d:  %s' % (l, ';  '.join(parts)))
say('   inside the laws (p(N) divided by the law\'s power, on lane ℓ; N with ℓᵇ dividing 24N − 1 and M\' = (24N − 1)/ℓᵇ with no square factor):')
for l, levels in ((5, (1, 2)), (7, (1, 2)), (11, (1, 2))):
    t = tab(l)
    for b in levels:
        A = l ** b; r0 = pow(24, -1, A); Ns = np.arange(r0, N + 1, A, dtype=np.int64); seq = p[Ns]; v = 0
        while np.all(seq % l ** (v + 1) == 0): v += 1
        Mp = (24 * Ns - 1) // A; core = sqf_but[l][Ns] & (Mp % l != 0); cz = ((seq // l ** v) % l) == 0; parts = []
        for cl in (1, -1):
            s = core & (t[Mp % l] == cl); sh = cz[s].mean(); se = math.sqrt((1 / l) * (1 - 1 / l) / s.sum())
            parts.append('M\' %s on lane %d: %.5f of %d (%s)' % ('a square' if cl == 1 else 'not a square', l, sh, int(s.sum()), 'all: the square rule' if sh == 1 else 'z = %+.1f' % ((sh - 1 / l) / se)))
        say('   ℓ = %2d, level %d (÷%d^%d):  %s' % (l, b, l, v, ';  '.join(parts)))

say(); say('2  WHAT IT DOES NOT DEPEND ON (Atkin\'s class, 24N − 1 with no square factor)')
for l in (5, 7):
    t = tab(l); z = (p % l == 0); base = (t[(-M) % l] == -1) & sqf
    by_w = ', '.join('%d: %.4f' % (w, z[base & (omega == w)].mean()) for w in range(1, 6) if (base & (omega == w)).sum() > 3000)
    by_c = ', '.join('≡ %d: %.4f / %.4f' % (cc, z[base & has[l][:, cc]].mean(), z[base & ~has[l][:, cc]].mean()) for cc in range(1, l))
    say('   ℓ = %d: by the number of prime factors of 24N − 1: %s' % (l, by_w))
    say('          with / without a prime factor in each class on lane %d: %s' % (l, by_c))
    shares = []
    for key, mod in ((nn, 8), (nn, 9), (M, l * l)) + tuple((M, q) for q in sp.primerange(5, 32) if q != l):
        r = key % mod
        for a_ in range(mod):
            s = base & (r == a_)
            if s.sum() > 5000: shares.append(z[s].mean())
    say('          in each of %d residue classes (N on lanes 8, 9; 24N − 1 on lanes %d and 5 to 31): from %.4f to %.4f;  by size of N (to 100,000 / 1,000,000 / 2,400,000): %s' % (
        len(shares), l * l, min(shares), max(shares), ', '.join('%.4f' % z[base & (nn >= lo) & (nn < hi)].mean() for lo, hi in ((1, 100_000), (100_000, 1_000_000), (1_000_000, N + 1)))))

say(); say('3  HOW DEEP THE EXTRA MULTIPLES GO (Atkin\'s class, 24N − 1 with no square factor): share of p(N) with exactly ℓ⁰, ℓ¹, ℓ², and ℓ³ or more')
for l in (5, 7):
    t = tab(l); base = (t[(-M) % l] == -1) & sqf; x = p[base]; v = np.zeros(len(x), dtype=np.int64); y = x.copy(); alive = np.ones(len(x), dtype=bool)
    for _ in range(4):
        d = alive & (y % l == 0); v[d] += 1; y[d] //= l; alive = d
    got = [(v == 0).mean(), (v == 1).mean(), (v == 2).mean(), (v >= 3).mean()]; ch = [1 - 1 / l, (1 - 1 / l) / l, (1 - 1 / l) / l ** 2, 1 / l ** 3]
    say('   ℓ = %d: %s  (chance: %s);  of the multiples, the share with ℓ² or more: %.4f (chance %.4f)' % (
        l, ' '.join('%.5f' % g for g in got), ' '.join('%.5f' % g for g in ch), (got[2] + got[3]) / (1 - got[0]), 1 / l))


def h(D):
    D = int(D); cnt = 0
    for a_ in range(1, math.isqrt(-D // 3) + 1):
        bs = np.arange(-a_ + 1, a_ + 1); bs = bs[(bs - D) % 2 == 0]; num = bs * bs - D; ok = num % (4 * a_) == 0
        bs = bs[ok]; cs = num[ok] // (4 * a_); keep = (cs >= a_) & ~((cs == a_) & (bs < 0)); bs = bs[keep]; cs = cs[keep]
        if len(bs): cnt += int((np.gcd(np.gcd(a_, np.abs(bs)), cs) == 1).sum())
    return cnt


say(); say('4  IS IT THE CLASS NUMBERS? For N ≤ 3,000 in Atkin\'s class for 5 with 24N − 1 = M with no square factor: p(N) on lane 5 against')
say('   the class number h of four quadratic fields. One fixed ratio (per class of M on lane 5) would mean p(N) ≡ C·h.')
assert [h(-23), h(-47), h(-71), h(-167), h(-20), h(-56)] == [3, 5, 7, 11, 2, 4]
rows = [(n_, 24 * n_ - 1, int(p[n_] % 5)) for n_ in range(1, 3001) if sqf[n_] and tab(5)[(-(24 * n_ - 1)) % 5] == -1]
for name, D in (('Q(√−M)', lambda m: -m), ('Q(√−5M)', lambda m: -5 * m), ('Q(√−3M)', lambda m: -12 * m), ('Q(√−15M)', lambda m: -60 * m)):
    hs = [h(D(m)) % 5 for _, m, _ in rows]; cnt = collections.Counter(); both = h0 = 0
    for (n_, m, pv), hv in zip(rows, hs):
        if hv: cnt[(pv * pow(hv, -1, 5)) % 5] += 1
        else: h0 += 1; both += pv == 0
    say('   %-9s %d numbers: p/h on lane 5 takes the values 0..4 %s times;  5 divides h for %d of them, and p(N) for %d of those' % (
        name, len(rows), [cnt[i] for i in range(5)], h0, both))
say(); say('seconds: %.0f' % (time.time() - t0))
