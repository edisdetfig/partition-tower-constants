# A dig, 13 Nov 2026: THE NON-ORDINARY ROOM IN WEIGHT ℓ − 1 (README "The speed of the towers"; the reading and its
# prediction are written there first). For every prime ℓ from 13 to 400, on lane ℓ: the cusp forms of level 1 and weight
# k = ℓ − 1, with the basis g_i = Δ^i·E4^a·E6^b (4a + 6b = k − 12i, i = 1 … d), which starts g_i = q^i + …; the matrix of
# T_ℓ (on lane ℓ, (T_ℓ g)(n) ≡ g(ℓn) since ℓ^(k−1) ≡ 0); the non-ordinary count = the dimension of T_ℓ's generalized
# kernel = d − rank(T_ℓ^d). PREDICTED: d − R_ℓ, with R_ℓ = ⌊(ℓ + 12)/24⌋. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/level1_ordinary_check.py [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
from flint import nmod_poly, nmod_mat
from sympy import primerange
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'level1_ordinary_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def dimS(k):
    if k < 12 or k % 2: return 0
    return k // 12 - (1 if k % 12 == 2 else 0)


t0 = time.time(); rows = []; bad = []
NMAXS = 400 * 40
SIG3 = [0] * NMAXS; SIG5 = [0] * NMAXS                     # σ₃ and σ₅ by a sieve
for dv in range(1, NMAXS):
    for mlt in range(dv, NMAXS, dv): SIG3[mlt] += dv ** 3; SIG5[mlt] += dv ** 5
for L in primerange(13, 401):
    k = L - 1; d = dimS(k); R = (L + 12) // 24
    Nl = L * d + 2; mul = lambda a, b: a.mul_low(b, Nl)
    s3 = [SIG3[n] % L for n in range(Nl)]; s5 = [SIG5[n] % L for n in range(Nl)]
    E4 = nmod_poly([1] + [240 * s3[n] % L for n in range(1, Nl)], L)
    E6 = nmod_poly([1] + [(-504 * s5[n]) % L for n in range(1, Nl)], L)
    eta = [0] * Nl; j = 0
    while True:                                               # ∏(1 − qⁿ), Euler's pentagonal series
        done = True
        for e, sgn in (((j * (3 * j - 1)) // 2, 1 if j % 2 == 0 else -1), ((j * (3 * j + 1)) // 2, 1 if j % 2 == 0 else -1)):
            if e < Nl: eta[e] = (eta[e] + sgn) % L; done = False
        if done: break
        j += 1
    eta[0] = 1
    Delta = mul(nmod_poly([0, 1], L), nmod_poly(eta, L).pow_trunc(24, Nl))
    G = []
    for i in range(1, d + 1):
        rest = k - 12 * i
        a = next(a for a in range(rest // 4, -1, -1) if (rest - 4 * a) % 6 == 0)
        b = (rest - 4 * a) // 6
        g = mul(mul(Delta.pow_trunc(i, Nl), E4.pow_trunc(a, Nl) if a else nmod_poly([1], L)), E6.pow_trunc(b, Nl) if b else nmod_poly([1], L))
        G.append(g)
    at = lambda p, i: int(p[i]) if i < p.length() else 0
    U = nmod_mat(d, d, [at(G[c], r + 1) for r in range(d) for c in range(d)], L)      # U[n][i] = g_i(n), n, i = 1 … d
    W = nmod_mat(d, d, [at(G[c], L * (r + 1)) for r in range(d) for c in range(d)], L)
    T = U.inv() * W
    P = T
    for _ in range(d): P = P * T
    ordinary = P.rank(); nono = d - ordinary
    ok = nono == d - R
    rows.append((L, d, R, ordinary, nono, ok))
    if not ok: bad.append(L)
say('Level 1, weight ℓ − 1, on lane ℓ: d = dim S_(ℓ−1), R = ⌊(ℓ + 12)/24⌋, the ordinary part (rank of T_ℓ^(d+1)) and the non-ordinary room')
for L, d, R, o, nn, ok in rows:
    say('   ℓ = %3d: d = %2d, R = %2d, ordinary %2d, not ordinary %2d  %s' % (L, d, R, o, nn, 'as PREDICTED (d − R)' if ok else 'NOT as predicted'))
say(); say('%d of %d primes as predicted%s;  seconds: %.0f' % (len(rows) - len(bad), len(rows), ('; not: %s' % bad) if bad else '', time.time() - t0))
