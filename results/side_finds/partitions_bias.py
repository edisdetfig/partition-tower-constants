# A dig, 11 Nov 2026: THE EXTRA ZEROS. Beyond every law (Ramanujan, Watson, the square rule), p(N) is a multiple of 5 for
# 20.5% of N, not 20%, and of 7 for 15.2%, not 14.3% (partitions_fractal.txt). Berndt and Ono's table (1999, n ≤ 100,000)
# shows the same and calls the densities "very little known". Here: how big the excess is for each prime from 2 to 31,
# whether it fades as n grows, where it lives (a scan of arithmetic progressions), and whether the fractal of
# partitions_fractal.py holds for 11 as well. p(n) for n up to 2,400,000 by Euler's pentagonal rule, whole numbers only,
# kept on the lane 11⁴·2·3·5·7·13·17·19·23·29·31. Shares, classes and members are COUNTED; z is DERIVED from the counts;
# Berndt and Ono's shares are QUOTED from their table.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_bias.py [--keep]
import os, sys, time, math
import numpy as np
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_bias.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); N = 2_400_000; PR = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]; L = 11 ** 3 * math.prod(PR)
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
say('p(n) for every n up to %d, kept on the lane 11⁴·2·3·5·7·13·17·19·23·29·31 (%.0f s); p(100) comes out whole: %s' % (N, time.time() - t0, int(p[100]) == 190569292))

say(); say('1  THE FRACTAL FOR 11, with the laws divided out (as partitions_fractal.txt shows for 5 and 7)')
h1 = (p[6::11] % 121) // 11 % 11; h3 = (p[721::1331] % 11 ** 4) // 1331 % 11; m = len(h3)
whole = not (p[6::11] % 11).any() and not (p[721::1331] % 1331).any()
lam = next((x for x in range(1, 11) if not ((h3 - x * h1[:m]) % 11).any()), None)
say('   every term a multiple of the law\'s power: %s;  p(1331·n + 721)/1331 ≡ %s · p(11·n + 6)/11 (mod 11), for every n ≤ %d' % (whole, lam if lam is not None else 'NO SINGLE CONSTANT', m - 1))

LAW = {5: (5, 4), 7: (7, 5), 11: (11, 6)}                      # the classes where p(n) is always a multiple of ℓ
BO = {2: 0.498, 3: 0.333, 5: 0.362, 7: 0.272, 11: 0.174, 13: 0.080}      # Berndt and Ono's shares of 0 for n ≤ 100,000 (copied from their table)
say(); say('2  THE SHARE OF ZEROS ON THE CLASSES WITHOUT A LAW, prime by prime, against 1/ℓ (z: how many standard errors above chance)')
blocks = 8; edges = np.linspace(0, N + 1, blocks + 1).astype(int)
for l in PR:
    n_all = np.arange(N + 1); free = np.ones(N + 1, dtype=bool)
    if l in LAW: free &= (n_all % LAW[l][0] != LAW[l][1])
    zero = (p % l == 0) & free; share = zero.sum() / free.sum(); z = (share - 1 / l) / math.sqrt((1 / l) * (1 - 1 / l) / free.sum())
    allshare = (p % l == 0).mean(); per = [zero[edges[i]:edges[i + 1]].sum() / max(1, free[edges[i]:edges[i + 1]].sum()) for i in range(blocks)]
    say('   ℓ = %2d: share %.5f against %.5f (z = %+6.1f); all n: %.4f (Berndt–Ono, n ≤ 100,000: %s);  by block of 300,000: %s' % (
        l, share, 1 / l, z, allshare, BO.get(l, '—'), ' '.join('%.4f' % x for x in per)))

say(); say('3  WHERE THE EXTRA ZEROS LIVE: every class B (mod A), A from 2 to 1,200, on the numbers without a law; the classes most above')
say('   chance, and any class where p is a multiple of ℓ for EVERY member (a law not counted so far)')
for l in (5, 7, 13):
    n_all = np.arange(N + 1); free = np.ones(N + 1, dtype=bool)
    if l in LAW: free &= (n_all % LAW[l][0] != LAW[l][1])
    zero = ((p % l == 0) & free).astype(np.int64); fr = free.astype(np.int64); best = []; laws = []
    for A in range(2, 1201):
        rows = (N + 1) // A; zc = zero[:rows * A].reshape(rows, A).sum(0); fc = fr[:rows * A].reshape(rows, A).sum(0)
        ok = fc >= 400
        sh = np.where(ok, zc / np.maximum(fc, 1), 0); zz = (sh - 1 / l) / np.sqrt((1 / l) * (1 - 1 / l) / np.maximum(fc, 1))
        for B in np.nonzero(ok & (zc == fc))[0]: laws.append((A, int(B), int(fc[B])))
        i = int(np.argmax(np.where(ok, zz, -1e9))); best.append((float(zz[i]), A, i, float(sh[i]), int(fc[i])))
    best.sort(reverse=True)
    say('   ℓ = %d: classes where every member is a multiple (beyond the known law): %s' % (l, laws[:10] if laws else 'none'))
    say('      most biased classes (z, A, B, share, members): ' + ';  '.join('%+.1f: %d·n + %d, %.3f of %d' % (zz, A, B, sh, fc) for zz, A, B, sh, fc in best[:8]))
say(); say('seconds: %.0f' % (time.time() - t0))
