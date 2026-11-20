# A dig, 11 Nov 2026: THE PARTITION LAWS FOR 5 AND 7, THREE MORE THINGS. p(n) for every n up to 2,400,000 by Euler's
# pentagonal rule, whole numbers only, kept on the lane 5⁸·7⁶·11². Every number printed is COUNTED.
#   1  The fractal, with the laws divided out. Folsom, Kent and Ono (2012, Thm 1.1) show that two levels apart the partition
#      numbers repeat up to a constant (their example: p(13³n + 1007) ≡ 6·p(13n + 6) mod 13); for 5, 7 and 11 their constant
#      is 0, because Ramanujan's laws make every term divisible. Here: divide by the power the laws give, then compare.
#   2  "The rest is chance". Watson's bound with the square rule on top gives, for every N, how many 5s and 7s p(N) must
#      have. How many more it has, against a coin with 1-in-5 (1-in-7) odds per extra factor.
#   3  Ramanujan's guesses for how often p(n) is a multiple of 5 and of 7 (his manuscript: "about 1/3", "about 1/4", from
#      n ≤ 200), against 2,400,000 values and against what laws-plus-chance predicts (9/25 and 13/49).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_fractal.py [--keep]
import os, sys, time, collections
import numpy as np
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_fractal.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); N = 2_400_000; E5, E7 = 8, 6; L = 5 ** E5 * 7 ** E7 * 11 ** 2
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
say('p(n) for every n up to %d, kept on the lane 5⁸·7⁶·11² (%.0f s); p(100) comes out whole: %s' % (N, time.time() - t0, int(p[100]) == 190569292))


def vl(x, l, cap):
    v = 0
    while v < cap and x % l == 0: x //= l; v += 1
    return v


say(); say('1  THE FRACTAL WITH THE LAWS DIVIDED OUT: Hᵦ(n) = p(ℓᵇ·n + δ) / ℓ^(laws) on lane ℓ, set beside H at the level two below')
for l, cap, rows in ((5, E5, [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (7, 7)]), (7, E7, [(1, 1), (2, 2), (3, 2), (4, 3), (5, 3), (6, 4)])):
    H = {}
    for bb, v in rows:
        A = l ** bb
        if v + 1 > cap: continue
        r0 = next(r for r in range(A) if (24 * r - 1) % A == 0); seq = p[r0::A] % (l ** (v + 1))
        if (seq % l ** v).any(): say('   ℓ = %d, b = %d: NOT every term a multiple of %d' % (l, bb, l ** v)); continue
        H[bb] = ((seq // l ** v) % l, r0, v)
    for bb in sorted(H):
        if bb - 2 not in H: continue
        h2, r2, v2 = H[bb]; h1, r1, v1 = H[bb - 2]; m = min(len(h1), len(h2))
        lam = next((x for x in range(1, l) if not ((h2[:m] - x * h1[:m]) % l).any()), None)
        say('   ℓ = %d: p(%d·n + %d)/%d ≡ %s · p(%d·n + %d)/%d  (mod %d), for every n ≤ %d' % (
            l, l ** bb, r2, l ** v2, ('%d' % lam if lam is not None else 'NO SINGLE CONSTANT'), l ** (bb - 2), r1, l ** v1, l, m - 1))

say(); say('2  "THE REST IS CHANCE": for every N ≤ %d, how many 5s (7s) p(N) has beyond what the laws force' % N)
say('   forced for 5: a = v₅(24N − 1) (Watson), plus 1 when a is even, at least 2, and (24N − 1)/5ᵃ is a square on lane 5 (the rule)')
say('   forced for 7: ⌊(b + 2)/2⌋ for b = v₇(24N − 1) ≥ 1 (Watson), plus 1 when b is odd and (24N − 1)/7ᵇ is a square on lane 7 (the rule)')
for l, cap in ((5, E5), (7, E7)):
    sq = {(x * x) % l for x in range(1, l)}; tab = collections.defaultdict(collections.Counter); skipped = 0
    for n in range(0, N + 1):
        m = 24 * n - 1; a = vl(m, l, 30) if m > 0 else 0; rest = m // l ** a
        if l == 5: forced = a + (1 if a >= 2 and a % 2 == 0 and rest % l in sq else 0)
        else: forced = (a + 2) // 2 + (1 if a % 2 == 1 and rest % l in sq else 0) if a >= 1 else 0
        if forced + 2 > cap: skipped += 1; continue
        x = int(p[n]) % l ** cap; v = vl(x, l, cap)
        if v < forced: tab['A LAW BROKEN'][a] += 1; continue
        e = min(v - forced, 2); tab[min(a, 3)][e] += 1
    say('   ℓ = %d   (N with too many forced factors to measure here: %d)' % (l, skipped))
    for a in sorted(k_ for k_ in tab if k_ != 'A LAW BROKEN'):
        cnt = tab[a]; tot = sum(cnt.values())
        say('      %s: %8d numbers;  no extra factor %.4f, one or more %.4f (chance: %.4f), two or more %.4f (chance: %.4f)' % (
            'v(24N − 1) = %d' % a if a < 3 else 'v(24N − 1) ≥ 3', tot, cnt[0] / tot, (cnt[1] + cnt[2]) / tot, 1 / l, cnt[2] / tot, 1 / l ** 2))
    say('      numbers where p(N) has FEWER factors than forced (a law broken): %d' % sum(tab['A LAW BROKEN'].values()))

say(); say('3  RAMANUJAN\'S GUESSES FOR HOW OFTEN p(n) IS A MULTIPLE OF 5 AND OF 7')
for l, guess, said, model in ((5, '69, 33, 34, 34, 30', '"about 1/3", the others "about 1/6" each', (9 / 25, 4 / 25)), (7, '50, 33, 22, 28, 23, 23, 21', '"about 1/4", the others "about 1/8" each', (13 / 49, 6 / 49))):
    small = collections.Counter(int(p[n]) % l for n in range(1, 201)); big = np.bincount(p[1:] % l, minlength=l) / N
    say('   ℓ = %d: n ≤ 200 (his count: %s): %s' % (l, guess, ', '.join(str(small[r]) for r in range(l))))
    say('      he guessed: multiples %s' % said)
    say('      n ≤ %d: multiples %.4f; the other remainders %s' % (N, big[0], ', '.join('%.4f' % x for x in big[1:])))
    say('      laws plus chance predict: multiples %.4f (= %s), each other remainder %.4f' % (model[0], '9/25' if l == 5 else '13/49', model[1]))
say(); say('seconds: %.0f' % (time.time() - t0))
