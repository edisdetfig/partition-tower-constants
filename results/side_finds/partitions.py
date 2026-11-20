# A dig, 11 Nov 2026: EVERY 125th PARTITION NUMBER FROM THE 74th AND THE 124th IS A MULTIPLE OF 125, and how far the
# pattern reaches. Found blind by the side study (LAWS2); it is in Ramanujan's unpublished manuscript (Berndt and Ono,
# Sém. Lothar. Combin. 42, 1999, eq. 21.3: p(125n−1), p(125n−26), p(125n−51) ≡ 0 mod 125).
#   1  p(n) for every n up to 400,000, kept on the lane 5⁵·7⁴·11² (Euler's pentagonal rule; whole numbers only).
#   2  For 5, 7 and 11 and each power A = ℓᵃ: every class r (mod A) on which p(A·m + r) is a multiple of ℓᵇ for every m
#      tried. Only the classes that are new at that power are shown (not just a class of the power below, carried up).
#   3  The rule, said once: on the classes the known laws leave open one level down, p(N) is a multiple of the next
#      power exactly where (24N − 1)/ℓᵏ is a square on lane ℓ (0 counts). Its forecasts, set beside what came out.
# Every number printed is COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions.py [--keep]
import os, sys, time
import numpy as np
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); N = 400_000; L = 5 ** 5 * 7 ** 4 * 11 ** 2
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
say('1  p(n) for n up to %d, kept on the lane %d (%.0f s). Check: p(100) on that lane = %d (p(100) = 190,569,292).' % (N, L, time.time() - t0, p[100]))
say('   (p(100) itself is smaller than the lane, so it must come out whole: %s)' % (int(p[100]) == 190569292))


def zero_classes(A, M):
    """Every class r mod A on which p(A·m + r) is a multiple of M for every m tried; with how many m were tried."""
    out = {}
    for r in range(A):
        v = p[r::A]
        if len(v) and not (v % M).any(): out[r] = len(v)
    return out


say(); say('2  THE LAWS, POWER BY POWER: the classes r (mod ℓᵃ) on which every p(ℓᵃ·m + r) tried is a multiple of ℓᵇ, new at that power')
found = {}
for l, top in ((5, 5), (7, 4), (11, 2)):
    for a in range(1, top + 1):
        A = l ** a
        for b in range(1, a + 1):
            z = zero_classes(A, l ** b); found[(l, a, b)] = z
            below = found.get((l, a - 1, b), {}); new = {r: n for r, n in z.items() if (r % (A // l)) not in below}
            if new: say('   ℓ = %2d, every %5d-th from %-28s a multiple of %d   (m tried: %s)' % (l, A, ', '.join(str(r) for r in sorted(new)), l ** b, min(new.values())))
say('   (%.0f s)' % (time.time() - t0))

say(); say('3  THE RULE: on the classes left open one level down, p(N) is a multiple of the next power exactly where (24N − 1)/ℓᵏ is a square on lane ℓ')
say('   (N runs over the class 24N ≡ 1 mod A, which the known law of that level covers; the forecast picks the next-level')
say('   classes where (24N − 1)/A is a square on lane ℓ. The rule is expected where Ramanujan\'s identity for that level comes')
say('   down to Jacobi\'s series, and not at the levels in between: those rows are a control.)')
sq = lambda l: {(x * x) % l for x in range(l)}
for l, A, M in ((5, 25, 125), (5, 125, 625), (5, 625, 3125), (7, 7, 49), (7, 49, 343), (7, 343, 343)):
    r0 = next(r for r in range(A) if (24 * r - 1) % A == 0)
    forecast = sorted(r0 + j * A for j in range(l) if ((24 * (r0 + j * A) - 1) // A) % l in sq(l))
    got = sorted(r for r in zero_classes(A * l, M) if r % A == r0)
    say('   ℓ = %d, from N ≡ %d (mod %d): forecast multiples of %d on every %d-th from %s;  found: %s   %s' % (
        l, r0, A, M, A * l, forecast, got, 'the same' if forecast == got else 'not the same'))
say(); say('seconds: %.0f' % (time.time() - t0))
