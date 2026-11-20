# A deeper dig, 11 Nov 2026: THE LAWS FOR 7 THAT NO SOURCE STATES (every 2401st partition number from the 929th, 1615th
# and 1958th is a multiple of 343), what they say about Ramanujan's one wrong guess about partitions, why they hold, and
# whether the pattern goes on. p(n) is worked out by Euler's pentagonal rule with whole numbers only, kept on the lane
# 5⁷·7⁶·11³. Every number printed is COUNTED; where a sentence says "derived", it is reasoning, not a run.
#   1  p(n) for every n up to 2,400,000.
#   2  The laws, power by power, for 5 (to 5⁷), 7 (to 7⁶) and 11 (to 11³): only the classes new at each power.
#   3  The square rule's forecasts at every level where it should apply, the next ones included, and at the controls.
#   4  Why: at the level below each set of laws, the generating function, divided by the power that divides it, is on
#      lane ℓ a constant times E(q^ℓ)^a · E(q)³ (Jacobi's series, which skips whole classes). Checked coefficient by
#      coefficient for every n available.
#   5  Ramanujan's conjecture for 7³ (if 7³ divides 24N − 1 then 7³ divides p(N); false, Chowla 1934: N = 243), sorted by
#      the square rule; the two OEIS lists (A340957 where it holds, A340757 where it fails) made again from scratch.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_deeper.py [--keep]
import os, sys, time
import numpy as np
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_deeper.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); N = 2_400_000; POW = {5: 7, 7: 6, 11: 3}; L = 5 ** 7 * 7 ** 6 * 11 ** 3
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
say('1  p(n) for every n up to %d, kept on the lane 5⁷·7⁶·11³ = %d (%.0f s). p(100) comes out whole: %s' % (N, L, time.time() - t0, int(p[100]) == 190569292))


def zero_classes(A, M):
    out = {}
    for r in range(A):
        v = p[r::A]
        if len(v) and not (v % M).any(): out[r] = len(v)
    return out


say(); say('2  THE LAWS, POWER BY POWER (only the classes new at each power; "m tried" is how many members each class has here)')
found = {}
for l in (5, 7, 11):
    for a in range(1, POW[l] + 1):
        A = l ** a
        if N // A < 8: break
        for b in range(1, min(a, POW[l]) + 1):
            z = zero_classes(A, l ** b); found[(l, a, b)] = z
            below = found.get((l, a - 1, b), {}); new = {r: n for r, n in z.items() if (r % (A // l)) not in below}
            if new: say('   ℓ = %2d, every %6d-th from %-36s a multiple of %-5d (m tried: %s)' % (l, A, ', '.join(str(r) for r in sorted(new)), l ** b, min(new.values())))

say(); say('3  THE SQUARE RULE: among N with 24N ≡ 1 (mod A), p(N) gains one more factor ℓ exactly where (24N − 1)/A is a square on lane ℓ')
sq = lambda l: {(x * x) % l for x in range(l)}
for l, A, M, what in ((5, 25, 125, 'rule'), (5, 625, 3125, 'rule'), (5, 15625, 78125, 'rule, next level'), (7, 7, 49, 'rule'), (7, 343, 343, 'rule'),
                      (7, 16807, 2401, 'rule, next level'), (5, 125, 625, 'control'), (5, 3125, 15625, 'control'), (7, 49, 343, 'control'), (7, 2401, 2401, 'control')):
    r0 = next(r for r in range(A) if (24 * r - 1) % A == 0)
    forecast = sorted(r0 + j * A for j in range(l) if ((24 * (r0 + j * A) - 1) // A) % l in sq(l))
    z = zero_classes(A * l, M); got = sorted(r for r in z if r % A == r0); tried = min([z[r] for r in got], default=0)
    say('   %-17s ℓ = %d, from N ≡ %d (mod %d): forecast multiples of %d on every %d-th from %s;  found %s (m tried: %d)  %s' % (
        what + ':', l, r0, A, M, A * l, forecast, got, tried, 'the same' if forecast == got else 'not the same'))

say(); say('4  WHY: Σ p(A·n + r₀) qⁿ, divided by the power ℓᵛ that divides every coefficient, set beside c·E(q^ℓ)^a·E(q)³ on lane ℓ')


def euler(l, nmax):
    e = np.zeros(nmax + 1, dtype=np.int64); kk = 0
    while True:
        done = True
        for k in ((kk, -kk) if kk else (0,)):
            g = k * (3 * k - 1) // 2
            if g * l <= nmax: e[g * l] += -1 if k % 2 else 1; done = False
        if done and kk > 0: break
        kk += 1
    return e


def jacobi3(nmax):
    e = np.zeros(nmax + 1, dtype=np.int64); k = 0
    while k * (k + 1) // 2 <= nmax: e[k * (k + 1) // 2] += (-1) ** k * (2 * k + 1); k += 1
    return e


def mul(x, y, l, nmax): return np.convolve(x % l, y % l)[:nmax + 1] % l


for l, A, v in ((5, 25, 2), (5, 625, 4), (5, 15625, 6), (7, 7, 1), (7, 343, 2), (7, 16807, 3)):
    r0 = next(r for r in range(A) if (24 * r - 1) % A == 0); nmax = min((N - r0) // A, 3000)      # 3,000 coefficients at most: the series are multiplied the plain way
    g = p[r0::A][:nmax + 1] % (l ** (v + 1)); whole = not (g % l ** v).any(); H = (g // l ** v) % l
    E3 = jacobi3(nmax) % l; El = euler(l, nmax) % l; match = None
    for a in range(0, 3 * l):
        base = E3.copy()
        for _ in range(a): base = mul(base, El, l, nmax)
        cc = int(H[0] * pow(int(base[0]), -1, l)) % l if base[0] % l else None
        if cc and not ((H - cc * base) % l).any(): match = (a, cc); break
    say('   ℓ = %d, A = %5d, r₀ = %5d: every coefficient a multiple of %d: %s;  on lane %d it is %s, for all n ≤ %d' % (
        l, A, r0, l ** v, whole, l, ('%d · E(q^%d)^%d · E(q)³' % (match[1], l, match[0])) if match else 'NOT of that shape', nmax))
say('   (derived) E(q)³ = Σ (−1)ᵏ(2k + 1) q^(k(k+1)/2) has, on lane ℓ, nothing on the places where 8n + 1 is not a square, and on the')
say('   place where 8n + 1 is 0 its coefficient 2k + 1 is a multiple of ℓ. E(q^ℓ)^a only moves things by whole steps of ℓ. So the')
say('   classes where Jacobi\'s series is silent are exactly the classes of the laws above.')

say(); say('5  RAMANUJAN\'S CONJECTURE FOR 7³, SORTED BY THE SQUARE RULE')
say('   He guessed: if 7³ divides 24N − 1, then 7³ divides p(N). It is false (Chowla, 1934: N = 243). Watson (1938) proved 7² there.')


def v7(n):
    x = int(p[n]) % 7 ** 6; v = 0
    while v < 6 and x % 7 == 0: x //= 7; v += 1
    return v


def v7_of(m):
    v = 0
    while m % 7 == 0: m //= 7; v += 1
    return v


ks = range(243, N + 1, 343); rows = {}
for k in ks:
    j = ((k - 243) // 343) % 7; q7 = ((24 * k - 1) // 343) % 7; rows.setdefault((j, q7), []).append(k)
for (j, q7), members in sorted(rows.items(), key=lambda x: x[0][0]):
    kind = 'square on lane 7' if q7 in (1, 2, 4) else ('0: 7⁴ divides 24N − 1 (Watson: 7³)' if q7 == 0 else 'NOT a square on lane 7')
    hold = sum(1 for k in members if v7(k) >= 3); four = sum(1 for k in members if v7(k) >= 4)
    say('   N ≡ %4d (mod 2401): (24N − 1)/343 leaves %d on lane 7, %-38s 7³ divides p(N): %6d of %d;  7⁴: %d' % (243 + 343 * j, q7, kind + ';', hold, len(members), four))
holds = [k for k in ks if v7(k) >= v7_of(24 * k - 1)]; fails = [k for k in ks if v7(k) < v7_of(24 * k - 1) and v7(k) < 6]
A340957 = [929, 1615, 1958, 3330, 4016, 4359, 5731, 6417, 6760, 8132, 8818, 9161, 10533, 11219, 11562, 12934, 13620, 13963, 15335, 16021, 16364, 17736, 18422,
           18765, 20137, 20823, 21166, 22538, 23224, 23567, 24939, 25282, 25625, 25968, 27340, 28026, 28369, 29741, 30427]       # the terms OEIS shows (copied from the entry, for comparison)
A340757 = [243, 586, 1272, 2301, 2644, 2987, 3673, 4702, 5045, 5388, 6074, 7103, 7446, 7789, 8475, 9504, 9847, 10190, 10876, 11905, 12248, 12591, 13277, 14306, 14649,
           14992, 15678, 16707, 17050, 17393, 18079, 19108, 19451, 19794, 20480, 21509, 21852, 22195, 22881, 23910]                   # the same for A340757
say('   OEIS A340957 (where it holds), made again: the first %d terms the same as OEIS shows: %s' % (len(A340957), holds[:len(A340957)] == A340957))
say('   OEIS A340757 (where it fails), made again: the first %d terms the same as OEIS shows: %s' % (len(A340757), fails[:len(A340757)] == A340757))
law = {929, 1615, 1958}; odd = [k for k in A340957 if k % 2401 not in law]
say('   of the %d terms OEIS shows for A340957, from the three laws: %d; from elsewhere (by chance): %s' % (len(A340957), len(A340957) - len(odd), odd))
say(); say('seconds: %.0f' % (time.time() - t0))
