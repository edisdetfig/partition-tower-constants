# first_step_ghost.py, 14 Nov 2026 (README "Why exactly 13, 17, 19 and 23 settle fast").
# The slow side needs the first-step families to have slope exactly 1 at the tower's weight −2. Their Serre weight is
# σ = Sym^(ℓ−7) ⊗ det, and Theorem 1.3 gives their slopes as those of the ghost series G_σ. Its zeros come from the counts
# d_z (level 6) and D_z (level 6ℓ) of a family with that one seed (multiplicity one), which are computed from GL2(F_ℓ):
# slope_rule.py's counting, repeated here.
# At w_(−2) (step x = −1): v(g_i) = Σ over the zeros of g_i of mult × (1 + v_ℓ(z + 1)). The first slope of the Newton
# polygon is v(g_1) when v(g_2) ≥ 2·v(g_1). Printed for every prime 29 … LMAX: v(g_1), v(g_2), the first slope, and the
# step where the family's second form arrives. DERIVED (representation theory and the ghost series).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/first_step_ghost.py [LMAX] [--keep]
import os, sys
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; LMAX = nums[0] if nums else 613; KEEP = '--keep' in sys.argv
OUT = os.path.join(HERE, 'first_step_ghost_%d.txt' % LMAX)
if KEEP and os.path.exists(OUT): sys.exit('stopped before computing anything: %s exists' % OUT)
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def is_prime(n): return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def vl(x, L):
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def jh_ind(s, t, LM):
    s, t = s % LM, t % LM
    if s == t: return [(0, t), (LM, t)]
    r = (s - t) % LM
    return [(r, t), (LM - r, s)]


def count_sym(n, sig, L):
    LM = L - 1; c = 0; j = 0
    while True:
        if n <= LM: return c + ((n, j % LM) == sig)
        c += sum(1 for g in jh_ind(n, 0, LM) if (g[0], (g[1] + j) % LM) == sig)
        if n < L + 1: return c
        n -= L + 1; j += 1


def count_ind(k, sig, L):
    LM = L - 1; n = k - 2
    full, rest = divmod(n + 1, LM)
    per = [sum(1 for g in jh_ind(n - i, i, LM) if g == sig) for i in range(LM)]
    return full * sum(per) + sum(per[:rest])


say('the first-step families (seed Sym^(ℓ−7) ⊗ det): the ghost series at the tower\'s weight −2, primes 29 … %d' % LMAX)
ok = total = 0; lines = []
for L in range(29, LMAX + 1):
    if not is_prime(L): continue
    sig = (L - 7, 1); zmax = L + 10
    d, D = [], []
    for z in range(zmax + 1):
        k = -2 + (L - 1) * (z + 1); d.append(count_sym(k - 2, sig, L)); D.append(count_ind(k, sig, L))
    def v(i):
        return sum(min(i - d[z], D[z] - d[z] - i) * (1 + vl(z + 1, L)) for z in range(zmax + 1) if d[z] < i < D[z] - d[z])
    v1, v2 = v(1), v(2)
    first = Fraction(v1) if v2 >= 2 * v1 else Fraction(v2, 2)
    second = next((z for z in range(zmax + 1) if d[z] >= 2), None)
    total += 1; ok += (first == 1)
    lines.append((L, d[:3], D[:3], v1, v2, first, second))
    if L <= 61 or first != 1:
        say('   ℓ = %3d: d = %s …, D = %s …; v(g_1) = %d, v(g_2) = %d; first slope at −2: %s; second form arrives at x = %s' % (
            L, d[:3], D[:3], v1, v2, first, second))
say('   first slope exactly 1 at %d of %d primes; v(g_2) ranges over %d … %d' % (
    ok, total, min(t[4] for t in lines), max(t[4] for t in lines)))
