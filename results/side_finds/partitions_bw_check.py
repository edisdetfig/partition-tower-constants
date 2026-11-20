# A dig, 12 Nov 2026: OUR CONSTANTS AGAINST WHAT BOYLAN–WEBB PRINTED. Boylan–Webb (Trans. AMS 365 (2013); QUOTED in
# research_notes/Partition tower literature check/primes_37_and_up.md) print, for ℓ = 53:
#     p(53n+42) ≡ 22 p(53³n+117861) + 25 p(53⁵n+331071432)                    (mod 53),
#     p(53³n+117861) ≡ 2672 p(53⁵n+331071432) + 2304 p(53⁷n+929979652371)       (mod 53²),
# and the periods c₃₇ = 36, c₄₁ = 10 (the order of the matrix that moves the tower two levels, on lane ℓ).
# PREDICTED before the run (derived): our constants' polynomial x² + c₁x + c₀ (partitions_flint_family_53.txt, from modular
# forms only) is their relation scaled: c₀·(1, −22, −25) ≡ (c₀, c₁, 1) on lane 53, and c₀·(1, −2672, −2304) ≡ (c₀, c₁, 1) on
# lane 53²; and the period is the order of the companion matrix of our polynomial on lane ℓ: 36 for 37, 10 for 41.
# COUNTED.   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_bw_check.py [--keep]
import os, sys, re, glob
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_bw_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def poly_of(L):
    best = None
    for f in glob.glob(os.path.join(HERE, 'partitions_flint_family_%d*.txt' % L)):
        mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\]', open(f, encoding='utf-8').read())
        if mt and (best is None or int(mt.group(1)) > best[0]): best = (int(mt.group(1)), [int(v) for v in mt.group(2).split(',')], os.path.basename(f))
    return best


K, poly, name = poly_of(53)
say('53: our polynomial (lowest first) from %s, on lane 53^%d' % (name, K))
for e, bw in ((1, (1, -22, -25)), (2, (1, -2672, -2304))):
    m = 53 ** e; ours = [c % m for c in poly]; scaled = [ours[0] * b % m for b in bw]
    say('   lane 53^%d: ours %s;  Boylan–Webb\'s relation times c₀: %s  →  %s' % (e, ours, scaled, 'the same' if scaled == ours else 'DIFFERENT'))
for L, c in ((37, 36), (41, 10)):
    K, poly, name = poly_of(L); c0, c1 = poly[0] % L, poly[1] % L
    A = [[(-c1) % L, (-c0) % L], [1, 0]]; P = [[1, 0], [0, 1]]; order = None
    for t in range(1, L * L):
        P = [[sum(P[i][l] * A[l][j] for l in range(2)) % L for j in range(2)] for i in range(2)]
        if P == [[1, 0], [0, 1]]: order = t; break
    say('%d: our polynomial on lane %d: x² + %d·x + %d (%s); the order of its companion matrix: %d;  Boylan–Webb\'s period c_%d = %d  →  %s' % (
        L, L, c1, c0, name, order, L, c, 'the same' if order == c else 'DIFFERENT'))
