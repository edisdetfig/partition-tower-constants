# A dig, 11 Nov 2026: THE SURFACES FOR 17, 19, 23, CHECKED FURTHER. partitions_abelian_*.txt found, for each of these primes, a
# weight-2 newform with a quadratic coefficient field (an abelian surface over Q) at level 2·3·ℓ² whose numbers on lane ℓ are
# the tower's: a_Q(g) ≡ u_Q = Q^(−(k−2)/2)·a_Q(f) at one prime 𝔩 above ℓ. Here each match is checked at every prime up to 1,000
# (not dividing 6ℓ), and the first coefficients are shown with their values on lane ℓ. PARI 2.15.4 through cypari2, in WSL.
# Coefficients from PARI; agreements COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_surfaces.py [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(3 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_surfaces.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); QMAX = 1000
CASES = ((17, 14, 64, -729, 1734, 'y^2 - 2'), (19, 16, -128, -2187, 2166, 'y^2 - y - 1'), (23, 20, 512, 19683, 3174, 'y^2 - 3'))
for l, k, a2, a3, N, field in CASES:
    f = None
    for F in pari.mfeigenbasis(pari('mfinit([6, %d], 0)' % k)):
        c = pari.mfcoefs(F, QMAX)
        if int(c[2]) == a2 and int(c[3]) == a3: f = [int(x) for x in c]
    b0 = (k - 2) // 2; Qs = [int(q) for q in pari('primes([5, %d])' % QMAX) if int(q) != l]
    u = {Q: f[Q] * pow(Q, -b0, l) % l for Q in Qs}
    mf = pari('mfinit([%d, 2], 0)' % N); B = pari.mfeigenbasis(mf); flds = [str(P) for P in pari.mffields(mf)]
    for i, F in enumerate(B):
        if flds[i] != field: continue
        c = pari.mfcoefs(F, QMAX); roots = [int(pari.lift(r)) for r in pari('polrootsmod(%s, %d)' % (field, l))]
        for r in roots:
            val = {Q: int(pari('lift(subst(lift(%s), y, Mod(%d, %d)))' % (c[Q], r, l))) for Q in Qs}
            ok = sum(1 for Q in Qs if val[Q] == u[Q])
            if ok < len(Qs): continue                                  # only the full matches (the twin forms match a twist)
            say('ℓ = %d: level %d = 2·3·%d², field %s, y ↦ %d on lane %d: agrees with the tower at %d of %d primes up to %d' % (
                l, N, l, field, r, l, ok, len(Qs), QMAX))
            say('   a_Q of the surface (y as above), and on lane %d / what the tower needs:  %s' % (l, ';  '.join(
                'a_%d = %s → %d / %d' % (Q, pari.lift(c[Q]), val[Q], u[Q]) for Q in Qs[:6])))
say(); say('seconds: %.0f' % (time.time() - t0))
