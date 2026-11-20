# A dig, 12 Nov 2026: THREE FORMS FROM 61 ON, THE FORMS' SIDE. The tower's class of S_{ℓ−3}^new(Γ0(6)) (a₂ of sign (2|ℓ),
# a₃ of sign (3|ℓ)) held one form for 13 to 31 and one conjugate pair for 37 to 47 (partitions_split.txt). Boylan–Webb
# (Trans. AMS 2013, quoted in reports/Partition tower literature check.md) computed the rank of the stable space to be
# ⌊(ℓ + 12)/24⌋ for every prime up to 1297 except 607. PREDICTED before the run: the class has dimension ⌊(ℓ + 12)/24⌋:
# 2 for 53 and 59, 3 for 61 to 83, 4 for 89.
# Then, for the laws (checked on partition numbers by partitions_three_check.py, which reads this file): for ℓ = 53 (two
# forms, the control) and ℓ = 61 (three), the polynomial whose roots are the class's numbers λ = (3|Q)·a_Q on lane ℓ, at
# Q = 5, 7, 11, 13 and at Q = ℓ, written lowest coefficient first and ending in 1 (the product of the characteristic
# polynomials of (3|Q)·a_Q over the class's Galois orbits, whole numbers, reduced mod ℓ). PARI 2.15.4 (cypari2) in WSL.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_three_forms.py [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(8 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_three_forms.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time()
kron = lambda a, q: int(pari.kronecker(a, q))


def tower_class(L):
    s2, s3 = kron(2, L), kron(3, L); mf = pari('mfinit([6, %d], 0)' % (L - 3)); out = []
    for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
        c = pari.mfcoefs(F, 3)
        if (int(pari.lift(c[2])) > 0) == (s2 > 0) and (int(pari.lift(c[3])) > 0) == (s3 > 0): out.append((F, P))
    return out


say('THE CENSUS (the tower\'s class of S_{ℓ−3}^new(Γ0(6)); predicted dimension ⌊(ℓ + 12)/24⌋)')
for L in (53, 59, 61, 67, 71, 73, 79, 83, 89):
    t = time.time(); cl = tower_class(L); degs = [int(pari.poldegree(P)) for F, P in cl]
    say('   ℓ = %d: dimension %d (Galois orbits of degree %s);  predicted %d  →  %s  (%.0f s)' % (
        L, sum(degs), degs, (L + 12) // 24, 'as predicted' if sum(degs) == (L + 12) // 24 else 'NOT as predicted', time.time() - t))
say(); say('THE POLYNOMIALS OF λ = (3|Q)·a_Q ON LANE ℓ (lowest coefficient first, ending in 1)')
for L in (53, 61):
    cl = tower_class(L); QS = (5, 7, 11, 13, L)
    cs = [(pari.mfcoefs(F, L), P) for F, P in cl]
    for Q in QS:
        poly = pari('1')
        for c, P in cs:
            a = kron(3, Q) * (c[Q] if int(pari.poldegree(P)) == 1 else pari.Mod(pari.lift(c[Q]), P))
            poly = poly * (pari('x') - a if int(pari.poldegree(P)) == 1 else pari.charpoly(a, 'x'))
        co = [int(pari.lift(pari.Mod(pari.polcoef(poly, i, 'x'), L))) for i in range(int(pari.poldegree(poly, 'x')) + 1)]
        say('   ℓ = %d, Q = %d: %s' % (L, Q, co))
say(); say('seconds: %.0f' % (time.time() - t0))
