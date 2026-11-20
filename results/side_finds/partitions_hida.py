# A dig, 11 Nov 2026: THE SECOND DIGIT OF THE CONSTANTS. partitions_split.txt found the first digit of each constant of
# partitions_constants.txt (13 to 31) to be (3|ℓ)·a_ℓ(f), f the tower's level-6 newform of weight ℓ − 3. Suppose the whole
# constant is the U_ℓ-number of f's ordinary ℓ-adic family (Hida) at weight −2, the integral-weight partner of the partitions'
# weight −1/2. Members of the family whose weights agree on lane (ℓ − 1)·ℓ agree on lane ℓ², and the classical weight that
# agrees with −2 there is k = (ℓ − 1)·ℓ − 2 (154 for 13, 270 for 17, 340 for 19).
# PREDICTED before the run: (3|ℓ)·a_ℓ(g) ≡ α_ℓ (mod ℓ²) for the level-6 newform g of weight k that belongs to f's family
# (a_Q(g) ≡ a_Q(f) on lane ℓ for every prime Q ≤ 60 not dividing 6ℓ). α_ℓ on lane ℓ² (partitions_constants.txt, digits lowest
# first): 13 → 6 + 3·13 = 45, 17 → 6 + 12·17 = 210, 19 → 10 + 2·19 = 48. For contrast: f itself on lane ℓ², and the other
# forms of weight k. PARI 2.15.4 (cypari2) in WSL; coefficients exact, ℓ-adic roots from polrootspadic. COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_hida.py [ℓ …] [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(6 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_hida.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); QM = 60
kron = lambda a, q: int(pari.kronecker(a, q))
ALPHA = {13: (6, 3), 17: (6, 12), 19: (10, 2), 23: (5, 20)}                # the first two digits, partitions_constants.txt
Ls = [int(a) for a in sys.argv[1:] if a.isdigit()] or [13, 17, 19]


def on_lane(x, L, e):                                                       # a rational or ℓ-adic integer on lane ℓᵉ
    if pari.type(x) == 't_PADIC': x = pari.truncate(x)
    return int(pari.lift(pari.Mod(x, L ** e)))


for L in Ls:
    t = time.time(); k0, k1 = L - 3, (L - 1) * L - 2; s2, s3 = kron(2, L), kron(3, L); want = ALPHA[L][0] + ALPHA[L][1] * L
    PR = [int(q) for q in pari('primes([5, %d])' % QM) if int(q) != L]
    f = None
    for F, P in zip(pari.mfeigenbasis(pari('mfinit([6, %d], 0)' % k0)), pari.mffields(pari('mfinit([6, %d], 0)' % k0))):
        c = pari.mfcoefs(F, max(QM, L))
        if int(pari.poldegree(P)) == 1 and (int(c[2]) > 0) == (s2 > 0) and (int(c[3]) > 0) == (s3 > 0): f = c
    say(); say('ℓ = %d: the constant on lane %d² is %d;  f (weight %d): (3|ℓ)·a_ℓ(f) on lane ℓ = %d, on lane ℓ² = %d' % (
        L, L, want, k0, kron(3, L) * on_lane(f[L], L, 1) % L, kron(3, L) * on_lane(f[L], L, 2) % (L * L)))
    mf = pari('mfinit([6, %d], 0)' % k1); B = pari.mfeigenbasis(mf); FL = pari.mffields(mf)
    say('   weight %d, level 6: %d newforms in %d Galois orbits, degrees %s  (%.0f s)' % (
        k1, int(pari.mfdim(mf)), len(B), [int(pari.poldegree(P)) for P in FL], time.time() - t))
    rows = []
    for i, (F, P) in enumerate(zip(B, FL)):
        c = pari.mfcoefs(F, max(QM, L)); deg = int(pari.poldegree(P))
        a2, a3 = pari.lift(c[2]), pari.lift(c[3])
        cls = ((int(a2) > 0) == (s2 > 0)) and ((int(a3) > 0) == (s3 > 0))
        if deg == 1: embs = [None]
        else:
            try: embs = list(pari.polrootspadic(P, L, 8))                  # the embeddings into the ℓ-adic integers
            except cypari2.PariError as err: say('   orbit %d: no ℓ-adic roots computed (%s)' % (i + 1, err)); continue
        for r in embs:
            ev = (lambda Q: c[Q]) if r is None else (lambda Q, r=r: pari.subst(pari.lift(c[Q]), 'y', r))
            try: vals = {Q: on_lane(ev(Q), L, 1) for Q in PR}
            except cypari2.PariError: continue
            same = all(vals[Q] == on_lane(f[Q], L, 1) for Q in PR)
            aL = kron(3, L) * on_lane(ev(L), L, 2) % (L * L)
            rows.append((same, cls, i + 1, deg, aL))
    for same, cls, i, deg, aL in rows:
        if same: say('   IN f\'s FAMILY: orbit %d (degree %d)%s: (3|ℓ)·a_ℓ on lane ℓ² = %d  →  %s' % (
            i, deg, '' if cls else ' (other sign class!)', aL, 'EQUAL to the constant' if aL == want else 'not equal (the constant is %d)' % want))
    others = sorted(set(aL for same, cls, i, deg, aL in rows if not same))
    say('   the %d other embeddings into the ℓ-adic integers: (3|ℓ)·a_ℓ on lane ℓ² takes %d values; the constant %d is %samong them' % (
        sum(1 for r in rows if not r[0]), len(others), want, '' if want in others else 'not '))
    say('   (%.0f s)' % (time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
