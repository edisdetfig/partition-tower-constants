# A dig, 11 Nov 2026: TESTING THE m = (ℓ ± 1)/2 RULE AT 37, 41, 43, 47. The ninth dig found, for every prime ℓ from 7 to 31,
# that the weight-2 newforms carrying the tower's rule on lane ℓ have coefficient field Q(ζ_m)⁺ with m = (ℓ ± 1)/2, so the
# smallest has dimension min φ((ℓ ± 1)/2)/2; set down BEFORE this run: 37 → 3 (m = 18), 41 → 4 (m = 20), 43 → 5 (m = 22),
# 47 → 4 (m = 24). From 37 on the tower's rule need not be one form (Ahlgren–Allen–Tang: for 37 the relevant space is two-
# dimensional), so every level-6 newform of weight ℓ − 3 in the tower's sign class is tested; the class is read off the
# primes 13 to 31 (partitions_newforms.txt: a₂ > 0 exactly when (2|ℓ) = +1, a₃ > 0 exactly when (3|ℓ) = +1).
# For each such form f (and each prime λ above ℓ of its field with residue degree 1), the target u_Q = Q^(−(k−2)/2)·a_Q(f)
# mod λ is compared, as in partitions_abelian.py, with every weight-2 newform of trivial character at the levels ℓ, 2ℓ, 3ℓ,
# 6ℓ, ℓ², 2ℓ², 3ℓ², 6ℓ² whose Galois orbit has dimension at most D (PARI's mfsplit with a dimension limit, so the large orbits
# are never split), twists by ±1, ±2, ±3, ±6 and those times ℓ, at every prime Q ≤ 300 not dividing 6ℓ. PARI 2.15.4 (cypari2).
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_abelian_more.py ℓ D [--all] [--keep]
#   (--all: every level-6 newform of weight ℓ − 3 is a target; the .txt is then named partitions_abelian_more_<ℓ>_all_<D>.txt)
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(6 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
args = [int(a) for a in sys.argv[1:] if a.isdigit()]; L, DLIM = args[0], args[1]
ALL = '--all' in sys.argv                     # every level-6 newform of weight ℓ − 3 as a target, not only the tower's sign class
OUT = os.path.join(HERE, 'partitions_abelian_more_%d%s.txt' % (L, '_all_%d' % DLIM if ALL else ''))
if KEEP and os.path.exists(OUT): OUT = OUT[:-4] + '_%d.txt' % DLIM      # a later run with a larger limit (added for 47 at 8)
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); QMAX = 300; k = L - 3; b0 = (k - 2) // 2
PR = [int(q) for q in pari('primes([5, %d])' % QMAX) if int(q) != L]
s2, s3 = int(pari.kronecker(2, L)), int(pari.kronecker(3, L))
real_cyclo = {}
for m in ((L - 1) // 2, (L + 1) // 2):
    real_cyclo[m] = pari('polredabs(minpoly(Mod(y + y^-1, polcyclo(%d, y)), y))' % m) if m > 2 else pari('y')
say('ℓ = %d: weight %d; the tower\'s sign class: a₂ %s 0, a₃ %s 0;  predicted smallest dimension %d (m = %d: %d, m = %d: %d);  orbits up to dimension %d' % (
    L, k, '>' if s2 > 0 else '<', '>' if s3 > 0 else '<', min(int(pari.eulerphi(m)) // 2 if m > 2 else 1 for m in real_cyclo),
    (L - 1) // 2, int(pari.eulerphi((L - 1) // 2)) // 2, (L + 1) // 2, int(pari.eulerphi((L + 1) // 2)) // 2, DLIM))

# ---- the level-6 forms of the tower's sign class, and their targets ----
mf6 = pari('mfinit([6, %d], 0)' % k); B6 = pari.mfeigenbasis(mf6); F6 = pari.mffields(mf6); targets = []
for i, F in enumerate(B6):
    c = pari.mfcoefs(F, QMAX); a2, a3 = int(pari.lift(c[2])), int(pari.lift(c[3]))
    if not ALL and ((a2 > 0) != (s2 > 0) or (a3 > 0) != (s3 > 0)): continue
    P = F6[i]; deg = int(pari.poldegree(P))
    if deg == 1: hs = [pari('t')]                                          # the field with ℓ elements (F_ℓ[t]/(t))
    else: hs = [pari.subst(pari.lift(h), 'y', pari('t')) for h in pari.factormod(P, L)[0]]
    say('   level-6 newform %d: field of degree %d (%s), a₂ = %+d, a₃ = %+d; its primes above %d have residue degrees %s' % (
        i + 1, deg, P, a2, a3, L, [int(pari.poldegree(h)) for h in hs] if deg > 1 else [1]))
    for h in hs:                                                           # one target per prime λ above ℓ, in the field F_λ
        T = pari('Mod(Mod(1, %d) * t, %s)' % (L, h))
        try:
            u = {Q: (pari.subst(pari.lift(c[Q]), 'y', T) + 0 * T) * pari('Mod(%d, %d)' % (pow(Q, -b0, L), L)) for Q in PR}
            targets.append(('f%d%s' % (i + 1, '' if deg == 1 else ' at the prime %s' % h), h, T, u))
        except cypari2.PariError: say('      (a denominator divisible by %d at this prime; skipped)' % L)
say('   %d targets' % len(targets))
# (the residue-field handling above, needed for 47 where ℓ is inert in the pair's field, was added after the runs for 37, 41,
#  43; for those the pair's primes above ℓ have degree 1 and the comparison is the same as before.)

# ---- the weight-2 search ----
found = []
for N in (L, 2 * L, 3 * L, 6 * L, L * L, 2 * L * L, 3 * L * L, 6 * L * L):
    t = time.time(); mf = pari('mfinit([%d, 2], 0)' % N); dim = int(pari.mfdim(mf))
    if dim == 0: continue
    vF, vK = pari.mfsplit(mf, DLIM); C = pari.mfcoefs(mf, QMAX); degs = [int(pari.poldegree(P)) for P in vK]
    for i in range(len(vF)):
        P = vK[i]; v = C * vF[i]; v = v * (v[1] ** -1) if degs[i] > 1 else v / v[1]       # normalized: a₁ = 1
        for name, h, T, u in targets:
            if degs[i] == 1: roots = [None]
            else:
                try: roots = list(pari('polrootsmod(%s, [%s, %d])' % (P, h, L)))     # embeddings of the field into F_λ
                except cypari2.PariError: continue
            for rho in roots:
                try:
                    val = {Q: (pari.subst(pari.lift(v[Q]), 'y', rho) if rho is not None else v[Q]) + 0 * T for Q in PR}
                except cypari2.PariError: continue
                for d in (1, -1, 2, -2, 3, -3, 6, -6, L, -L, 2 * L, -2 * L, 3 * L, -3 * L, 6 * L, -6 * L):
                    if all(val[Q] == int(pari.kronecker(d, Q)) * u[Q] for Q in PR):
                        cyc = [m for m, Pm in real_cyclo.items() if degs[i] == int(pari.poldegree(Pm)) and pari.polredabs(P) == Pm]
                        found.append((N, degs[i], str(P), name, d, cyc))
    say('   level %5d: %4d newforms; orbits of dimension ≤ %d: %s  (%.0f s)' % (N, dim, DLIM, ' '.join(map(str, sorted(degs))) or 'none', time.time() - t))
for N, dg, P, name, d, cyc in sorted(found, key=lambda x: (x[1], x[0])):
    say('   MATCH: level %d, dimension %d, field %s%s, target %s, twist d = %d' % (
        N, dg, P, ' = Q(ζ_%d)⁺' % cyc[0] if cyc else '', name, d))
if not found: say('   no match among the orbits of dimension ≤ %d' % DLIM)
say('   smallest dimension found: %s' % (min(f[1] for f in found) if found else 'none ≤ %d' % DLIM))
say(); say('seconds: %.0f' % (time.time() - t0))
