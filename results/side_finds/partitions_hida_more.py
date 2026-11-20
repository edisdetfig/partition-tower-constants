# A dig, 12 Nov 2026: THE CONSTANTS OF 47 AND 61. 47's tower class is one conjugate pair in which 47 stays prime (inert), so
# its two constants live in the field with 47² elements, not on lane 47; 61's class is one cubic orbit (three constants).
# Method: for every prime 𝔭 of the class's coefficient field above ℓ (found by PARI's factorpadic; residue degree f), the
# family members at the weights k = ℓ − 3 + (ℓ − 1)·x are the level-6 newforms with a prime of the same degree whose
# characteristic polynomials of a_Q on lane ℓ agree with 𝔭's for every prime Q ≤ 60 (Q ∤ 6ℓ), with a_ℓ a unit there. For
# each, the unit root u of z² − a_ℓ·z + ℓ^(k−1) is taken in the ℓ-adic completion, and the coefficients of its
# characteristic polynomial (degree f; their symmetric functions are ℓ-adic integers that vary smoothly along the family,
# Hida) are carried to x = −1 (weight −2) by Newton's series, as in partitions_hida_digits.py.
# PREDICTED before the run: on lane ℓ the constants' polynomial equals the one the partition numbers already showed (47:
# x² − 29x + 19, partitions_split.txt; 61: 34 + 15x + 31x² + x³, partitions_three_forms.txt); partitions_three_more_check.py
# then checks it deeper on partition numbers. Weights: 47 up to 274 (6 weights), 61 up to 298 (5), below PARI's crash zone.
# PARI 2.15.4 (cypari2) in WSL. No partition number is used. COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_hida_more.py [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(12 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_hida_more.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); QM = 60; FPREC = 150
CASES = {47: 5, 61: 4}
kron = lambda a, q: int(pari.kronecker(a, q))


def tower_class(L):
    s2, s3 = kron(2, L), kron(3, L); mf = pari('mfinit([6, %d], 0)' % (L - 3)); out = []
    for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
        c = pari.mfcoefs(F, max(QM, L))
        if (int(pari.lift(c[2])) > 0) == (s2 > 0) and (int(pari.lift(c[3])) > 0) == (s3 > 0): out.append((c, P))
    return out


def local(c, P, L, PR, k, prec):
    """For an orbit (coefficients c, field P): per prime above L: (f, signature on lane L, charpoly of the unit root or None)."""
    out = []
    if int(pari.poldegree(P)) == 1:
        hs = [None]
    else:
        fa = pari.factorpadic(P, L, FPREC)
        if any(int(e) != 1 for e in fa[1]): return out                       # not squarefree over Q_ℓ: skipped
        hs = list(fa[0])
    for h in hs:
        img = (lambda Q: c[Q] + pari('O(%d^%d)' % (L, FPREC))) if h is None else (lambda Q, h=h: pari.Mod(pari.lift(c[Q]), h))
        cpz = lambda z: [z] if h is None else None
        try:
            sig = []
            for Q in PR:
                if h is None: sig.append((int(pari.lift(pari.Mod(c[Q], L))),))
                else: sig.append(tuple(int(pari.lift(pari.Mod(pari.polcoef(pari.charpoly(img(Q), 'x'), i, 'x'), L))) for i in range(int(pari.poldegree(h)))))
            sig = tuple(sig); f = 1 if h is None else int(pari.poldegree(h))
            a = img(L); na = a if h is None else pari.norm(a)
            if int(pari.valuation(na, L)) != 0: out.append((f, sig, None)); continue
            z = a
            for _ in range(prec + 4): z = a - pari('%d^%d' % (L, k - 1)) / z                      # the unit root
            if h is None: cp = [-z, pari('1')]
            else:
                cpol = pari.charpoly(z, 'x'); cp = [pari.polcoef(cpol, i, 'x') for i in range(f + 1)]
            out.append((f, sig, [int(pari.lift(pari.Mod(x, L ** prec))) for x in cp]))
        except cypari2.PariError: continue
    return out


for L, X in CASES.items():
    t = time.time(); PREC = X + 8; PR = [int(q) for q in pari('primes([5, %d])' % QM) if int(q) != L]
    ref = []
    for c, P in tower_class(L): ref += [(f, sig) for f, sig, cp in local(c, P, L, PR, L - 3, PREC)]
    say(); say('ℓ = %d: the class\'s primes above ℓ have residue degrees %s' % (L, [f for f, s in ref]))
    series = [[] for _ in ref]; ok = True
    for x in range(X + 1):
        k = L - 3 + (L - 1) * x; mf = pari('mfinit([6, %d], 0)' % k); found = [[] for _ in ref]
        for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
            c = pari.mfcoefs(F, max(QM, L))
            for f, sig, cp in local(c, P, L, PR, k, PREC):
                for i, (fr, sr) in enumerate(ref):
                    if f == fr and sig == sr and cp is not None: found[i].append(cp)
        if any(len(fd) != 1 for fd in found):
            say('   weight %d: members found per prime %s; stopped' % (k, [len(fd) for fd in found])); ok = False; break
        for i in range(len(ref)): series[i].append(found[i][0])
    if not ok: continue
    say('   weights %s  (%.0f s)' % (', '.join(str(L - 3 + (L - 1) * x) for x in range(X + 1)), time.time() - t))
    M = L ** PREC; s3 = kron(3, L); total = [1]
    for i, (f, sig) in enumerate(ref):
        coefs = []; vals = []
        for j in range(f + 1):
            u = [series[i][x][j] for x in range(X + 1)]
            d = [sum((-1) ** (n - m) * int(pari.binomial(n, m)) * u[m] for m in range(n + 1)) % M for n in range(X + 1)]
            vals.append([int(pari.valuation(v, L)) if v else PREC for v in d])
            coefs.append(sum((-1) ** n * d[n] for n in range(X + 1)) % M)
        al = [coefs[j] * s3 ** (f - j) % M for j in range(f + 1)]                 # the polynomial of α = (3|ℓ)·u(−2)
        say('   prime %d (degree %d): differences of the coefficients divisible by ℓ to %s' % (i + 1, f, vals[:f]))
        if f == 1:
            a1 = (-al[0]) % M
            say('      its constant, digits lowest first: %s' % [(a1 // L ** e) % L for e in range(X + 2)])
        elif f == 2:
            s1, s2 = (-al[1]) % M, al[0] % M                                      # α + α′ and α·α′
            r = next(r for r in range(2, L) if kron(r, L) == -1); disc = (s1 * s1 - 4 * s2) % M
            cc = int(pari.lift(pari.Mod(pari.sqrt(pari('(%d + O(%d^%d)) / %d' % (disc, L, PREC, r))), M)))
            inv2 = pow(2, -1, M); A, B = s1 * inv2 % M, cc * inv2 % M
            say('      its two constants are A ± B·√%d (%d has no square root on lane %d):' % (r, r, L))
            say('         A, digits lowest first: %s' % [(A // L ** e) % L for e in range(X + 2)])
            say('         B, digits lowest first: %s' % [(B // L ** e) % L for e in range(X + 2)])
            say('         their sum and product, digits: %s;  %s' % ([(s1 // L ** e) % L for e in range(X + 2)], [(s2 // L ** e) % L for e in range(X + 2)]))
        else:
            say('      the polynomial of its constants (lowest coefficient first): %s' % al)
        new = [0] * (len(total) + f)
        for a in range(len(total)):
            for b in range(f + 1): new[a + b] = (new[a + b] + total[a] * al[b]) % M
        total = new
    say('   the constants\' polynomial on lane ℓ (lowest first): %s;  on lane ℓ^%d: %s' % ([v % L for v in total], X + 1, [v % L ** (X + 1) for v in total]))
say(); say('seconds: %.0f' % (time.time() - t0))
