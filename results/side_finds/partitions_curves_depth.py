# A dig, 12 Nov 2026: HOW DEEP THE ELLIPTIC CURVES GO. On lane ℓ each small tower's steps were found to follow a weight-2
# object, twisted by a power of Q: for 5 the Eisenstein shape 1 + Q (with Q²), for 7 Jameson's E7 = 42a1 (Q¹; on lane 49 Q¹⁹),
# for 11 her E11 = 726d1 (Q³), for 13 the curve 1014b1 (Q⁴; partitions_curve13*). The towers themselves are the weight −2
# members of ℓ-adic families (partitions_hida_steps, _seven_steps, _law_steps). A weight-2 object twisted by Q^j acts like
# weight 2 + 2j, so to agree with weight −2 on lane ℓᵉ it needs 2 + 2j ≡ −2 modulo (ℓ − 1)·ℓ^(e−1): on lane ℓ that gives the
# exponents 1, 3, 4 for 7, 11, 13 (and 2 for 5); on lane ℓ² it gives 19 (mod 21), 53 (mod 55), 76 (mod 78), 8 (mod 10).
# Here each family member's a_Q at weight −2 (Newton's series over the family's weights, forms only) is compared with
# s·(d|Q)·Q^j·a_Q(E) for every sign s, quadratic twist d ∈ {±1, ±2, ±3, ±6}, and exponent j modulo (ℓ − 1)·ℓ², at every
# prime Q from 5 to 100 not dividing 6ℓ; the score of a choice is the number of ℓ-adic digits it gets right at every Q.
# PREDICTED before the run: E7 scores 2 (partitions_curve343 found no exponent reaching lane 343); the Eisenstein shape for
# 5 and E11 score 1 (the README: no shape on lane 25; E11 fails on lane 121); 1014b1 for 13: not tested before — if the
# weight rule is the whole story it scores 2 with j ≡ 76 (mod 78). (A trial before this kept run: 13's curve scored 1.)
# PARI 2.15.4 (cypari2) in WSL. COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_curves_depth.py [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(12 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_curves_depth.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); QM = 60; QT = 100; K = 3
kron = lambda a, q: int(pari.kronecker(a, q))


def rational_form(k, a2, a3):
    m = pari('mfinit([6, %d], 0)' % k)
    for F, P in zip(pari.mfeigenbasis(m), pari.mffields(m)):
        c = pari.mfcoefs(F, QM)
        if int(pari.poldegree(P)) == 1 and int(c[2]) == a2 and int(c[3]) == a3: return [int(x) for x in c]


# ℓ: (reference system on lane ℓ, slope, first weight, weight step, number of weights, the weight-2 object)
CASES = {13: ('form', (10, -16, 81), 0, 10, 12, 11, [1, 0, 1, -10482, 722308], '1014b1'),
         7: ('form', (16, 128, -2187), 1, 40, 42, 7, [1, 1, 1, -4, 5], '42a1 (E7)'),
         11: ('form', (28, -8192, 1594323), 2, 108, 110, 3, [1, 0, 1, -14, 20], '726d1 (E11)'),
         5: ('eis', None, 2, 18, 20, 13, None, 'the Eisenstein shape 1 + Q')}
for L in (5, 7, 11, 13):
    kind, fa, slope, k0, step, nw, ainv, name = CASES[L]; t = time.time(); PREC = 2 * nw + 10
    PR = [int(q) for q in pari('primes([5, %d])' % QM) if int(q) != L]
    QS = [int(q) for q in pari('primes([5, %d])' % QT) if int(q) != L]
    if kind == 'eis': ref = {Q: Q * Q * (1 + Q) % L for Q in PR}
    else: c = rational_form(*fa); ref = {Q: c[Q] % L for Q in PR}
    A = {Q: [] for Q in QS}; ok = True
    for x in range(nw):
        k = k0 + step * x; m = pari('mfinit([6, %d], 0)' % k); found = []
        for G, PG in zip(pari.mfeigenbasis(m), pari.mffields(m)):
            g = pari.mfcoefs(G, QT); deg = int(pari.poldegree(PG))
            try: embs = [None] if deg == 1 else list(pari.polrootspadic(PG, L, 400))
            except cypari2.PariError: continue
            for r in embs:
                ev = (lambda Q: g[Q]) if r is None else (lambda Q, r=r: pari.subst(pari.lift(g[Q]), 'y', r))
                try:
                    if all(int(pari.lift(pari.Mod(ev(Q), L))) == ref[Q] for Q in PR) and int(pari.valuation(ev(L), L)) == slope:
                        found.append({Q: int(pari.lift(pari.Mod(ev(Q), L ** PREC))) for Q in QS})
                except cypari2.PariError: continue
        if len(found) != 1: say('ℓ = %d, weight %d: %d members; stopped' % (L, k, len(found))); ok = False; break
        for Q in QS: A[Q].append(found[0][Q])
    if not ok: continue
    M = L ** PREC; AQ = {}
    for Q in QS:
        d = [sum((-1) ** (j - i) * int(pari.binomial(j, i)) * A[Q][i] for i in range(j + 1)) % M for j in range(nw)]
        AQ[Q] = sum((-1) ** j * d[j] for j in range(nw)) % M                    # a_Q of the family at weight −2
    if kind == 'eis': aE = {Q: 1 + Q for Q in QS}
    else: E = pari.ellinit(ainv); aE = {Q: int(pari.ellap(E, Q)) for Q in QS}
    MK = L ** K; per = (L - 1) * L ** (K - 1)
    def score(s, dd, j):
        best = K
        for Q in QS:
            v = s * kron(dd, Q) * pow(Q, j, MK) * aE[Q] - AQ[Q]
            e = next((e for e in range(K) if v % L ** (e + 1)), K); best = min(best, e)
            if best == 0: break
        return best
    res = sorted(((score(s, dd, j), s, dd, j) for s in (1, -1) for dd in (1, -1, 2, -2, 3, -3, 6, -6) for j in range(per)),
                 key=lambda z: -z[0])
    top = res[0][0]; wins = [(s, dd, j) for sc, s, dd, j in res if sc == top]
    pj = [j for j in range(per) if (2 + 2 * j - (-2)) % ((L - 1) * L) == 0]                # the weight rule on lane ℓ²
    say(); say('ℓ = %d, %s against the weight −2 member (%d weights from %d, step %d; %.0f s):' % (L, name, nw, k0, step, time.time() - t))
    say('   the best score: %d digit%s right at all %d primes Q; reached by %d choices, e.g. %s' % (
        top, '' if top == 1 else 's', len(QS), len(wins), ', '.join('sign %+d, d = %d, Q^%d' % w for w in wins[:4])))
    def two(s, dd, j):                                                         # primes Q right on lane ℓ² too
        return sum(1 for Q in QS if (s * kron(dd, Q) * pow(Q, j, MK) * aE[Q] - AQ[Q]) % (L * L) == 0)
    say('   the weight rule on lane ℓ² asks for j ≡ %s (mod %d); best score, and the most primes Q right on lane ℓ²: %s' % (
        ', '.join(map(str, pj[:4])), per, ';  '.join('Q^%d: %d, %d of %d' % (j, max(score(s, dd, j) for s in (1, -1) for dd in (1, -1, 2, -2, 3, -3, 6, -6)),
        max(two(s, dd, j) for s in (1, -1) for dd in (1, -1, 2, -2, 3, -3, 6, -6)), len(QS)) for j in pj[:4])))
say(); say('seconds: %.0f' % (time.time() - t0))
