# A dig, 11 Nov 2026: THE STEPS ALONG Q, BEYOND LANE ℓ. partitions_hida_digits.txt: the constant along ℓ is (3|ℓ)·u(−2), the
# U_ℓ-number of the tower form's ordinary ℓ-adic family at weight −2. If the deep tower IS that family's member at weight −2,
# the steps along every other prime Q follow the same member: c(Q²M) ≡ μ_ε·c(M), μ_ε = (3|Q)·A_Q − ε·(−3|Q)·Q⁻², where A_Q is
# a_Q of the family at weight −2 (Newton's series over the weights ℓ − 3 + (ℓ − 1)·x, as for the constants), ε = (M|Q),
# c(M) = p(N), M = 24N − 1. On lane ℓ this rule was known (the form f itself; for 13 the curve 1014b1); what governs the steps
# beyond lane ℓ was open (README, "Still open after this turn").
# PREDICTED before the run: on exact partition numbers at level b (ℓᵇ dividing 24N − 1), the ratio c(Q²M)/c(M) agrees with
# μ_ε in about b digits, one more per level, as the constants did; f itself agrees in 1 digit, as a rule. Q = 5, 7, 11, both
# classes, the first suitable N at each level. No partition number is used for μ_ε. PARI 2.15.4 (cypari2) in WSL. COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_hida_steps.py [ℓ …] [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(8 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_hida_steps.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); QM = 60; X = 10; PREC = X + 8; K = 10; QS = (5, 7, 11)
BMAX = {13: 7, 17: 6, 19: 6}
kron = lambda a, q: int(pari.kronecker(a, q))
Ls = [int(a) for a in sys.argv[1:] if a.isdigit()] or [13, 17, 19]
PC = {}


def pnum(N):
    if N not in PC: PC[N] = int(pari.numbpart(N))
    return PC[N]


def lane(x, L, e):
    if pari.type(x) == 't_PADIC': x = pari.truncate(x)
    return int(pari.lift(pari.Mod(x, L ** e)))


def agree(a, b, L, cap):                                                   # how many digits a and b share on lane ℓ
    return next((e for e in range(cap) if (a - b) % L ** (e + 1)), cap)


for L in Ls:
    t = time.time(); k0 = L - 3; PR = [int(q) for q in pari('primes([5, %d])' % QM) if int(q) != L]
    s2, s3 = kron(2, L), kron(3, L); f = None; mf = pari('mfinit([6, %d], 0)' % k0)
    for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
        c = pari.mfcoefs(F, max(QM, L))
        if int(pari.poldegree(P)) == 1 and (int(c[2]) > 0) == (s2 > 0) and (int(c[3]) > 0) == (s3 > 0): f = [int(x) for x in c]
    A = {Q: [] for Q in QS}
    for x in range(X + 1):
        k = k0 + (L - 1) * x; m = pari('mfinit([6, %d], 0)' % k); found = []
        for G, PG in zip(pari.mfeigenbasis(m), pari.mffields(m)):
            g = pari.mfcoefs(G, max(QM, L)); deg = int(pari.poldegree(PG))
            try: embs = [None] if deg == 1 else list(pari.polrootspadic(PG, L, PREC + 4))
            except cypari2.PariError: continue
            for r in embs:
                ev = (lambda Q: g[Q]) if r is None else (lambda Q, r=r: pari.subst(pari.lift(g[Q]), 'y', r))
                try:
                    if all(lane(ev(Q), L, 1) == f[Q] % L for Q in PR) and lane(ev(L), L, 1) != 0:
                        found.append({Q: lane(ev(Q), L, PREC) for Q in QS})
                except cypari2.PariError: continue
        if len(found) != 1: say('ℓ = %d, weight %d: %d members of f\'s family found; stopped' % (L, k, len(found))); break
        for Q in QS: A[Q].append(found[0][Q])
    if len(A[QS[0]]) < X + 1: continue
    M = L ** PREC; pred, vals = {}, {}
    for Q in QS:
        d = [sum((-1) ** (j - i) * int(pari.binomial(j, i)) * A[Q][i] for i in range(j + 1)) % M for j in range(X + 1)]
        vals[Q] = [int(pari.valuation(v, L)) if v else PREC for v in d]
        aQ = sum((-1) ** j * d[j] for j in range(X + 1)) % M                    # a_Q of the family at weight −2
        for e in (1, -1):
            dl = e * kron(-3, Q) * pow(Q, -2, L ** K) % L ** K
            pred[Q, e] = ((kron(3, Q) * aQ - dl) % L ** K, (kron(3, Q) * f[Q] - dl) % L ** K)
    say(); say('ℓ = %d: the family at weights %d to %d (%.0f s); differences of a_Q divisible by ℓ to %s' % (
        L, k0, k0 + (L - 1) * X, time.time() - t, '; '.join('Q = %d: %s' % (Q, vals[Q]) for Q in QS)))
    say('   predicted multipliers μ_ε on lane ℓ^%d (weight −2): %s' % (K, ';  '.join(
        'Q = %d, ε = %+d: %d' % (Q, e, pred[Q, e][0]) for Q in QS for e in (1, -1))))
    for b in range(1, BMAX[L] + 1):
        db = pow(24, -1, L ** b); rows = []
        for Q in QS:
            for e in (1, -1):
                j = 0
                while True:
                    Ns = db + L ** b * j; Mx = 24 * Ns - 1; j += 1
                    if Mx % Q == 0 or kron(Mx, Q) != e or pnum(Ns) % L == 0: continue
                    r = pnum(Q * Q * Ns - (Q * Q - 1) // 24) * pow(pnum(Ns), -1, L ** K) % L ** K; break
                rows.append('Q = %d, ε = %+d: %d / %d' % (Q, e, agree(r, pred[Q, e][0], L, K), agree(r, pred[Q, e][1], L, K)))
        say('   level %d: digits shared by the measured ratio with weight −2 / with f itself:  %s  (%.0f s)' % (
            b, ';  '.join(rows), time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
