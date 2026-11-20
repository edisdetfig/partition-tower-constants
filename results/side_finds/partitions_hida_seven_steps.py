# A dig, 11 Nov 2026: 7'S STEPS ALONG Q, BEYOND LANE 49. partitions_curve49/343 found 7's multipliers on lane 49 to be E7's,
# λ ≡ (3|Q)·Q¹⁹·a_Q(E7) (E7 = 42a1), and that this fails on lane 343. partitions_hida_seven.txt found 7's constant (13 of 14
# digits) as the number of the slope-1 family of 7's tower system at weight −2 (weights 40 + 42·x). If that member is the
# deep tower, its a_Q give the steps along Q on every lane: μ_ε = (3|Q)·A_Q − ε·(−3|Q)·Q⁻², A_Q = a_Q of the family at
# weight −2 (Newton's series over the weights 40 to 292, forms only).
# PREDICTED before the run: at deeper levels the measured ratio c(Q²M)/c(M) shares more and more digits with μ_ε, past
# lane 49, while E7's rule stays at 2 digits. Data: exact p(N) (numbpart), Watson's factor 7^(⌊b/2⌋+1) divided out of both
# p(N) and p(N*) at level b (24N − 1 = Q²·(24N* − 1), 7ᵇ dividing 24N* − 1), Q = 5, 11, 13, both classes, the first N* at
# each level whose p(N*) carries exactly Watson's factor. PARI 2.15.4 (cypari2) in WSL. COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_hida_seven_steps.py [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(8 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_hida_seven_steps.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); L = 7; QM = 60; X = 6; PREC = 16; K = 12; QS = (5, 11, 13); BMAX = 9   # (a trial capped at 8 digits)
PR = [int(q) for q in pari('primes([5, %d])' % QM) if int(q) != L]
kron = lambda a, q: int(pari.kronecker(a, q))
PC = {}


def pnum(N):
    if N not in PC: PC[N] = int(pari.numbpart(N))
    return PC[N]


def lane(x, e):
    if pari.type(x) == 't_PADIC': x = pari.truncate(x)
    return int(pari.lift(pari.Mod(x, L ** e)))


def agree(a, b, cap): return next((e for e in range(cap) if (a - b) % L ** (e + 1)), cap)


def v7(n):
    v = 0
    while n % 7 == 0: n //= 7; v += 1
    return v


f = None
mf = pari('mfinit([6, 16], 0)')
for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
    c = pari.mfcoefs(F, QM)
    if int(pari.poldegree(P)) == 1 and int(c[2]) == 128 and int(c[3]) == -2187: f = [int(x) for x in c]
E = [int(x) for x in pari.mfcoefs(pari.mfeigenbasis(pari('mfinit([42, 2], 0)'))[0], QM)]   # 42a1, the only class of level 42
A = {Q: [] for Q in QS}
for x in range(X + 1):
    k = 40 + 42 * x; m = pari('mfinit([6, %d], 0)' % k); found = []
    for G, PG in zip(pari.mfeigenbasis(m), pari.mffields(m)):
        g = pari.mfcoefs(G, QM); deg = int(pari.poldegree(PG))
        try: embs = [None] if deg == 1 else list(pari.polrootspadic(PG, L, 400))
        except cypari2.PariError: continue
        for r in embs:
            ev = (lambda Q: g[Q]) if r is None else (lambda Q, r=r: pari.subst(pari.lift(g[Q]), 'y', r))
            try:
                if all(lane(ev(Q), 1) == f[Q] % L for Q in PR) and int(pari.valuation(ev(7), L)) == 1:
                    found.append({Q: lane(ev(Q), PREC) for Q in QS})
            except cypari2.PariError: continue
    if len(found) != 1: say('weight %d: %d members of slope 1; stopped' % (k, len(found))); sys.exit()
    for Q in QS: A[Q].append(found[0][Q])
M = L ** PREC; pred = {}
for Q in QS:
    d = [sum((-1) ** (j - i) * int(pari.binomial(j, i)) * A[Q][i] for i in range(j + 1)) % M for j in range(X + 1)]
    aQ = sum((-1) ** j * d[j] for j in range(X + 1)) % M
    say('Q = %d: differences of a_Q divisible by 7 to %s' % (Q, [int(pari.valuation(v, L)) if v else PREC for v in d]))
    for e in (1, -1):
        dl = e * kron(-3, Q) * pow(Q, -2, L ** K) % L ** K
        pred[Q, e] = ((kron(3, Q) * aQ - dl) % L ** K, (kron(3, Q) * pow(Q, 19, L ** K) * E[Q] - dl) % L ** K)
for b in range(1, BMAX + 1):
    t = time.time(); db = pow(24, -1, L ** b); w = b // 2 + 1; rows = []
    for Q in QS:
        for e in (1, -1):
            j = 0
            while True:
                Ns = db + L ** b * j; Mx = 24 * Ns - 1; j += 1
                if Mx % Q == 0 or kron(Mx, Q) != e or v7(pnum(Ns)) != w: continue
                N = Q * Q * Ns - (Q * Q - 1) // 24
                r = (pnum(N) // 7 ** w) * pow(pnum(Ns) // 7 ** w, -1, L ** K) % L ** K; break
            rows.append('Q = %d, ε = %+d: %d / %d' % (Q, e, agree(r, pred[Q, e][0], K), agree(r, pred[Q, e][1], K)))
    say('level %d (Watson\'s factor 7^%d out): digits shared with weight −2 / with E7\'s rule:  %s  (%.0f s)' % (b, w, ';  '.join(rows), time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
