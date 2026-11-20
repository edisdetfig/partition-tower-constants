# A dig, 12 Nov 2026: 5'S AND 11'S STEPS ALONG Q, BEYOND LANE ℓ. The README ("Still open"): for 5 the 1 + Q rule holds on lane
# 5 only, and on lanes 25 and 125 "no shape χ(Q)·Qᵃ·(1 + Q^c) fits"; for 11 the rule fails on lane 121. partitions_hida_laws.txt
# found 5's and 11's constants as the numbers of their slope-2 families at weight −2. If those members are the deep towers,
# their a_Q give the steps on every lane: μ_ε = (3|Q)·A_Q − ε·(−3|Q)·Q⁻², A_Q = a_Q of the family at weight −2 (Newton's
# series over the weights −2 + (ℓ − 1)·ℓ·(x + 1), forms only), as for 7 (partitions_hida_seven_steps.txt).
# PREDICTED before the run: at deeper levels the measured ratio c(Q²M)/c(M) shares more and more digits with μ_ε, while the
# lane-ℓ rule (5: a_Q = Q²(1 + Q); 11: the weight-28 form) stays near 1. Data: exact p(N), the laws' factor ℓᵇ (Watson for
# 5, Atkin for 11) divided out of p(N) and p(N*) at level b (24N − 1 = Q²·(24N* − 1), ℓᵇ dividing 24N* − 1), Q = 7, 11 or 13,
# both classes, the first N* at each level whose p(N*) carries exactly ℓᵇ. PARI 2.15.4 (cypari2) in WSL for the forms;
# partition numbers from sympy (PARI's numbpart was found wrong once, partitions_twins_43_sympy.py). COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2>:<folder with sympy> python3 -B -s .../partitions_hida_law_steps.py [ℓ …] [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(12 * 10 ** 9)
sys.set_int_max_str_digits(0)
from sympy.functions.combinatorial.numbers import partition
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
NAMED = [a for a in sys.argv[1:] if a.isdigit()]
fh = open(os.path.join(HERE, 'partitions_hida_law_steps%s.txt' % ''.join('_' + a for a in NAMED)), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); QM = 60
CONF = {5: dict(X=12, QS=(7, 11, 13), B=10, K=12), 11: dict(X=2, QS=(5, 7, 13), B=7, K=8)}
Ls = [int(a) for a in NAMED] or [5, 11]
kron = lambda a, q: int(pari.kronecker(a, q))
PC = {}


def pnum(N):
    if N not in PC: PC[N] = int(partition(N))
    return PC[N]


def val(n, L):
    v = 0
    while n % L == 0: n //= L; v += 1
    return v


for L in Ls:
    cf = CONF[L]; X, QS, B, K = cf['X'], cf['QS'], cf['B'], cf['K']; PREC = 2 * X + 12
    PR = [int(q) for q in pari('primes([5, %d])' % QM) if int(q) != L]; t = time.time()
    if L == 5: ref = {Q: Q * Q * (1 + Q) for Q in PR}
    else:
        m28 = pari('mfinit([6, 28], 0)')
        for F, P in zip(pari.mfeigenbasis(m28), pari.mffields(m28)):
            c = pari.mfcoefs(F, QM)
            if int(pari.poldegree(P)) == 1 and int(c[2]) == -8192 and int(c[3]) == 1594323: ref = {Q: int(c[Q]) for Q in PR}
    A = {Q: [] for Q in QS}
    for x in range(X + 1):
        k = -2 + (L - 1) * L * (x + 1); m = pari('mfinit([6, %d], 0)' % k); found = []
        for G, PG in zip(pari.mfeigenbasis(m), pari.mffields(m)):
            g = pari.mfcoefs(G, QM); deg = int(pari.poldegree(PG))
            try: embs = [None] if deg == 1 else list(pari.polrootspadic(PG, L, 400))
            except cypari2.PariError: continue
            for r in embs:
                ev = (lambda Q: g[Q]) if r is None else (lambda Q, r=r: pari.subst(pari.lift(g[Q]), 'y', r))
                try:
                    if all(int(pari.lift(pari.Mod(ev(Q), L))) == ref[Q] % L for Q in PR) and int(pari.valuation(ev(L), L)) == 2:
                        found.append({Q: int(pari.lift(pari.Mod(ev(Q), L ** PREC))) for Q in QS})
                except cypari2.PariError: continue
        if len(found) != 1: say('ℓ = %d, weight %d: %d members of slope 2; stopped' % (L, k, len(found))); found = None; break
        for Q in QS: A[Q].append(found[0][Q])
    if found is None: continue
    M = L ** PREC; pred = {}; say(); say('ℓ = %d: the family at weights %d to %d (%.0f s)' % (L, -2 + (L - 1) * L, -2 + (L - 1) * L * (X + 1), time.time() - t))
    for Q in QS:
        d = [sum((-1) ** (j - i) * int(pari.binomial(j, i)) * A[Q][i] for i in range(j + 1)) % M for j in range(X + 1)]
        aQ = sum((-1) ** j * d[j] for j in range(X + 1)) % M
        say('   Q = %d: differences of a_Q divisible by ℓ to %s' % (Q, [int(pari.valuation(v, L)) if v else PREC for v in d]))
        for e in (1, -1):
            dl = e * kron(-3, Q) * pow(Q, -2, L ** K) % L ** K
            pred[Q, e] = ((kron(3, Q) * aQ - dl) % L ** K, (kron(3, Q) * ref[Q] - dl) % L ** K)
    agree = lambda a, b: next((e for e in range(K) if (a - b) % L ** (e + 1)), K)
    for b in range(1, B + 1):
        t = time.time(); db = pow(24, -1, L ** b); rows = []
        for Q in QS:
            for e in (1, -1):
                j = 0
                while True:
                    Ns = db + L ** b * j; Mx = 24 * Ns - 1; j += 1
                    if Mx % Q == 0 or kron(Mx, Q) != e or val(pnum(Ns), L) != b: continue
                    N = Q * Q * Ns - (Q * Q - 1) // 24
                    r = (pnum(N) // L ** b) * pow(pnum(Ns) // L ** b, -1, L ** K) % L ** K; break
                rows.append('Q = %d, ε = %+d: %d / %d' % (Q, e, agree(r, pred[Q, e][0]), agree(r, pred[Q, e][1])))
        say('   level %d (the factor %d^%d out): digits shared with weight −2 / with the lane-ℓ rule:  %s  (%.0f s)' % (
            b, L, b, ';  '.join(rows), time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
