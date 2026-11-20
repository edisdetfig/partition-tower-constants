# A dig, 11 Nov 2026: WHAT THE SPLIT OF THE TOWER'S RULE TELLS (37, 41, 43, 47). From 37 on the tower (ℓ dividing 24N − 1) is
# no longer one Hecke eigenform on lane ℓ (partitions_tower37.txt). Its class holds one conjugate pair f, f′ of level-6 newforms
# of weight ℓ − 3 over a real quadratic field, so each prime Q ≠ ℓ gives two numbers (3|Q)·a_Q(f), (3|Q)·a_Q(f′). Their sum t and
# product n are whole numbers. With c(M) = p(N) on lane ℓ, M = 24N − 1, ε = (M|Q), δ = −(−3|Q)·Q⁻², exact values (numbpart):
#   1  ALONG Q: c(Q⁴M) ≡ (t + εδ)·c(Q²M) − (Q⁻³ + δ² + tεδ + n)·c(M), DERIVED from the step rule as T² − tT + n = 0.
#   2  LAWS NEED BOTH NUMBERS: t ≡ −2εδ, n ≡ δ² gives p(N) ≡ 0 on a class (modulo ℓQ³); t ≡ n ≡ 0 gives the classes modulo ℓQ⁴.
#      Ahlgren–Allen–Tang list 6599, 7541, 9547 for 37 (Q < 10,000). The smallest members are computed.
#   3  ALONG ℓ: H_b(n) = p(ℓᵇn + δ_b) = c(ℓᵇ·m). PREDICTED before this run (the Shimura correspondence, with the twist (3|Q)
#      found for every Q ≠ ℓ): H₅ ≡ T·H₃ − N·H₁ with T ≡ (3|ℓ)·Tr a_ℓ(f), N ≡ Norm a_ℓ(f). For one form (13 to 31, the
#      control): H₃ ≡ (3|ℓ)·a_ℓ(f)·H₁, so the first digit of each constant of partitions_constants.txt would be a Hecke number.
# PARI 2.15.4 (cypari2) in WSL. Counts COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_split.py [--quick] [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(4 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv; QUICK = '--quick' in sys.argv
fh = open(os.path.join(HERE, 'partitions_split.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); QBIG = 1000 if QUICK else 10000
kron = lambda a, q: int(pari.kronecker(a, q))
md = lambda x, l: int(pari.lift(pari.Mod(x, l)))                 # a rational number (denominator prime to l) on lane l
PC = {}


def pm(N, l):
    if N not in PC: PC[N] = pari.numbpart(N)
    return int(PC[N] % l)


def sign_class(L):                                               # the level-6 newforms of weight ℓ − 3 in the tower's class
    s2, s3 = kron(2, L), kron(3, L); mf = pari('mfinit([6, %d], 0)' % (L - 3)); out = []
    for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
        c = pari.mfcoefs(F, 3)
        if (int(pari.lift(c[2])) > 0) == (s2 > 0) and (int(pari.lift(c[3])) > 0) == (s3 > 0): out.append((F, P))
    return out


say('THE SPLIT TOWERS (the pair f, f′; t and n are the trace and norm of (3|Q)·a_Q(f), on lane ℓ)')
for L in (37, 41, 43, 47):
    t = time.time(); cls = sign_class(L); assert len(cls) == 1
    F, P = cls[0]; assert int(pari.poldegree(P)) == 2
    roots = [int(pari.lift(r)) for r in pari('polrootsmod(%s, %d)' % (P, L))]; split = len(roots) == 2
    cs = pari.mfcoefs(F, QBIG)
    tn = lambda Q: (kron(3, Q) * md(pari.trace(cs[Q]), L) % L, md(pari.norm(cs[Q]), L))
    lam = lambda Q: [kron(3, Q) * int(pari('lift(subst(%s, y, Mod(%d, %d)))' % (pari.lift(cs[Q]), r, L))) % L for r in roots]
    say(); say('ℓ = %d: field %s, discriminant %s; ℓ %s' % (L, P, pari.factor(pari.nfdisc(P)),
        'splits: two rules on lane %d, from y ↦ %s' % (L, roots) if split else 'stays whole (inert): the two rules live in the field with %d² elements' % L))
    say('   at small Q:  %s' % ';  '.join('Q = %d: %s t = %d, n = %d' % (Q, 'λ = %s,' % lam(Q) if split else '', *tn(Q)) for Q in (5, 7, 11, 13)))
    if not split:
        out = [Q for Q in [int(q) for q in pari('primes([5, 300])')] if Q != L and md(pari.polcoef(pari.lift(cs[Q]), 1), L) != 0]
        say('   the numbers leave the field with %d elements at %d of the %d primes 5 ≤ Q ≤ 300 (first: %s)' % (
            L, len(out), int(pari.primepi(300)) - 3, out[:5]))

    # ---- 1. along Q ----
    d1 = pow(24, -1, L); ok2 = tot2 = ok1 = 0; rows = []
    for Q in (5, 7, 11, 13):
        dl = -kron(-3, Q) * pow(Q, -2, L) % L; tq, nq = tn(Q); Q3 = pow(Q, -3, L); got = j = okq = 0
        while got < 12:
            Ns = L * j + d1; M = 24 * Ns - 1; j += 1
            if M % Q == 0: continue
            got += 1; e = kron(M, Q)
            c0, c1, c2 = pm(Ns, L), pm(Q * Q * Ns - (Q * Q - 1) // 24, L), pm(Q ** 4 * Ns - (Q ** 4 - 1) // 24, L)
            good = (c2 - (tq + e * dl) * c1 + (Q3 + dl * dl + tq * e * dl + nq) * c0) % L == 0; ok2 += good; okq += good; tot2 += 1
            if split: ok1 += any((c1 - (lm + e * dl) * c0) % L == 0 for lm in lam(Q))
        rows.append('Q = %d: %d of 12' % (Q, okq))
    say('   1  along Q, the two-step law: %d of %d chains right (%s)%s' % (ok2, tot2, ', '.join(rows),
        ';  one rule alone (either λ): %d of %d' % (ok1, tot2) if split else ''))

    # ---- 2. laws need both numbers ----
    acc, acc0 = [], []
    for Q in [int(q) for q in pari('primes([5, %d])' % (QBIG - 1))]:
        if Q == L: continue
        tq, nq = tn(Q); dl = -kron(-3, Q) * pow(Q, -2, L) % L
        for e in (1, -1):
            if tq == (-2 * e * dl) % L and nq == dl * dl % L: acc.append((Q, kron(L, Q) * e))   # class written as (n|Q), M = ℓn
        if tq == 0 and nq == 0: acc0.append(Q)
    say('   2  primes Q < %d with both numbers in place for a class (24N − 1 = Q²·ℓ·n, class (n|Q)): %s' % (
        QBIG, ', '.join('%d ((n|Q) = %+d)' % qe for qe in acc) or 'none'))
    say('      with both numbers 0 (24N − 1 = ℓ·Q³·m, Q ∤ m): %s   (%.0f s so far)' % (', '.join(map(str, acc0)) or 'none', time.time() - t))
    if L == 37: say('      Ahlgren–Allen–Tang (QUOTED): "we have a congruence (1.6) when Q = 6599, 7541, and 9547"; their examples: 6599 with (n|6599) = −1, 7541 with (n|7541) = 1')
    for Q, e in acc[:3]:
        Q2 = Q * Q; mem, oth = [], []; j = 0
        while len(mem) < 2 or len(oth) < 1:
            Ns = L * j + d1; M = 24 * Ns - 1; j += 1
            if M % Q == 0: continue
            ce = kron(L, Q) * kron(M, Q); N = Q2 * Ns - (Q2 - 1) // 24
            if ce == e and len(mem) < 2: mem.append((N, pm(N, L)))
            elif ce == -e and len(oth) < 1: oth.append((N, pm(N, L)))
        say('      Q = %d: members %s;  the other class: %s  (p(N) on lane %d, exact)' % (
            Q, ', '.join('p(%d) ≡ %d' % x for x in mem), ', '.join('p(%d) ≡ %d' % x for x in oth), L))
    for Q in [q for q in acc0 if q < 1000][:2]:
        m = (-pow(L * Q ** 3, -1, 24)) % 24; mem = []
        while len(mem) < 2:
            if m % Q: N = (L * Q ** 3 * m + 1) // 24; mem.append((N, pm(N, L)))
            m += 24
        say('      Q = %d: members %s  (exact)' % (Q, ', '.join('p(%d) ≡ %d' % x for x in mem)))

    # ---- 3. along ℓ ----
    H = {b: [pm(L ** b * n + pow(24, -1, L ** b), L) for n in range(10)] for b in (1, 3, 5)}
    sols = [(T, Nn) for T in range(L) for Nn in range(L) if all((H[5][n] - T * H[3][n] + Nn * H[1][n]) % L == 0 for n in range(10))]
    A13 = [A for A in range(L) if all((H[3][n] - A * H[1][n]) % L == 0 for n in range(10))]
    A35 = [A for A in range(L) if all((H[5][n] - A * H[3][n]) % L == 0 for n in range(10))]
    Tp, Np = kron(3, L) * md(pari.trace(cs[L]), L) % L, md(pari.norm(cs[L]), L)
    say('   3  along ℓ (n = 0 to 9): H₁ = %s, H₃ = %s, H₅ = %s' % (H[1], H[3], H[5]))
    say('      H₅ ≡ T·H₃ − N·H₁ holds for (T, N) = %s;  PREDICTED from a_%d of the pair: T = %d, N = %d%s' % (
        sols or 'none', L, Tp, Np, ' (the two numbers: %s)' % lam(L) if split else ''))
    say('      one constant from level 1 (H₃ ≡ A·H₁): %s;  from level 3 (H₅ ≡ A·H₃): %s' % (A13 or 'none', A35 or 'none'))
    d3 = pow(24, -1, L ** 3); rows = []
    for Q in (5, 7, 11, 13):
        dl = -kron(-3, Q) * pow(Q, -2, L) % L; g = {1: set(), -1: set()}; cnt = {1: 0, -1: 0}; j = 0
        while min(cnt.values()) < 4 and j < 40:
            Ns = L ** 3 * j + d3; M = 24 * Ns - 1; j += 1
            if M % Q == 0: continue
            c0 = pm(Ns, L)
            if c0 == 0: continue
            e = kron(M, Q); cnt[e] += 1; g[e].add(pm(Q * Q * Ns - (Q * Q - 1) // 24, L) * pow(c0, -1, L) % L)
        if min(cnt.values()) < 4: rows.append('Q = %d: too few values' % Q); continue
        if len(g[1]) == 1 and len(g[-1]) == 1:
            l_p, l_m = (min(g[1]) - dl) % L, (min(g[-1]) + dl) % L
            which = [i + 1 for i, lm in enumerate(lam(Q))] if split else []
            which = [i for i in which if lam(Q)[i - 1] == l_p]
            rows.append('Q = %d: one each, λ = %d%s%s' % (Q, l_p, '' if l_p == l_m else ' / %d (the classes disagree)' % l_m,
                ', rule %s' % which if split else ''))
        else: rows.append('Q = %d: %d/%d values' % (Q, len(g[1]), len(g[-1])))
    say('      along Q at levels ≥ 3 (ℓ³ dividing 24N − 1), one multiplier per class?  %s' % ';  '.join(rows))
    say('      (%.0f s)' % (time.time() - t))

say(); say('THE CONTROL: ONE FORM (13 to 31). Is the constant of the levels, H₃ ≡ A·H₁ (n = 0 to 9), equal to (3|ℓ)·a_ℓ(f)?')
say('   (the first digit of each constant in partitions_constants.txt: 13 → 6, 17 → 6, 19 → 10, 23 → 5, 29 → 3, 31 → 30)')
for L in (13, 17, 19, 23, 29, 31):
    cls = sign_class(L)
    H1 = [pm(L * n + pow(24, -1, L), L) for n in range(10)]; H3 = [pm(L ** 3 * n + pow(24, -1, L ** 3), L) for n in range(10)]
    A = [A for A in range(L) if all((H3[n] - A * H1[n]) % L == 0 for n in range(10))]
    pr = [(kron(3, L) * int(pari.mfcoefs(F, L)[L]) % L) if int(pari.poldegree(P)) == 1 else 'degree %d' % int(pari.poldegree(P)) for F, P in cls]
    say('   ℓ = %2d: A = %s;  (3|ℓ)·a_ℓ(f) on lane ℓ = %s (forms in the class: %d);  %s' % (
        L, A or 'none', pr, len(cls), 'equal' if len(A) == 1 and pr == A else 'NOT equal'))
say(); say('seconds: %.0f' % (time.time() - t0))
