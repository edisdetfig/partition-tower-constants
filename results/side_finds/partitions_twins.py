# A dig, 11 Nov 2026: THE TWO CONSTANTS OF A SPLIT TOWER. partitions_split.txt found, on lane ℓ, that from 37 on the levels two
# apart obey a two-number law H_(b+4) ≡ T·H_(b+2) − N·H_b, T and N the sum and product of (3|ℓ)·a_ℓ of the twin forms f, f′.
# partitions_hida_digits.txt found each single constant (13 to 31) to be (3|ℓ)·u(−2), u the U_ℓ-number of the form's ordinary
# ℓ-adic family. Here, for 37, 41, 43 (ℓ splits in the twins' field, so each twin is an ordinary form over the ℓ-adic
# integers): both constants α₁, α₂ from the twins' families (weights ℓ − 3 + (ℓ − 1)·x, x = 0 to 8, Newton's series at x = −1),
# and the law with T = α₁ + α₂, N = α₁·α₂ on lanes beyond ℓ, on exact partition numbers H_b(0) = p(δ_b), b = 1 to 7.
# PREDICTED before the run (derived: on lane ℓ the level-1 tower lies in the twins' two-dimensional space, so the other
# eigen-parts are divisible by ℓ there and pick up at least one more factor per two levels): from levels 1, 3, 5 the law
# holds on lane ℓ, from levels 3, 5, 7 on lane ℓ² at least (levels 2, 4, 6, which obey the same law: on lane ℓ at least).
# (A trial on 37 alone, before this kept run, gave ℓ¹ and ℓ² for the odd levels.)
# No partition number is used for α₁, α₂. PARI 2.15.4 (cypari2) in WSL. COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_twins.py [ℓ …] [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(8 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
NAMED = [a for a in sys.argv[1:] if a.isdigit()]                          # a run for chosen primes keeps its own file
fh = open(os.path.join(HERE, 'partitions_twins%s.txt' % ''.join('_' + a for a in NAMED)), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


X = next((int(a[4:]) for a in sys.argv[1:] if a.startswith('--x=')), 8)  # weights x = 0 to X (43 ran with 4: PARI ran out of
t0 = time.time(); QM = 60; PREC = X + 8                                   # its 8 GB stack at weight 376, crashed at 8 MB C stack)
kron = lambda a, q: int(pari.kronecker(a, q))
Ls = [int(a) for a in sys.argv[1:] if a.isdigit()] or [37, 41, 43]


def lane(x, L, e):
    if pari.type(x) == 't_PADIC': x = pari.truncate(x)
    return int(pari.lift(pari.Mod(x, L ** e)))


def val(x, L, cap): return int(pari.valuation(x, L)) if x else cap


for L in Ls:
    t = time.time(); k0 = L - 3; PR = [int(q) for q in pari('primes([5, %d])' % QM) if int(q) != L]
    s2, s3 = kron(2, L), kron(3, L); mf = pari('mfinit([6, %d], 0)' % k0); twin = None
    for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
        c = pari.mfcoefs(F, max(QM, L))
        if (int(pari.lift(c[2])) > 0) == (s2 > 0) and (int(pari.lift(c[3])) > 0) == (s3 > 0): twin = (c, P)
    c, P = twin; roots = list(pari.polrootspadic(P, L, PREC + 4)); assert len(roots) == 2
    targets = [{Q: lane(pari.subst(pari.lift(c[Q]), 'y', r), L, 1) for Q in PR} for r in roots]
    alphas = []
    for i, tg in enumerate(targets):
        u = []
        for x in range(X + 1):
            k = k0 + (L - 1) * x; m = pari('mfinit([6, %d], 0)' % k); found = []
            for G, PG in zip(pari.mfeigenbasis(m), pari.mffields(m)):
                g = pari.mfcoefs(G, max(QM, L)); deg = int(pari.poldegree(PG))
                try: embs = [None] if deg == 1 else list(pari.polrootspadic(PG, L, PREC + 4))
                except cypari2.PariError: continue
                for r in embs:
                    ev = (lambda Q: g[Q]) if r is None else (lambda Q, r=r: pari.subst(pari.lift(g[Q]), 'y', r))
                    try:
                        if all(lane(ev(Q), L, 1) == tg[Q] for Q in PR) and lane(ev(L), L, 1) != 0:
                            a = ev(L) + pari('O(%d^%d)' % (L, PREC)); z = a
                            for _ in range(PREC + 2): z = a - pari('%d^%d' % (L, k - 1)) / z
                            found.append(lane(z, L, PREC))
                    except cypari2.PariError: continue
            if len(found) != 1: say('ℓ = %d, twin %d, weight %d: %d members found; stopped' % (L, i + 1, k, len(found))); break
            u.append(found[0])
        if len(u) < X + 1: break
        M = L ** PREC
        diffs = [sum((-1) ** (j - q) * int(pari.binomial(j, q)) * u[q] for q in range(j + 1)) % M for j in range(X + 1)]
        held = sum(int(pari.binomial(X, j)) * diffs[j] for j in range(X)) % M
        agree = next((e for e in range(PREC) if (held - u[X]) % L ** (e + 1)), PREC)
        al = kron(3, L) * sum((-1) ** j * diffs[j] for j in range(X + 1)) % M; alphas.append(al)
        say('ℓ = %d, twin %d (y ↦ %d on lane %d): differences divisible by ℓ to %s; the member at x = %d foretold to %d digits' % (
            L, i + 1, lane(roots[i], L, 1), L, [val(d, L, PREC) for d in diffs], X, agree))
        say('   its constant (3|ℓ)·u(−2), digits lowest first: %s' % [(al // L ** e) % L for e in range(X + 1)])
    if len(alphas) < 2: continue
    K = X + 1; MK = L ** K; a1, a2 = alphas[0] % MK, alphas[1] % MK; T, N = (a1 + a2) % MK, a1 * a2 % MK
    say('   T = α₁ + α₂ and N = α₁·α₂ on lane %d: %d, %d  (%.0f s for the forms)' % (L, T % L, N % L, time.time() - t))
    H = {}
    for b in range(1, 8):
        tb = time.time(); d = pow(24, -1, L ** b); g = pari.numbpart(d); H[b] = int(g) % MK
        say('   p(δ_%d) = p(%d): %d digits, exact  (%.0f s)' % (b, d, int(pari.logint(g, 10)) + 1, time.time() - tb))
    for b in (1, 2, 3):
        r = (H[b + 4] - T * H[b + 2] + N * H[b]) % MK
        one = [val((H[b + 2] - al * H[b]) % MK, L, K) for al in (a1, a2)]
        say('   levels %d, %d, %d: H_(b+4) − T·H_(b+2) + N·H_b is divisible by ℓ^%d;  one constant alone, H_(b+2) − αᵢ·H_b: by ℓ^%s' % (
            b, b + 2, b + 4, val(r, L, K), one))
say(); say('seconds: %.0f' % (time.time() - t0))
