# A dig, 12 Nov 2026: THE CONSTANTS OF 5 AND 11, THE OTHER PRIMES WITH LAWS. 7's constant is the number of its slope-1 family
# at weight −2 (partitions_hida_seven.txt, 13 of 14 digits). For 5 and 11 the laws (Watson, Atkin) add two factors of ℓ per
# two levels, so the family should have slope 2: its U_ℓ-number divisible by ℓ exactly twice, and α_ℓ = (3|ℓ)·u(−2)/ℓ².
# The tower's system on lane ℓ: for 5, a_Q ≡ Q²(1 + Q) (the README's "1 + Q" rule, λ = (3|Q)·a_Q); for 11, the weight-28
# level-6 form with a₂ = −8192, a₃ = +3¹³ (partitions_newforms.txt, matched with no extra factor).
#   1  THE SURVEY: every weight k ≡ −2 (mod ℓ − 1) up to 98 (5) or 148 (11): the slopes of the members of the system, and
#      (3|ℓ)·a_ℓ/ℓ² on lane ℓ for the members of slope 2. A scratch survey had shown this digit repeating with the weight
#      modulo (ℓ − 1)·ℓ and equal to α_ℓ's first digit at the weights ≡ −2 there.
#   2  THE DIGITS: Newton's series over k = −2 + (ℓ − 1)·ℓ·(x + 1) (5: 18 to 258; 11: 108, 218, 328, PARI's limit here),
#      at x = −1 (weight −2). PREDICTED before the run: they agree with the settled digits of partitions_constants.txt
#      (5: 33 digits; 11: 9, 3, 2, 0) as far as the series reaches. (A trial on 5, before this kept run: 19 digits.)
# No partition number is used. PARI 2.15.4 (cypari2) in WSL. COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_hida_laws.py [ℓ …] [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(12 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_hida_laws.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); QM = 60
SETTLED = {5: [4, 2, 3, 1, 0, 2, 0, 3, 0, 1, 1, 0, 1, 4, 0, 1, 3, 2, 3, 0, 3, 1, 0, 4, 3, 0, 4, 1, 2, 4, 2, 3, 0],
           11: [9, 3, 2, 0]}
XS = {5: 12, 11: 2}
SURVEY = {5: range(18, 99, 4), 11: range(28, 149, 10)}
Ls = [int(a) for a in sys.argv[1:] if a.isdigit()] or [5, 11]
kron = lambda a, q: int(pari.kronecker(a, q))


def members(L, ref, PR, k, prec):                          # (slope, a_ℓ as an ℓ-adic number) of each member of the system
    m = pari('mfinit([6, %d], 0)' % k); out = []
    for G, PG in zip(pari.mfeigenbasis(m), pari.mffields(m)):
        g = pari.mfcoefs(G, QM); deg = int(pari.poldegree(PG))
        try: embs = [None] if deg == 1 else list(pari.polrootspadic(PG, L, 400))
        except cypari2.PariError: continue
        for r in embs:
            ev = (lambda Q: g[Q]) if r is None else (lambda Q, r=r: pari.subst(pari.lift(g[Q]), 'y', r))
            try:
                if all(int(pari.lift(pari.Mod(ev(Q), L))) == ref[Q] for Q in PR):
                    a = ev(L) + pari('O(%d^%d)' % (L, prec)); out.append((int(pari.valuation(a, L)), a))
            except cypari2.PariError: continue
    return out


for L in Ls:
    t = time.time(); PR = [int(q) for q in pari('primes([5, %d])' % QM) if int(q) != L]; X = XS[L]; PREC = 2 * X + 12
    lane = lambda x, e: int(pari.lift(pari.Mod(x, L ** e)))
    if L == 5: ref = {Q: Q * Q * (1 + Q) % 5 for Q in PR}
    else:
        m28 = pari('mfinit([6, 28], 0)')
        for F, P in zip(pari.mfeigenbasis(m28), pari.mffields(m28)):
            c = pari.mfcoefs(F, QM)
            if int(pari.poldegree(P)) == 1 and int(c[2]) == -8192 and int(c[3]) == 1594323: ref = {Q: int(c[Q]) % L for Q in PR}
    say(); say('ℓ = %d' % L)
    rows = []
    for k in SURVEY[L]:
        ms = members(L, ref, PR, k, 8)
        rows.append('%d: %s' % (k, ', '.join('slope %d%s' % (v, ' → %d' % (kron(3, L) * lane(a / L ** 2, 1) % L) if v == 2 else '')
                                             for v, a in sorted(ms, key=lambda z: z[0])) or 'none'))
    say('   1  the survey (weights ≡ −2 mod %d; for slope 2, (3|ℓ)·a_ℓ/ℓ² on lane ℓ):  %s' % (L - 1, ';  '.join(rows)))
    step = (L - 1) * L; u = []; rows = []
    for x in range(X + 1):
        k = -2 + step * (x + 1); ms = members(L, ref, PR, k, PREC + 2); two = [a for v, a in ms if v == 2]
        rows.append('%d: slopes %s' % (k, sorted(v for v, a in ms)))
        if len(two) != 1: say('   weight %d: %d members of slope 2; stopped' % (k, len(two))); break
        a = two[0]; z = a
        for _ in range(PREC + 4): z = a - pari('%d^%d' % (L, k - 1)) / z            # the slope-2 root of z² − a·z + ℓ^(k−1)
        u.append(lane(z / L ** 2, PREC))
    say('   2  the weights used: %s  (%.0f s)' % (';  '.join(rows), time.time() - t))
    if len(u) < X + 1: continue
    M = L ** PREC
    diffs = [sum((-1) ** (j - i) * int(pari.binomial(j, i)) * u[i] for i in range(j + 1)) % M for j in range(X + 1)]
    vals = [int(pari.valuation(d, L)) if d else PREC for d in diffs]
    w = kron(3, L) * sum((-1) ** j * diffs[j] for j in range(X + 1)) % M          # x = −1 is weight −2
    nd = min(PREC - 2, len(SETTLED[L]) + 2); dg = [(w // L ** e) % L for e in range(nd)]
    same = next((i for i in range(min(nd, len(SETTLED[L]))) if dg[i] != SETTLED[L][i]), min(nd, len(SETTLED[L])))
    say('      differences of u/ℓ² divisible by ℓ to the powers %s' % vals)
    if X >= 2:
        held = sum(int(pari.binomial(X, j)) * diffs[j] for j in range(X)) % M
        say('      from the first %d weights the series foretells the last (weight %d) to %d digits' % (
            X, -2 + step * (X + 1), next((e for e in range(PREC) if (held - u[X]) % L ** (e + 1)), PREC)))
    say('      (3|ℓ)·u(−2)/ℓ², digits lowest first: %s' % dg)
    say('      the constant from the partition numbers:   %s  →  the first %d digits agree' % (SETTLED[L][:nd], same))
say(); say('seconds: %.0f' % (time.time() - t0))
