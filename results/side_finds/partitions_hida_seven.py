# A dig, 11 Nov 2026: THE CONSTANT OF 7, A PRIME WITH LAWS. For 13 to 31 each constant is (3|ℓ)·u(−2), u the U_ℓ-number of
# the tower form's ordinary family (partitions_hida_digits.txt). For 7, Watson's laws add one factor of 7 per two levels, so
# the unnormalized two-level step is 7·α₇ and the family must have slope 1 (U_7-number divisible by 7 exactly once). By the
# Gouvêa–Mazur picture (my reading) a slope-1 family is one analytic piece only over weights that agree modulo 6·7 = 42, so
# the weights to use are k = 40 + 42·x (40 ≡ −2 modulo 42, also the weight where E7's lane-49 law sat: 2 + 2·19).
# Each member: the level-6 newform in 7's tower system (a_Q ≡ a_Q(f) on lane 7 for every prime Q ≤ 60 not dividing 42, f the
# weight-16 form with a₂ = +128, a₃ = −2187, partitions_newforms.txt) whose a_7 is divisible by 7 exactly once.
# PREDICTED before the run: (3|7)·u(−2)/7 agrees with the 14 settled digits of α₇ (partitions_constants.txt: 5, 6, 1, 1, 2,
# 1, 3, 1, 0, 6, 5, 3, 1, 1) as far as Newton's series reaches. A scratch look found the first digit at weights 40 and 82.
# Also shown: the members at the weights 16 to 46, where (3|7)·a_7/7 on lane 7 climbs by one per step of 6.
# No partition number is used. PARI 2.15.4 (cypari2) in WSL. COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_hida_seven.py [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(8 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
DIG = next((int(a[9:]) for a in sys.argv[1:] if a.startswith('--digits=')), 0)   # a later run: print this many digits
fh = open(os.path.join(HERE, 'partitions_hida_seven%s.txt' % ('_%d' % DIG if DIG else '')), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); L = 7; QM = 60; X = 6; PREC = X + 10   # weights 40 to 292: in the trials PARI crashed at weight 334
SETTLED = [5, 6, 1, 1, 2, 1, 3, 1, 0, 6, 5, 3, 1, 1]
PR = [int(q) for q in pari('primes([5, %d])' % QM) if int(q) != L]
kron = lambda a, q: int(pari.kronecker(a, q))


def lane(x, e):
    if pari.type(x) == 't_PADIC': x = pari.truncate(x)
    return int(pari.lift(pari.Mod(x, L ** e)))


f = None
mf = pari('mfinit([6, 16], 0)')
for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
    c = pari.mfcoefs(F, QM)
    if int(pari.poldegree(P)) == 1 and int(c[2]) == 128 and int(c[3]) == -2187: f = [int(x) for x in c]


def members(k):                                       # (slope, a_7/7^slope on lane 7^PREC) of each member of 7's system
    m = pari('mfinit([6, %d], 0)' % k); out = []
    for G, PG in zip(pari.mfeigenbasis(m), pari.mffields(m)):
        g = pari.mfcoefs(G, QM); deg = int(pari.poldegree(PG))
        try: embs = [None] if deg == 1 else list(pari.polrootspadic(PG, L, 400))   # high: the coefficients' denominators
        except cypari2.PariError: continue                                               # hold large powers of 7 (weight 334)
        for r in embs:
            ev = (lambda Q: g[Q]) if r is None else (lambda Q, r=r: pari.subst(pari.lift(g[Q]), 'y', r))
            try:
                if all(lane(ev(Q), 1) == f[Q] % L for Q in PR):
                    a7 = ev(7); v = int(pari.valuation(a7, L))
                    out.append((v, lane(a7 / L ** v, PREC)))
            except cypari2.PariError: continue
    return out


say('the members of 7\'s system at the weights 16 to 46: (3|7)·a_7/7 on lane 7, and the slope')
say('   ' + ';  '.join('weight %d: %s' % (k, ', '.join('%d (slope %d)' % ((kron(3, L) * w) % L, v) for v, w in members(k)) or 'none')
                       for k in (16, 22, 28, 34, 40, 46)))
u = []
for x in range(X + 1):
    k = 40 + 42 * x; t = time.time(); ms = members(k); one = [w for v, w in ms if v == 1]
    say('   weight %d: %d members of the system, slopes %s  (%.0f s)' % (k, len(ms), sorted(v for v, w in ms), time.time() - t))
    if len(one) != 1: say('   not exactly one member of slope 1; stopped'); break
    u.append(one[0])                                  # u/7 = a_7/7 (the slope-1 root equals a_7 to 7^(k−2))
if len(u) == X + 1:
    M = L ** PREC
    diffs = [sum((-1) ** (j - i) * int(pari.binomial(j, i)) * u[i] for i in range(j + 1)) % M for j in range(X + 1)]
    vals = [int(pari.valuation(d, L)) if d else PREC for d in diffs]
    held = sum(int(pari.binomial(X, j)) * diffs[j] for j in range(X)) % M
    agree = next((e for e in range(PREC) if (held - u[X]) % L ** (e + 1)), PREC)
    w = kron(3, L) * sum((-1) ** j * diffs[j] for j in range(X + 1)) % M; dg = [(w // L ** e) % L for e in range(DIG or X + 1)]
    if DIG: say('(--digits=%d: the differences shrink faster than one power per step, so more digits are shown; how far they' % DIG)
    if DIG: say(' agree with the partition numbers is measured here, not predicted)')
    same = next((i for i in range(len(dg)) if dg[i] != SETTLED[i]), len(dg))
    say(); say('differences of a_7/7 divisible by 7 to the powers %s' % vals)
    say('from x = 0 to %d the series foretells the member at x = %d (weight %d) to %d digits' % (X - 1, X, 40 + 42 * X, agree))
    say('(3|7)·u(−2)/7, digits lowest first: %s' % dg)
    say('the constant of 7 from the partition numbers:    %s  →  the first %d digits agree' % (SETTLED[:len(dg)], same))
say(); say('seconds: %.0f' % (time.time() - t0))
