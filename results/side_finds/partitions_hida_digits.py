# A dig, 11 Nov 2026: THE CONSTANTS FROM MODULAR FORMS ALONE. partitions_split.txt and partitions_hida.txt found the first two
# digits of each constant of partitions_constants.txt in f's ordinary ℓ-adic family (f: the tower's level-6 newform of weight
# ℓ − 3). Suppose the constant is (3|ℓ)·u(−2), u(k) the family's U_ℓ-number at weight k. On the weights k = ℓ − 3 + (ℓ − 1)·x
# Hida's theory makes u a power series in (1 + ℓ)^((ℓ−1)x) − 1 with ℓ-adic integer coefficients, so the n-th difference of u in
# x is divisible by ℓⁿ, and Newton's series from x = 0, …, X gives u at x = −1 (weight −2) right to about X + 1 digits.
# PREDICTED before the run: it agrees with every settled digit of partitions_constants.txt (lowest first; 13: 6, 3, 1, 2, 5, 8;
# 17: 6, 12, 12, 7, 7; 19: 10, 2, 10, 8, 10; 23: 5, 20, 8, 4; 29: 3, 14; 31: 30, 6). Checks: the divisibility of the
# differences, and the series from x = 0 to X − 1 foretells the member at x = X.
# u(k) is the unit root of z² − a_ℓ·z + ℓ^(k−1), a_ℓ of the one member of f's family at weight k (a_Q ≡ a_Q(f) on lane ℓ for
# every prime Q ≤ 60 not dividing 6ℓ, and a_ℓ a unit). No partition number is used. PARI 2.15.4 (cypari2) in WSL. COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_hida_digits.py [ℓ …] [--keep]
import os, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(8 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
X = next((int(a[4:]) for a in sys.argv[1:] if a.startswith('--x=')), 10)   # weights x = 0 to X (a later run used 12)
fh = open(os.path.join(HERE, 'partitions_hida_digits%s.txt' % ('' if X == 10 else '_x%d' % X)), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); QM = 60; PREC = X + 8
kron = lambda a, q: int(pari.kronecker(a, q))
SETTLED = {13: [6, 3, 1, 2, 5, 8], 17: [6, 12, 12, 7, 7], 19: [10, 2, 10, 8, 10], 23: [5, 20, 8, 4], 29: [3, 14], 31: [30, 6]}
Ls = [int(a) for a in sys.argv[1:] if a.isdigit()] or [13, 17, 19, 23, 29, 31]


def lane(x, L, e):                                                          # an ℓ-adic integer (or rational) on lane ℓᵉ
    if pari.type(x) == 't_PADIC': x = pari.truncate(x)
    return int(pari.lift(pari.Mod(x, L ** e)))


def digits(v, L, n): return [(v // L ** i) % L for i in range(n)]


for L in Ls:
    t = time.time(); k0 = L - 3; PR = [int(q) for q in pari('primes([5, %d])' % QM) if int(q) != L]
    s2, s3 = kron(2, L), kron(3, L); f = None
    mf = pari('mfinit([6, %d], 0)' % k0)
    for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
        c = pari.mfcoefs(F, max(QM, L))
        if int(pari.poldegree(P)) == 1 and (int(c[2]) > 0) == (s2 > 0) and (int(c[3]) > 0) == (s3 > 0): f = [lane(c[Q], L, 1) if Q in PR else 0 for Q in range(QM + 1)]
    u, rows = [], []
    for x in range(X + 1):
        k = k0 + (L - 1) * x; mf = pari('mfinit([6, %d], 0)' % k); found = []
        for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
            c = pari.mfcoefs(F, max(QM, L)); deg = int(pari.poldegree(P))
            try: embs = [None] if deg == 1 else list(pari.polrootspadic(P, L, PREC + 4))
            except cypari2.PariError: continue
            for r in embs:
                ev = (lambda Q: c[Q]) if r is None else (lambda Q, r=r: pari.subst(pari.lift(c[Q]), 'y', r))
                try:
                    if all(lane(ev(Q), L, 1) == f[Q] for Q in PR) and lane(ev(L), L, 1) != 0:
                        a = ev(L) + pari('O(%d^%d)' % (L, PREC)); z = a
                        for _ in range(PREC + 2): z = a - pari('%d^%d' % (L, k - 1)) / z     # the unit root, by iteration
                        found.append(lane(z, L, PREC))
                except cypari2.PariError: continue
        if len(found) != 1: say('ℓ = %d, weight %d: %d members of f\'s family found; stopped' % (L, k, len(found))); break
        u.append(found[0]); rows.append('%d' % k)
    if len(u) < X + 1: continue
    M = L ** PREC
    diffs = [sum((-1) ** (j - i) * int(pari.binomial(j, i)) * u[i] for i in range(j + 1)) % M for j in range(X + 1)]
    val = [int(pari.valuation(d, L)) if d else PREC for d in diffs]
    at = lambda xx, n: sum(int(pari.binomial(xx, j)) * diffs[j] for j in range(n)) % M     # Newton's series with n terms
    held = at(X, X); agree = next((i for i in range(PREC) if (held - u[X]) % L ** (i + 1)), PREC)
    w = kron(3, L) * at(-1, X + 1) % M; dg = digits(w, L, X + 1); st = SETTLED[L]
    say(); say('ℓ = %d: weights %s (%.0f s)' % (L, ', '.join(rows), time.time() - t))
    say('   the differences of u are divisible by ℓ to the powers %s (needed: at least 0, 1, 2, …)' % val)
    say('   from x = 0 to %d the series foretells the member at x = %d (weight %d) to %d digits' % (X - 1, X, k0 + (L - 1) * X, agree))
    say('   (3|ℓ)·u(−2), digits lowest first: %s' % dg)
    say('   the constant from the partition numbers (partitions_constants.txt): %s  →  %s' % (st,
        'ALL %d settled digits agree' % len(st) if dg[:len(st)] == st else 'differ at digit %d' % next(i for i in range(len(st)) if dg[i] != st[i])))
say(); say('seconds: %.0f' % (time.time() - t0))
