# A dig, 12 Nov 2026: THE LOOPS OF EVERY PRIME TO 700. partitions_loops.txt found the tower's loop on lane ℓ (the order of
# its constants there) as long as it can be for every one-constant prime from 7 to 29: the constant is a primitive root,
# six in a row, about a 1-in-300 chance if the constants are random; 5 and 31 have the shortest loop instead (−1). From 37
# on, the twins' constants were primitive roots less often than chance (1 of 8).
# PREDICTED before the run (derived): the constants on lane ℓ behave like random units, so the number of primitive roots
# among the constants on lane ℓ is close to Σ φ(ℓ − 1)/(ℓ − 1), and among the pieces of degree 2 the number that generate
# the field with ℓ² elements is close to Σ φ(ℓ² − 1)/(ℓ² − ℓ); the run from 7 to 29 was luck.
# Method (forms6_flint.py, FLINT only): the class at weight ℓ − 3 on lane ℓ (exactly, class_padic, where forms of lower
# level interfere), M = (3|ℓ)·T_ℓ on it (its roots are the constants on lane ℓ, STATEMENT.md A3), the pieces of its
# characteristic polynomial; for a piece of degree 1 its root's order in F_ℓ^×; for degree 2 its root's order in
# F_(ℓ²)^×; for the whole tower, when every piece has degree ≤ 2, the loop c₁ = the lcm. COUNTED (orders), DERIVED (shares).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_loops_survey.py LO HI [--workers=4] [--keep]
import os, sys, time, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def primes_of(n):
    out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: out.append(n)
    return out


def order_mod(g, L):
    """The order of x in (F_ℓ[x]/g)^×, g irreducible of degree 1 or 2 over F_ℓ."""
    from flint import nmod_poly
    q1 = L ** g.degree() - 1; x = nmod_poly([0, 1], L); o = q1
    for r in primes_of(q1):
        while o % r == 0 and pow(x, o // r, g) == 1: o //= r
    return o


def one(L):
    import flint
    from flint import nmod_mat, nmod_poly
    import forms6_flint as F
    flint.ctx.threads = 2
    try:
        import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
    except Exception: pass
    t = time.time(); k = L - 3; n = k + 1; dexp = (L + 12) // 24; N = 4
    H = F.hecke_matrices(k, L ** N, [('U', 2), ('U', 3), ('T', L)])
    U2, U3, TL = (F.reduce_mat(H[op], L) for op in (('U', 2), ('U', 3), ('T', L)))
    c2 = leg(2, L) * pow(2, (k - 2) // 2, L) % L; c3 = leg(3, L) * pow(3, (k - 2) // 2, L) % L
    K1, r1 = F.generalized_kernel(U2 - F.identity(n, L) * c2); X1, R1 = F.normalized_basis(K1, r1, L)
    K2, r0 = F.generalized_kernel(F.restrict(U3, X1, R1, r1) - F.identity(r1, L) * c3)
    Xc, Rc = F.normalized_basis(X1 * nmod_mat(r1, r0, [int(K2[i, j]) for i in range(r1) for j in range(r0)], L), r0, L)
    if r0 != dexp:
        W, RW, r, Wn, Rn, d, Np = F.class_padic(H, k, L, N, leg(2, L), leg(3, L))
        Xc, Rc, r0 = F.reduce_mat(Wn, L), Rn, d
    M = F.restrict(TL, Xc, Rc, r0) * (leg(3, L) % L)
    fac = M.charpoly().factor()[1]
    pieces = []
    for g, e in fac:
        if g.degree() <= 2: pieces.append((g.degree(), e, order_mod(g, L), int(g.coeffs()[0]) if g.degree() == 1 else None))
        else: pieces.append((g.degree(), e, None, None))
    return L, r0, dexp, pieces, time.time() - t


if __name__ == '__main__':
    import multiprocessing as mp
    from sympy import primerange
    lo, hi = [int(a) for a in sys.argv[1:] if a.isdigit()][:2]
    WK = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--workers=')), 4); KEEP = '--keep' in sys.argv
    fh = open(os.path.join(HERE, 'partitions_loops_survey_%d_%d.txt' % (lo, hi)), 'x', encoding='utf-8') if KEEP else None

    def say(*a):
        line = ' '.join(str(x) for x in a); print(line, flush=True)
        if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())

    PS = list(primerange(max(lo, 13), hi + 1)); t0 = time.time()
    say('primes %d to %d (%d), %d workers; the constants on lane ℓ = the roots of (3|ℓ)·T_ℓ on the class' % (lo, hi, len(PS), WK))
    n1 = p1 = 0; e1 = 0.0; n2 = p2 = 0; e2 = 0.0; full = 0; nfull = 0
    with mp.Pool(WK) as pool:
        for L, d, dexp, pieces, dt in pool.imap(one, PS):
            parts = []; loop = 1; known = True
            for deg, e, o, root in sorted(pieces):
                if deg == 1:
                    n1 += 1; e1 += sum(1 for a in range(1, L) if math.gcd(a, L - 1) == 1) / (L - 1)
                    if o == L - 1: p1 += 1
                    parts.append('%d (order %d of %d%s)' % ((-root) % L, o, L - 1, ', primitive' if o == L - 1 else ''))
                elif deg == 2:
                    q1 = L * L - 1; ph = q1
                    for r in primes_of(q1): ph = ph // r * (r - 1)
                    n2 += 1; e2 += ph / (L * L - L)
                    if o == q1: p2 += 1
                    parts.append('a pair (order %d of %d%s)' % (o, q1, ', primitive' if o == q1 else ''))
                else: parts.append('a piece of degree %d' % deg); known = False
                if o: loop = math.lcm(loop, o)
                if e > 1: parts[-1] += ' ×%d' % e
            most = 1
            for deg, e, o, root in pieces: most = math.lcm(most, L ** deg - 1)
            if known:
                nfull += 1; full += loop == most
            say('   ℓ = %3d: %2d form%s (⌊(ℓ+12)/24⌋ = %d): %s%s  (%.0f s)' % (
                L, d, ' ' if d == 1 else 's', dexp, '; '.join(parts),
                ('  →  the loop c₁ = %d (%s)' % (loop, 'the longest possible' if loop == most else 'longest possible %d' % most)) if known else '', dt))
    say()
    say('constants on lane ℓ: %d; primitive roots: %d counted, %.1f expected (derived)' % (n1, p1, e1))
    say('pairs (degree 2): %d; generating the field with ℓ² elements: %d counted, %.1f expected (derived)' % (n2, p2, e2))
    say('towers whose pieces all have degree ≤ 2: %d; with the longest possible loop: %d' % (nfull, full))
    say('seconds: %.0f' % (time.time() - t0))
