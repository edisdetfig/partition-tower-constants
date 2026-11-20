# A dig, 12 Nov 2026: WHICH PRIMES BREAK THE STAIRCASE? Boylan–Webb computed the rank of the partition tower's stable space to
# be ⌊(ℓ + 12)/24⌋ for every prime up to 1297 except 607. partitions_flint_607.txt refuted my first reading (a non-ordinary
# form): 607's class has all 25 forms, none non-ordinary, but T_607 on the class has a REPEATED root on lane 607 (two
# conjugate forms share their number at 607), while 599, 601, 613 have none (cells/rep607.py, scratch).
# PREDICTED before this sweep (derived): the levels' rank drops below ⌊(ℓ + 12)/24⌋ exactly where T_ℓ on the tower's class
# has, on lane ℓ, a zero root (a non-ordinary form) or a repeated root (two constants equal on lane ℓ, so the tower's
# iterates span one direction fewer). By a 1-in-ℓ chance per prime the expected number of such primes up to 1297 is about
# Σ 1/ℓ ≈ 1; Boylan–Webb found exactly one, 607. So the sweep should flag 607 and no other prime.
# Method (forms6_flint.py, FLINT only, modulo ℓ): the class at weight ℓ − 3 (joint generalized eigenspace of U₂, U₃ at the
# class's Atkin–Lehner values); T_ℓ on it; its characteristic polynomial: the multiplicity of 0 and whether it is
# squarefree; the class's dimension against ⌊(ℓ + 12)/24⌋. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_exceptions.py LO HI [--workers=6] [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def one(L):
    import flint
    from flint import nmod_mat
    import forms6_flint as F
    flint.ctx.threads = 2
    try:
        import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
    except Exception: pass
    t = time.time(); k = L - 3; n = k + 1
    H = F.hecke_matrices(k, L, [('U', 2), ('U', 3), ('T', L)])
    c2 = leg(2, L) * pow(2, (k - 2) // 2, L) % L; c3 = leg(3, L) * pow(3, (k - 2) // 2, L) % L
    K1, r1 = F.generalized_kernel(H[('U', 2)] - F.identity(n, L) * c2); X1, R1 = F.normalized_basis(K1, r1, L)
    M3 = F.restrict(H[('U', 3)], X1, R1, r1); K2, r0 = F.generalized_kernel(M3 - F.identity(r1, L) * c3)
    Xc, Rc = F.normalized_basis(X1 * nmod_mat(r1, r0, [int(K2[i, j]) for i in range(r1) for j in range(r0)], L), r0, L)
    chi = F.restrict(H[('T', L)], Xc, Rc, r0).charpoly(); fac = chi.factor()[1]
    zero = sum(e for f, e in fac if f.degree() == 1 and int(f.coeffs()[0]) == 0)
    rep = sum(1 for f, e in fac if e > 1)
    return L, r0, (L + 12) // 24, zero, rep, sorted((f.degree(), e) for f, e in fac), time.time() - t


if __name__ == '__main__':
    import multiprocessing as mp
    lo, hi = [int(a) for a in sys.argv[1:] if a.isdigit()][:2]
    W = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--workers=')), 6); KEEP = '--keep' in sys.argv
    fh = open(os.path.join(HERE, 'partitions_flint_exceptions_%d_%d.txt' % (lo, hi)), 'x', encoding='utf-8') if KEEP else None

    def say(*a):
        line = ' '.join(str(x) for x in a); print(line, flush=True)
        if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())

    from sympy import primerange
    PS = list(primerange(max(lo, 13), hi + 1)); t0 = time.time()
    say('primes %d to %d (%d of them), %d workers' % (lo, hi, len(PS), W))
    flagged = []
    with mp.Pool(W) as pool:
        for L, r0, dexp, zero, rep, fac, dt in pool.imap(one, PS):
            mark = []
            if r0 != dexp: mark.append('CLASS %d ≠ %d' % (r0, dexp))
            if zero: mark.append('ZERO ROOT ×%d' % zero)
            if rep: mark.append('REPEATED ROOT')
            if mark: flagged.append(L)
            say('   ℓ = %4d: class %2d (⌊(ℓ+12)/24⌋ = %2d); T_ℓ on lane ℓ: zero roots %d, repeated factors %d  %s  (%.0f s)' % (
                L, r0, dexp, zero, rep, ' ← ' + ', '.join(mark) if mark else '', dt))
    say(); say('flagged: %s  (Boylan–Webb\'s exception up to 1297: 607)' % (flagged or 'none'))
    say('seconds: %.0f' % (time.time() - t0))
