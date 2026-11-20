# A dig, 12 Nov 2026: WHICH PRIMES BREAK THE STAIRCASE, with the class taken EXACTLY. The sweep on lane ℓ alone
# (partitions_flint_exceptions_37_700.txt) flagged 607, Boylan–Webb's one exception up to 1297, and also 499: there the
# class found on lane ℓ had 22 forms instead of ⌊(499 + 12)/24⌋ = 21, and T_499 on it a repeated root.
# PREDICTED before this run (derived): 499's extra form is a form of lower level (or an Eisenstein series) whose U₂, U₃
# values agree with the class's on lane 499 only (a level-raising congruence, as for 17 and 29), and since it is congruent
# to one of the class's forms it repeats that form's root of T_499. Taken exactly, 499's class has 21 forms and T_499 on it
# no repeated root and no zero root: 499 is not an exception, as Boylan–Webb found. 607 stays flagged (its class on lane ℓ
# already has the right 25 forms). Over all primes 37 to 1297 the exact test flags 607 alone.
# Method (forms6_flint.py, FLINT only): M_{ℓ−3}(Γ0(6)) modulo ℓ^N; the joint generalized eigenspace of U₂, U₃ for the
# class's values, lifted over the ℓ-adic integers (Fitting's splitting); inside it the forms whose U₂, U₃ are EXACTLY
# (2|ℓ)·2^((k−2)/2), (3|ℓ)·3^((k−2)/2) (an ℓ-adic kernel). T_ℓ on that class modulo ℓ: zero roots (non-ordinary forms),
# repeated factors; T_ℓ on the rest (the forms that agree only on lane ℓ). The ℓ-adic step runs only where the class on
# lane ℓ is larger than ⌊(ℓ + 12)/24⌋ (elsewhere it is the reduction of the exact class; --padic-all runs it everywhere).
# COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_exceptions_new.py (--range LO HI | ℓ …)
#          [--workers=6] [--N=4] [--padic-all] [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def factors(chi):
    return sorted((f.degree(), e) for f, e in chi.factor()[1])


def one(job):
    L, N, ALL = job
    import flint
    from flint import nmod_mat
    import forms6_flint as F
    flint.ctx.threads = 2
    try:
        import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
    except Exception: pass
    t = time.time(); k = L - 3; n = k + 1; dexp = (L + 12) // 24; m = L ** N
    H = F.hecke_matrices(k, m, [('U', 2), ('U', 3), ('T', L)])
    U2, U3, TL = (F.reduce_mat(H[op], L) for op in (('U', 2), ('U', 3), ('T', L)))
    c2 = leg(2, L) * pow(2, (k - 2) // 2, L) % L; c3 = leg(3, L) * pow(3, (k - 2) // 2, L) % L
    K1, r1 = F.generalized_kernel(U2 - F.identity(n, L) * c2); X1, R1 = F.normalized_basis(K1, r1, L)
    M3 = F.restrict(U3, X1, R1, r1); K2, r0 = F.generalized_kernel(M3 - F.identity(r1, L) * c3)
    Xc, Rc = F.normalized_basis(X1 * nmod_mat(r1, r0, [int(K2[i, j]) for i in range(r1) for j in range(r0)], L), r0, L)
    chi0 = F.restrict(TL, Xc, Rc, r0).charpoly()
    res = dict(L=L, r0=r0, dexp=dexp, lane=factors(chi0))
    if r0 != dexp or ALL:
        out = F.class_padic(H, k, L, N, leg(2, L), leg(3, L))
        W, RW, r, Wn, Rn, d, Np = out
        res.update(r=r, d=d, Np=Np)
        if Wn is not None:
            chin = F.restrict(TL, F.reduce_mat(Wn, L), Rn, d).charpoly()
            chiW = F.restrict(TL, F.reduce_mat(W, L), RW, r).charpoly()
            res.update(new=factors(chin), old=factors(chiW // chin) if r > d else [],
                       shared=(r > d and (chiW // chin).gcd(chin).degree() > 0))
            chi = chin
        else: chi = None
    else: chi = chi0; res.update(d=r0)
    if chi is not None:
        fac = chi.factor()[1]
        res.update(zero=sum(e for f, e in fac if f.degree() == 1 and int(f.coeffs()[0]) == 0), rep=sum(1 for f, e in fac if e > 1))
    res['sec'] = time.time() - t
    return res


if __name__ == '__main__':
    import multiprocessing as mp
    nums = [int(a) for a in sys.argv[1:] if a.isdigit()]
    W_ = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--workers=')), 6)
    N = next((int(a[4:]) for a in sys.argv[1:] if a.startswith('--N=')), 4)
    ALL = '--padic-all' in sys.argv; KEEP = '--keep' in sys.argv
    from sympy import primerange
    if "--range" in sys.argv:
        PS = list(primerange(max(nums[0], 13), nums[1] + 1)); tag = '%d_%d' % (nums[0], nums[1])
    else:
        PS = nums; tag = '_'.join(str(x) for x in nums)
    fh = open(os.path.join(HERE, 'partitions_flint_exceptions_new_%s%s.txt' % (tag, '_all' if ALL else '')), 'x', encoding='utf-8') if KEEP else None

    def say(*a):
        line = ' '.join(str(x) for x in a); print(line, flush=True)
        if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())

    t0 = time.time(); say('primes: %d (%d to %d), %d workers, modulo ℓ^%d%s' % (len(PS), PS[0], PS[-1], W_, N, ', the ℓ-adic step everywhere' if ALL else ''))
    flagged = []
    with mp.Pool(W_) as pool:
        for R in pool.imap(one, [(L, N, ALL) for L in PS]):
            L = R['L']; mark = []
            line = '   ℓ = %4d: on lane ℓ the class has %2d forms (⌊(ℓ+12)/24⌋ = %2d)' % (L, R['r0'], R['dexp'])
            if 'r' in R:
                line += '; exactly: %d of the %d (to ℓ^%d)' % (R['d'], R['r'], R['Np'])
                if R.get('old'): line += ', the other %d with T_ℓ factors %s%s' % (R['r'] - R['d'], R['old'], ' (shared with the class: the congruent copy)' if R['shared'] else ' (NOT shared)')
            if R.get('d') != R['dexp']: mark.append('CLASS %s ≠ %d' % (R.get('d'), R['dexp']))
            if R.get('zero'): mark.append('ZERO ROOT ×%d' % R['zero'])
            if R.get('rep'): mark.append('REPEATED ROOT')
            if mark: flagged.append(L)
            line += '; T_ℓ on the class: zero roots %s, repeated factors %s' % (R.get('zero'), R.get('rep'))
            say(line + ('  ← ' + ', '.join(mark) if mark else '') + '  (%.0f s)' % R['sec'])
    say(); say('flagged: %s  (Boylan–Webb\'s exception up to 1297: 607)' % (flagged or 'none'))
    say('seconds: %.0f' % (time.time() - t0))
