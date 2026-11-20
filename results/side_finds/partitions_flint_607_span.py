# A dig, 12 Nov 2026: WHY 607 LOSES ONE. At 607 the tower's class has all 25 forms, none non-ordinary, but T_607 on it
# has a repeated root on lane 607 (partitions_flint_607_599_601_607_613.txt, partitions_flint_exceptions_new_*.txt).
# The tower's levels act on the class, on lane ℓ, as T_ℓ does; so the space the levels fill is spanned by v, Tv, T²v, …
# for a starting vector v, and its dimension is at most the degree of T_ℓ's MINIMAL polynomial on the class (derived).
# PREDICTED before the run: at 607 the two forms that share the root are told apart by the other Hecke numbers (T₅, T₇,
# …), so T_607 acts on them as one number, the minimal polynomial has degree 24, not 25, and the levels fill at most 24
# directions: Boylan–Webb's rank 24 = ⌊(607 + 12)/24⌋ − 1. Controls 599, 601, 613: degree = the class's dimension.
# Method (forms6_flint.py, FLINT only): the class at weight ℓ − 3 (exactly, class_padic, when the class on lane ℓ is larger
# than ⌊(ℓ + 12)/24⌋), T_ℓ and T₅ on it modulo ℓ; minimal and characteristic polynomials; for each repeated factor f^e of
# the characteristic polynomial, the number of independent eigenvectors (nullity of f(T_ℓ) / deg f) and whether T₅ splits
# that piece. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_607_span.py ℓ … [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat
import forms6_flint as F
flint.ctx.threads = 4
PS = [int(a) for a in sys.argv[1:] if a.isdigit()] or [599, 601, 607, 613]; KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_flint_607_span_%s.txt' % '_'.join(map(str, PS))), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


t0 = time.time()
for L in PS:
    t = time.time(); k = L - 3; n = k + 1; dexp = (L + 12) // 24; N = 4
    H = F.hecke_matrices(k, L ** N, [('U', 2), ('U', 3), ('T', 5), ('T', L)])
    U2, U3, T5, TL = (F.reduce_mat(H[op], L) for op in (('U', 2), ('U', 3), ('T', 5), ('T', L)))
    c2 = leg(2, L) * pow(2, (k - 2) // 2, L) % L; c3 = leg(3, L) * pow(3, (k - 2) // 2, L) % L
    K1, r1 = F.generalized_kernel(U2 - F.identity(n, L) * c2); X1, R1 = F.normalized_basis(K1, r1, L)
    K2, r0 = F.generalized_kernel(F.restrict(U3, X1, R1, r1) - F.identity(r1, L) * c3)
    Xc, Rc = F.normalized_basis(X1 * nmod_mat(r1, r0, [int(K2[i, j]) for i in range(r1) for j in range(r0)], L), r0, L)
    how = 'on lane ℓ'
    if r0 != dexp:
        W, RW, r, Wn, Rn, d, Np = F.class_padic(H, k, L, N, leg(2, L), leg(3, L))
        Xc, Rc, r0 = F.reduce_mat(Wn, L), Rn, d; how = 'exactly (%d of the %d on lane ℓ)' % (d, r)
    Tl = F.restrict(TL, Xc, Rc, r0); T5c = F.restrict(T5, Xc, Rc, r0)
    chi = Tl.charpoly(); mu = Tl.minpoly(); fac = chi.factor()[1]
    say('ℓ = %d (weight %d): the class %s has %d forms (⌊(ℓ + 12)/24⌋ = %d); T_ℓ on it: characteristic polynomial of degree %d, '
        'MINIMAL polynomial of degree %d  →  the levels fill at most %d directions  (%.0f s)' % (
            L, k, how, r0, dexp, chi.degree(), mu.degree(), mu.degree(), time.time() - t))
    for f, e in fac:
        if e < 2: continue
        G = F.eval_poly(f, Tl); nul = r0 - G.rank()
        K, kd = G.nullspace()                                  # the piece where f(T_ℓ) = 0; T₅ on it
        Xk, Rk = F.normalized_basis(K, kd, L); T5k = F.restrict(T5c, Xk, Rk, kd)
        say('   repeated factor of degree %d, %d times: %d independent eigenvectors (%s); T₅ on them: characteristic polynomial %s, '
            'minimal of degree %d (%s)' % (f.degree(), e, nul // f.degree(), 'T_ℓ acts there as one number' if nul // f.degree() == e
                                           else 'a Jordan block', [int(c) for c in T5k.charpoly().coeffs()], T5k.minpoly().degree(),
                                           'T₅ tells them apart' if T5k.minpoly().degree() == kd else 'T₅ does not tell them apart'))
say(); say('seconds: %.0f' % (time.time() - t0))
