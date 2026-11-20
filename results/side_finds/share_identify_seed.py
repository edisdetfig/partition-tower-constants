# share_identify_seed.py, 14 Nov 2026 (README "Why depth one? Every family starts as a unit"): which system starts deeper?
# At a prime where share_spectrum.py finds a family one digit deeper than its slope predicts, this names it. Grid S × S at
# level b (shift as in share_spectrum.py), Smith form with directions (share_identify.py). For every seed Sym^r ⊗ det^m on
# the row with slope ≤ SMAX, each of its systems' Hecke numbers on the grid's primes: Q ↦ (3|Q)·Q^m·a_Q(f), f the forms of
# weight r + 2 at level 6 with the twisted class signs (rational systems one by one; others as their span). Printed: for each
# system, the smallest depth d whose directions D_d contain its numbers. A unit start predicts d = its depth from the slope.
# COUNTED. Partition numbers exact (FLINT), workers at idle priority.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/share_identify_seed.py L S SHIFT LEVEL SMAX
#          [--prec=M (work modulo ℓ^M, default 6)] [--procs=N] [--keep]
import os, sys, time, math
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import flint
import forms6_flint as F
import share_test as ST
import share_spectrum as SS
from share_identify import smith_dirs, rank_mod


def systems_of_seed(L, r, m, qs):
    LM = L - 1; j = (L - 5 - r - 2 * m) // LM
    s2, s3 = ST.jacobi(2, L) ** ((1 + j) % 2), ST.jacobi(3, L) ** ((1 + j) % 2)
    w = r + 2; N = 4
    H = F.hecke_matrices(w, L ** N, [('T', q) for q in qs] + [('U', 2), ('U', 3)])
    res = F.class_padic(H, w, L, N, s2, s3)
    if res is None or res[3] is None: return []
    W, RW, r_, Wn, Rn, d, Np = res
    mp = L ** Np; n = w + 1
    T = {q: F.reduce_mat(F.restrict(F.mat_mod(n, n, [int(v) % mp for v in H[('T', q)].entries()], mp), Wn, Rn, d), L) for q in qs}
    twist = {q: ST.jacobi(3, q) * pow(q, m, L) % L for q in qs}
    # split by T at the first prime where the char poly has distinct linear factors; else keep the span (entry functions)
    for q0 in qs:
        fac = T[q0].charpoly().factor()[1]
        if all(len(f.coeffs()) == 2 and e == 1 for f, e in fac) and len(fac) == d:
            out = []
            for f, e in fac:
                lam = (-int(f.coeffs()[0])) % L
                A = T[q0] - flint.nmod_mat(d, d, [lam if a == b else 0 for a in range(d) for b in range(d)], L)
                ker = A.nullspace()[0]
                v = [int(ker[i, 0]) for i in range(d)]
                if not any(v): v = [int(ker[i, 1]) for i in range(d)]
                k = next(i for i in range(d) if v[i])
                vec = []
                for q in qs:
                    Tv = [sum(int(T[q][i, jj]) * v[jj] for jj in range(d)) % L for i in range(d)]
                    vec.append(twist[q] * Tv[k] * pow(v[k], -1, L) % L)
                out.append(('rational, T_%d ≡ %d' % (q0, lam), [vec]))
            return out
    funcs = [[twist[q] * int(T[q][a, b]) % L for q in qs] for a in range(d) for b in range(d)]
    return [('span of %d forms' % d, funcs)]


def main():
    L, S, SHIFT, b, SMAX = (int(x) for x in sys.argv[1:6])
    opt = {a.split('=')[0]: a.split('=')[1] for a in sys.argv[1:] if a.startswith('--') and '=' in a}
    KEEP = '--keep' in sys.argv; procs = int(opt.get('--procs', 6)); M = int(opt.get('--prec', 6))
    out = os.path.join(HERE, 'share_identify_seed_%d_grid%d_s%d_level%d%s.txt' % (L, S, SHIFT, b, ('_prec%d' % M) if M != 6 else ''))
    if KEEP and os.path.exists(out): sys.exit('stopped before computing anything: %s exists' % out)
    fh = open(out, 'x', encoding='utf-8') if KEEP else None

    def say(*a):
        line = ' '.join(str(x) for x in a); print(line, flush=True)
        if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())

    t0 = time.time(); P = L ** M
    qs = []; q = 5; k = 0
    while len(qs) < S:
        if ST.is_prime(q) and q != L:
            if k >= SHIFT: qs.append(q)
            k += 1
        q += 1
    cls = (23 * pow(L ** b, -1, 24)) % 24; ts = []; t = cls; k = 0
    while len(ts) < S:
        if t % L and ST.squarefree(t):
            if k >= SHIFT: ts.append(t)
            k += 1
        t += 24
    args = sorted({(L ** b * t * Q * Q + 1) // 24 for t in ts for Q in qs} | {(L ** b * t + 1) // 24 for t in ts}, reverse=True)
    with Pool(procs, initializer=ST.init_worker) as pool:
        val = dict(pool.imap_unordered(ST.part_mod, [(n, P) for n in args], chunksize=2))
    c = lambda Mm: val[(Mm + 1) // 24]
    mat = [[(c(Q * Q * L ** b * t) + ST.jacobi(L ** b * t, Q) * ST.jacobi(-3 % Q, Q) * pow(Q * Q, -1, P) * c(L ** b * t)) % P for Q in qs] for t in ts]
    sd = smith_dirs(mat, L, M)
    say('ℓ = %d, level %d, grid %d × %d shifted by %d (Q = %s): depths %s' % (L, b, S, S, SHIFT, qs, [v for v, d in sd]))
    sp = SS.spectrum(L); seen = set()
    for s, r, m, kind in sp:
        if s is None or s > SMAX or (r, m) in seen: continue
        seen.add((r, m))
        pred = SS.depth(s, b)
        for label, vecs in systems_of_seed(L, r, m, qs):
            found = next((dd for dd in range(M + 1) if rank_mod([d for v, d in sd if v <= dd] + vecs, L) == rank_mod([d for v, d in sd if v <= dd], L)), None)
            say('   seed Sym^%d ⊗ det^%d (slope %s%s), %s: in the directions of depth ≤ %s; a unit start predicts %d%s' % (
                r, m, s, ', ' + kind if kind else '', label, found, pred, '' if found == pred else '   <-- DIFFERENT'))
    say('(%.0f s)' % (time.time() - t0))


if __name__ == '__main__':
    main()
