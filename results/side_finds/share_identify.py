# share_identify.py, 14 Nov 2026 (README "The share condition, tested on partition numbers": which families give the
# factors). share_test.py sorts the tower's families at level 3 by the depth of their share: units, then ℓ². This script
# names them. In M_3 (share_test.py's grid, modulo ℓ⁵) each invariant factor comes with a direction in the space of
# functions of Q (a row of V⁻¹ in the Smith form U·M·V = D). Modulo ℓ:
#   D_0 = the directions of the unit factors,  D_2 = the directions of the factors of valuation ≤ 2.
# PREDICTED (derived): D_0 is spanned by the loud notes' Hecke numbers, Q ↦ (3|Q)·a_Q(f) for the class forms f of weight
# ℓ − 3; and D_2 = D_0 + the first-step systems' numbers, Q ↦ (3|Q)·Q·a_Q(g) for the class forms g of weight ℓ − 5 (θ-twisted:
# a_Q ↦ Q·a_Q). Both spans are read from the Hecke matrices on the classes (forms6_flint), matrix entry by entry, modulo ℓ.
# COUNTED (ranks). Partition numbers exact (FLINT), workers at idle priority.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/share_identify.py LMIN LMAX [--procs=N] [--keep]
import os, sys, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import share_test as ST


def smith_dirs(A, L, m):
    """Smith form of A over ℤ/ℓ^m with the column transform tracked: returns [(valuation, direction mod ℓ)], smallest first."""
    P = L ** m; A = [[x % P for x in row] for row in A]; nr, nc = len(A), len(A[0])
    W = [[1 if i == j else 0 for j in range(nc)] for i in range(nc)]; out = []
    for k in range(min(nr, nc)):
        best = None
        for i in range(k, nr):
            for j in range(k, nc):
                if A[i][j]:
                    v = ST.vl(A[i][j], L, m)
                    if best is None or v < best[0]: best = (v, i, j)
        if best is None:
            out += [(m, [x % L for x in W[kk]]) for kk in range(k, min(nr, nc))]; break
        v, i, j = best
        A[k], A[i] = A[i], A[k]
        for row in A: row[k], row[j] = row[j], row[k]
        W[k], W[j] = W[j], W[k]
        u = (A[k][k] // L ** v) % P; ui = pow(u, -1, P)
        for i2 in range(k + 1, nr):
            if A[i2][k]:
                f = (A[i2][k] // L ** v) * ui % P
                A[i2] = [(a - f * b) % P for a, b in zip(A[i2], A[k])]
        for j2 in range(k + 1, nc):
            if A[k][j2]:
                f = (A[k][j2] // L ** v) * ui % P
                A[k][j2] = 0
                W[k] = [(a + f * b) % P for a, b in zip(W[k], W[j2])]
        out.append((v, [x % L for x in W[k]]))
    return out


def rank_mod(vecs, L):
    rows = [list(v) for v in vecs if any(x % L for x in v)]; r = 0
    if not rows: return 0
    nc = len(rows[0])
    for c in range(nc):
        piv = next((i for i in range(r, len(rows)) if rows[i][c] % L), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]; inv = pow(rows[r][c], -1, L)
        rows[r] = [x * inv % L for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] % L:
                f = rows[i][c]; rows[i] = [(a - f * b) % L for a, b in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows): break
    return r


def class_functions(L, w, qs, twist_by_q):
    """Hecke numbers of the class at weight w (signs (2|ℓ), (3|ℓ)) as functions on the primes qs, modulo ℓ: every matrix
    entry of (3|Q)·Q^e·T_Q on the class (e = 1 for the first-step twist). Returns (the class's dimension, the list of
    functions)."""
    import forms6_flint as F
    N = 4; m = L ** N
    s2, s3 = ST.jacobi(2, L), ST.jacobi(3, L)
    H = F.hecke_matrices(w, m, [('T', q) for q in qs] + [('U', 2), ('U', 3)])
    res = F.class_padic(H, w, L, N, s2, s3)
    if res is None or res[3] is None: return 0, []
    W, RW, r, Wn, Rn, d, Np = res
    mp = L ** Np; n = w + 1
    funcs = {}
    for q in qs:
        Tq = F.mat_mod(n, n, [int(v) % mp for v in H[('T', q)].entries()], mp)
        Tc = F.restrict(Tq, Wn, Rn, d)
        sgn = ST.jacobi(3, q) * (q if twist_by_q else 1)
        for a in range(d):
            for b in range(d):
                funcs.setdefault((a, b), []).append(sgn * int(Tc[a, b]) % L)
    return d, list(funcs.values())


def main():
    nums = [int(a) for a in sys.argv[1:] if a.isdigit()]
    opt = {a.split('=')[0]: a.split('=')[1] for a in sys.argv[1:] if a.startswith('--') and '=' in a}
    KEEP = '--keep' in sys.argv; procs = int(opt.get('--procs', 4))
    primes = [p for p in range(nums[0], nums[1] + 1) if ST.is_prime(p) and p >= 29]
    out = os.path.join(HERE, 'share_identify_%d-%d.txt' % (nums[0], nums[1]))
    if KEEP and os.path.exists(out): sys.exit('stopped before computing anything: %s exists' % out)
    fh = open(out, 'x', encoding='utf-8') if KEEP else None

    def say(*a):
        line = ' '.join(str(x) for x in a); print(line, flush=True)
        if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())

    say('WHICH FAMILIES GIVE THE FACTORS AT LEVEL 3 (share_test.py\'s grid S = R + K + 2, modulo ℓ⁵; spans modulo ℓ)')
    say('   D_0: directions of the unit factors; D_2: of the factors of valuation ≤ 2. Loud: (3|Q)·a_Q(f), weight ℓ − 3;')
    say('   first step: (3|Q)·Q·a_Q(g), weight ℓ − 5. Predicted: D_0 = loud, D_2 = loud + first step.')
    t0 = time.time(); good = 0; tot = 0
    with Pool(procs, initializer=ST.init_worker) as pool:
        for L in primes:
            t1 = time.time()
            R, K = (L + 12) // 24, L // 24; S = R + K + 2
            tq, qs = ST.grid(L, S); P = L ** ST.M5
            args = sorted({(L ** 3 * t * Q * Q + 1) // 24 for t in tq[3] for Q in qs} | {(L ** 3 * t + 1) // 24 for t in tq[3]}, reverse=True)
            val = dict(pool.imap_unordered(ST.part_mod, [(n, P) for n in args], chunksize=1))
            c = lambda M: val[(M + 1) // 24]
            rows = []
            for t in tq[3]:
                M0 = L ** 3 * t
                rows.append([(c(Q * Q * M0) + ST.jacobi(M0, Q) * ST.jacobi(-3 % Q, Q) * pow(Q * Q, -1, P) * c(M0)) % P for Q in qs])
            sd = smith_dirs(rows, L, ST.M5)
            D0 = [d for v, d in sd if v == 0]; D2 = [d for v, d in sd if v <= 2]
            dl, loud = class_functions(L, L - 3, qs, False)
            df, fs = class_functions(L, L - 5, qs, True)
            rl, rf = rank_mod(loud, L), rank_mod(fs, L)
            r0, r2 = rank_mod(D0, L), rank_mod(D2, L)
            r0l = rank_mod(D0 + loud, L); r2lf = rank_mod(D2 + loud + fs, L); rlf = rank_mod(loud + fs, L)
            ok = (r0 == rl == r0l == dl) and (r2 == rlf == r2lf == dl + df)
            good += ok; tot += 1
            say('ℓ = %3d: valuations %s; classes: loud %d forms (span %d), first step %d forms (span %d), together span %d; '
                'D_0 rank %d, D_0 + loud %d; D_2 rank %d, D_2 + loud + first step %d  ->  %s  (%.0f s)' % (
                    L, [v for v, d in sd], dl, rl, df, rf, rlf, r0, r0l, r2, r2lf,
                    'D_0 = loud and D_2 = loud + first step' if ok else 'NOT AS PREDICTED', time.time() - t1))
    say(); say('as predicted at %d of %d primes  (%.0f s)' % (good, tot, time.time() - t0))


if __name__ == '__main__':
    main()
