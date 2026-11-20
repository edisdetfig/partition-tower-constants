# A dig, 12 Nov 2026: WHERE THE RECIPE STOPS. Every tower so far came from a modular object (η-quotients, spt, j, 1/E₄), and
# in every one the data alone found a short law on lane ℓ: Σ cᵢ·H_(b+i)(n) ≡ 0 with a few fixed constants (level1_check.py,
# colors_check, …). The reason (Serre, Hida): U_ℓ acts on ℓ-adic modular forms, and on lane ℓ the iterates of a modular
# form stay in a space of bounded size, so a law of bounded degree must hold.
# Here: two generating functions with NO modular origin, partitions into squares ∏(1 − q^(k²))⁻¹ and into cubes
# ∏(1 − q^(k³))⁻¹. Tower H_b(n) = c(ℓᵇn), n = 0…11.
# PREDICTED before the run: no law of degree ≤ 6 on lane ℓ, for any ℓ, from any starting level (the same search that found
# the law in every modular case). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/boundary_powers.py N ℓ … [--keep]
import os, sys, time, json
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly
flint.ctx.threads = 8
HERE = os.path.dirname(os.path.abspath(__file__))
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; N, PS = nums[0], nums[1:]; KEEP = '--keep' in sys.argv
tag = 'N%d_%s' % (N, '_'.join(map(str, PS)))
fh = open(os.path.join(HERE, 'boundary_powers_%s.txt' % tag), 'x', encoding='utf-8') if KEEP else None
out = {'N': N, 'series': {}}


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def product_tree(polys, M, n):
    while len(polys) > 1:
        nxt = []
        for i in range(0, len(polys) - 1, 2): nxt.append(polys[i].mul_low(polys[i + 1], n))
        if len(polys) % 2: nxt.append(polys[-1])
        polys = nxt
    return polys[0]


def solve_mod(rows, rhs, p):
    m = len(rows[0]) if rows else 0
    A = [list(r) + [b] for r, b in zip(rows, rhs)]; piv = []; r = 0
    for c in range(m):
        k = next((i for i in range(r, len(A)) if A[i][c] % p), None)
        if k is None: continue
        A[r], A[k] = A[k], A[r]; inv = pow(A[r][c], -1, p)
        A[r] = [x * inv % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] % p:
                f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if any(all(x % p == 0 for x in row[:m]) and row[m] % p for row in A): return False, None, r
    sol = [0] * m
    for i, c in enumerate(piv): sol[c] = A[i][m]
    return True, sol, r


t0 = time.time()
for name, e in (('squares', 2), ('cubes', 3)):
    out['series'][name] = {}
    for L in PS:
        t = time.time(); K = 1
        while L ** (K + 1) < 2 ** 30: K += 1
        M = L ** K; polys = []; k = 1
        while k ** e < N:
            P = nmod_poly([], M); P[k ** e] = M - 1; P[0] = 1; polys.append(P); k += 1
        F = product_tree(polys, M, N).inverse_series_trunc(N)
        first = [int(F[i]) for i in range(12)]
        lv = {}
        for b in range(0, 40):
            if L ** b >= N: break
            lv[b] = [int(F[L ** b * n]) for n in range(12) if L ** b * n < N]
        del F
        out['series'][name][str(L)] = {'K': K, 'levels': {str(b): v for b, v in lv.items()}}
        top = max(lv)
        res = []
        for b0 in (1, 2):
            found = None; neq = 0
            for d in range(1, 7):
                rows, rhs = [], []
                for b in sorted(lv):
                    if b < b0 or b + d > top: continue
                    nn = min(len(lv[b + i]) for i in range(d + 1))
                    for n in range(nn):
                        rows.append([lv[b + i][n] % L for i in range(d)]); rhs.append((-lv[b + d][n]) % L)
                if len(rows) <= d: break
                ok, sol, rk = solve_mod(rows, rhs, L); neq = len(rows)
                if ok: found = (d, sol + [1], len(rows), rk == d); break
            res.append('from level %d: %s' % (b0, ('a law of degree %d, %s, from %d equations (unique: %s)' % found) if found
                                                   else 'no law of degree ≤ %d (the last search had %d equations)' % (d, neq)))
        say('%s, ℓ = %d (%.0f s): c(0..11) = %s; levels 0..%d;  %s' % (name, L, time.time() - t, first, top, ';  '.join(res)))
if KEEP:
    json.dump(out, open(os.path.join(HERE, 'boundary_powers_%s.json' % tag), 'x', encoding='utf-8'))
say('seconds: %.0f' % (time.time() - t0))
