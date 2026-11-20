# A dig, 12 Nov 2026: THE SAME DETECTOR ON MODULAR FUNCTIONS NEVER LOOKED AT. boundary_powers.py runs the data-only search
# for a tower law (Σ cᵢ·H_(b+i)(n) ≡ 0 on lane ℓ, degree ≤ 6, unknown cᵢ solved over F_ℓ) on partitions into squares and cubes,
# which have no modular origin. Here the identical search runs on three functions that ARE modular but that no part of this
# study has touched, and no forms are computed for them:
#   distinct parts  ∏(1 + qⁿ) = η(2τ)/η(τ)·q^(−1/24)   (weight 0), tower Q(ℓᵇn + δ_b), 24δ_b ≡ −1 (mod ℓᵇ);
#   overpartitions  ∏(1 + qⁿ)/(1 − qⁿ) = η(2τ)/η(τ)²   (weight −1/2), tower p̄(ℓᵇn);
#   T_2B            (η(τ)/η(2τ))²⁴ + 24 = q⁻¹ + 276q − 2048q² + …, a Monster moonshine function, the Hauptmodul of Γ0(2)
#                   (weight 0, j's relative of level 2), tower c(ℓᵇn).
# PREDICTED before the run: a law of low degree on lane ℓ exists for every prime, as for every modular case so far, in
# contrast to squares and cubes. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/group_members.py N ℓ … [--keep]
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
fh = open(os.path.join(HERE, 'group_members_%s.txt' % tag), 'x', encoding='utf-8') if KEEP else None
out = {'N': N, 'series': {}}


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def euler(M, step, n):
    P = nmod_poly([], M); terms = [(0, 1)]; j = 1
    while step * j * (3 * j - 1) // 2 < n:
        terms.append((step * j * (3 * j - 1) // 2, (-1) ** j))
        if step * j * (3 * j + 1) // 2 < n: terms.append((step * j * (3 * j + 1) // 2, (-1) ** j))
        j += 1
    for e, c in sorted(terms, reverse=True): P[e] = c % M
    return P


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
for name in ('distinct parts', 'overpartitions', 'T_2B'):
    out['series'][name] = {}
    for L in PS:
        t = time.time(); K = 1
        while L ** (K + 1) < 2 ** 30: K += 1
        M = L ** K; NN = N + 1
        E1, E2 = euler(M, 1, NN), euler(M, 2, NN)
        if name == 'distinct parts':
            F = E2.mul_low(E1.inverse_series_trunc(NN), NN); c = lambda m: int(F[m])
            dl = lambda b: (-pow(24, -1, L ** b)) % L ** b
        elif name == 'overpartitions':
            F = E2.mul_low(E1.pow_trunc(2, NN).inverse_series_trunc(NN), NN); c = lambda m: int(F[m]); dl = lambda b: 0
        else:                                                         # q·T_2B = ∏(1 − qⁿ)²⁴/(1 − q²ⁿ)²⁴ + 24q
            F = E1.pow_trunc(24, NN).mul_low(E2.pow_trunc(24, NN).inverse_series_trunc(NN), NN)
            c = lambda m: (int(F[m + 1]) + (24 if m == 0 else 0)) % M; dl = lambda b: 0
        del E1, E2
        first = [c(m) if c(m) <= M // 2 else c(m) - M for m in range(8)]
        lv = {}
        for b in range(0, 40):
            if L ** b >= N: break
            d0 = dl(b); vals = [c(L ** b * n + d0) for n in range(12) if L ** b * n + d0 < N]
            if not vals: break
            lv[b] = vals
        del F
        out['series'][name][str(L)] = {'K': K, 'levels': {str(b): v for b, v in lv.items()}}
        top = max(lv); res = []
        for b0 in (1, 2):
            found = None; neq = 0; dd = 0
            for d in range(1, 7):
                rows, rhs = [], []
                for b in sorted(lv):
                    if b < b0 or b + d > top: continue
                    nn = min(len(lv[b + i]) for i in range(d + 1))
                    for n in range(nn):
                        rows.append([lv[b + i][n] % L for i in range(d)]); rhs.append((-lv[b + d][n]) % L)
                if len(rows) <= d: break
                ok, sol, rk = solve_mod(rows, rhs, L); neq = len(rows); dd = d
                if ok: found = (d, sol + [1], len(rows), rk == d); break
            res.append('from level %d: %s' % (b0, ('a law of degree %d, %s, from %d equations (unique: %s)' % found) if found
                                                   else 'no law of degree ≤ %d (the last search had %d equations)' % (dd, neq)))
        say('%s, ℓ = %d (%.0f s): first values %s; levels 0..%d;  %s' % (name, L, time.time() - t, first, top, ';  '.join(res)))
if KEEP:
    json.dump(out, open(os.path.join(HERE, 'group_members_%s.json' % tag), 'x', encoding='utf-8'))
say('seconds: %.0f' % (time.time() - t0))
