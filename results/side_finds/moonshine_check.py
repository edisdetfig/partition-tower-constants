# A dig, 12 Nov 2026: MORE MOONSHINE, T_2A AND T_3C (README "More moonshine"; predictions written there first).
#   T_2A = (η(τ)/η(2τ))²⁴ + 24 + 4096·(η(2τ)/η(τ))²⁴ = q⁻¹ + 4372q + 96256q² + …   tower c(ℓᵇn);
#   T_3C = E₄(3τ)/η(3τ)⁸ = Σ a(m)·q^(3m−1), a(m) the coefficients of E₄/E₁⁸ = 1 + 248q + 4124q² + …;
#          tower a(ℓᵇn + δ_b), 3δ_b ≡ 1 (mod ℓᵇ).
# Per prime, on lane ℓ:
#   (1) the data-only law (Σ cᵢH_(b+i) ≡ 0, the least degree ≤ 6, cᵢ solved over 𝔽_ℓ, from level 1 and level 2);
#   (2) the forms' laws, tested on every level and n.
#       - T_2A: the level-1 cusp forms times the level-2 newforms of one W₂ class, for each rule (w₂ = +1, −1, (2|ℓ), −(2|ℓ));
#         "dies" where that has no unit root.
#       - T_3C: the ordinary part of T_ℓ (ℓ ≡ 1 mod 3) or of T_ℓ∘T_ℓ (ℓ ≡ 2 mod 3) on η¹⁶·M_(ℓ−9) (forms_eta3).
# COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/moonshine_check.py N ℓ … [--keep] [--wheel=…]
import os, sys, time, json
WHEEL = next((a[8:] for a in sys.argv[1:] if a.startswith('--wheel=')), None)
if WHEEL: sys.path.insert(0, WHEEL)
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(1 if WHEEL else 0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception:
    try: os.nice(19)
    except Exception: pass
import flint
from flint import nmod_poly, nmod_mat
flint.ctx.threads = 8
import forms6_flint as F
import forms1_flint as F1
import forms2_flint as F2
import forms_eta3 as FE3
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; N, PS = nums[0], nums[1:]; KEEP = '--keep' in sys.argv
tag = 'N%d_%s' % (N, '_'.join(map(str, PS)))
fh = open(os.path.join(HERE, 'moonshine_check_%s.txt' % tag), 'x', encoding='utf-8') if KEEP else None
out = {'N': N, 'series': {}}


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


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


def data_law(lv, L, b0):
    top = max(lv); neq = 0; dd = 0
    if all(v % L == 0 for b in lv if b >= b0 for v in lv[b][1:]): return 'the tower is ≡ 0 from level %d (dies)' % b0
    for d in range(1, 7):
        rows, rhs = [], []
        for b in sorted(lv):
            if b < b0 or b + d > top: continue
            nn = min(len(lv[b + i]) for i in range(d + 1))
            for n in range(1, nn):
                rows.append([lv[b + i][n] % L for i in range(d)]); rhs.append((-lv[b + d][n]) % L)
        if len(rows) <= d: break
        ok, sol, rk = solve_mod(rows, rhs, L); neq = len(rows); dd = d
        if ok: return 'degree %d, %s (lowest first), from %d equations, unique: %s' % (d, sol + [1], len(rows), rk == d)
    return 'none of degree ≤ %d (%d equations)' % (dd, neq)


def test_law(lv, Q, L, step, b0=1):
    """Σ Q_i·H_(b+step·i)(n) ≡ 0 on lane ℓ, for all b ≥ b0 and n ≥ 1. Returns (holds, checks, first failing level)."""
    D = len(Q) - 1; top = max(lv); checks = 0
    for b in sorted(lv):
        if b < b0 or b + step * D > top: continue
        nn = min(len(lv[b + step * i]) for i in range(D + 1))
        for n in range(1, nn):
            checks += 1
            if sum(Q[i] * lv[b + step * i][n] for i in range(D + 1)) % L: return False, checks, b
    return True, checks, None


def unit_part(cp, L):
    """The polynomial (lowest first) with the zero roots removed, on lane ℓ, and how many there were."""
    cp = [c % L for c in cp]; z = next(i for i, c in enumerate(cp) if c)
    return cp[z:], z


def level1_cusp(k, L):
    if F1.dim_mk(k) == 0: return [1]
    T = F1.hecke_matrices1(k, L, [L])[L]; d = T.nrows()
    cp = nmod_poly([int(c) for c in T.charpoly().coeffs()], L)
    q, r = divmod(cp, nmod_poly([L - 1, 1], L))                    # the Eisenstein root 1 + ℓ^(k−1) ≡ 1
    assert r == 0
    return [int(c) for c in q.coeffs()]


def level2_new(k, L, s2):
    K, d = F2.class_basis_exact2(k, s2)
    if d == 0: return [1]
    n = k // 4 + 1
    H = F2.hecke_matrices2(k, L, [('T', L)])[('T', L)]
    X = nmod_mat(n, d, [int(K[i, j]) % L for i in range(n) for j in range(d)], L)
    if len(F.pivot_rows(X, L)) < d: return None
    Xn, Rn = F.normalized_basis(X, d, L)
    return [int(c) for c in F.restrict(H, Xn, Rn, d).charpoly().coeffs()]


t0 = time.time()
for L in PS:
    K = 1
    while L ** (K + 1) < 2 ** 30: K += 1
    M = L ** K; NN = N + 2
    E1, E2 = euler(M, 1, NN), euler(M, 2, NN)
    # T_2A: q·T_2A = (E₁/E₂)²⁴ + 24q + 4096q²·(E₂/E₁)²⁴; c(m) = the coefficient of q^(m+1)
    R = E1.pow_trunc(24, NN).mul_low(E2.pow_trunc(24, NN).inverse_series_trunc(NN), NN)
    S = R.inverse_series_trunc(NN)
    qT = R + nmod_poly([0, 24], M) + S.mul_low(nmod_poly([0, 0, 4096 % M], M), NN)
    del R, S
    c2a = lambda m: int(qT[m + 1]) % M
    first2a = [c2a(m) if c2a(m) <= M // 2 else c2a(m) - M for m in range(-1, 4)]
    lv2 = {}
    for b in range(0, 40):
        if L ** b >= N: break
        lv2[b] = [c2a(L ** b * n) for n in range(12) if L ** b * n < N]
    del qT
    # T_3C: E₄/E₁⁸, E₄ = E₁¹⁶/E₂⁸ + 256q·E₂¹⁶/E₁⁸
    E1i8 = E1.pow_trunc(8, NN).inverse_series_trunc(NN)
    E4 = E1.pow_trunc(16, NN).mul_low(E2.pow_trunc(8, NN).inverse_series_trunc(NN), NN) + \
        E2.pow_trunc(16, NN).mul_low(E1i8, NN).mul_low(nmod_poly([0, 256 % M], M), NN)
    A3 = E4.mul_low(E1i8, NN); del E4, E1i8, E1, E2
    first3c = [int(A3[m]) for m in range(5)]
    lv3 = {}
    for b in range(0, 40 if L != 3 else 0):                           # 3 is T_3C's own level: no tower there
        if L ** b >= N: break
        d0 = pow(3, -1, L ** b) if b else 0
        lv3[b] = [int(A3[L ** b * n + d0]) % M for n in range(12) if L ** b * n + d0 < N]
    del A3
    out['series'][str(L)] = {'K': K, 'T_2A': {str(b): v for b, v in lv2.items()}, 'T_3C': {str(b): v for b, v in lv3.items()}}
    say(); say('ℓ = %d, lane ℓ^%d, n < %d: T_2A c(−1..3) %s, levels 0..%d;  T_3C a(0..4) %s, levels 0..%s' % (
        L, K, N, first2a, max(lv2), first3c, max(lv3) if lv3 else '— (ℓ = 3)'))
    # T_2A
    k = L - 1
    say('   T_2A, the data alone: from level 1: %s;  from level 2: %s' % (data_law(lv2, L, 1), data_law(lv2, L, 2)))
    if L > 3:
        c1 = level1_cusp(k, L); news = {s: level2_new(k, L, s) for s in (1, -1)}
        say('   T_2A, the forms at weight %d on lane ℓ: level-1 cusp forms %s; newforms on Γ0(2) with w₂ = +1: %s, w₂ = −1: %s' % (
            k, c1, news[-1], news[1]))
        for hn, w in (('w₂ = +1', 1), ('w₂ = −1', -1), ('w₂ = (2|ℓ)', leg(2, L)), ('w₂ = −(2|ℓ)', -leg(2, L))):
            nw = news[-w]                                             # U₂ = s2·2^((k−2)/2) with s2 = −w₂
            if nw is None: say('      %s: the newform basis degenerates on lane ℓ' % hn); continue
            full = [int(c) for c in (nmod_poly(c1, L) * nmod_poly(nw, L)).coeffs()]
            up, z = unit_part(full, L)
            if len(up) == 1: res = test_law(lv2, [1], L, 1, 1); law = 'dies'
            else: res = test_law(lv2, up, L, 1, 1); law = str(up)
            res2 = test_law(lv2, up if len(up) > 1 else [1], L, 1, 2)
            say('      %-12s law %s%s: from level 1 %s (%d checks);  from level 2 %s (%d checks)' % (
                hn, law, ' (×x^%d dropped)' % z if z else '', 'HOLDS' if res[0] else 'fails at level %s' % res[2], res[1],
                'HOLDS' if res2[0] else 'fails at level %s' % res2[2], res2[1]))
    # T_3C
    if L == 3: continue
    say('   T_3C, the data alone: from level 1: %s;  from level 2: %s' % (data_law(lv3, L, 1), data_law(lv3, L, 2)))
    if L % 3:
        C, step, ok = FE3.step_operator(k, L, L)
        d = len(C)
        if d == 0:
            law = [1]; z = 0; cp = []
        else:
            cp = [int(c) for c in nmod_mat(d, d, [v % L for row in C for v in row], L).charpoly().coeffs()]
            law, z = unit_part(cp, L)
        res = test_law(lv3, law, L, step, 1); res2 = test_law(lv3, law, L, step, 2)
        say('   T_3C, the forms: η¹⁶·M_%d (dimension %d); %s on lane ℓ: %s; images in the space: %s;  law %s in %s: from level 1 %s (%d checks);  from level 2 %s (%d checks)' % (
            k - 8, d, 'T_ℓ' if step == 1 else 'T_ℓ∘T_ℓ', cp, ok, 'dies' if len(law) == 1 else law,
            'one-level steps' if step == 1 else 'two-level steps', 'HOLDS' if res[0] else 'fails at level %s' % res[2], res[1],
            'HOLDS' if res2[0] else 'fails at level %s' % res2[2], res2[1]))
if KEEP: json.dump(out, open(os.path.join(HERE, 'moonshine_check_%s.json' % tag), 'x', encoding='utf-8'))
say(); say('seconds: %.0f' % (time.time() - t0))
