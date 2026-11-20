# A dig, 12 Nov 2026: A MOCK THETA FUNCTION'S TOWERS, THE DATA ALONE (NEXT.md item 11). Ramanujan's third-order
#   f(q) = Σ_(n≥0) q^(n²)/(−q;q)_n² = 1 + q − 2q² + 3q³ − 3q⁴ + 3q⁵ − 5q⁶ + …   (α(n) its coefficients)
# is the holomorphic part of a weight-1/2 harmonic Maass form (q^(−1/24)·f), not a modular form. Its tower is
# α(ℓᵇn + δ_b) with 24δ_b ≡ 1 (mod ℓᵇ), the index of p(n).
# Computed by the Appell–Lerch identity f(q) = (1 + 4·Σ_(n≥1) (−1)ⁿ q^(n(3n+1)/2)/(1 + qⁿ))/(q;q)_∞. The identity is checked
# here against the definition for n < 400 (exact integers).
# Per prime, on lane ℓ: the least law Σ cᵢH_(b+i) ≡ 0 (degree ≤ 6, cᵢ over 𝔽_ℓ) from level 1 and from level 2, and whether
# the tower dies. No forms; no prediction beyond the project's standing one (a modular object has a law).
# COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/mock_theta_tower.py N ℓ … [--keep]
#          (WSL, n < 10⁸: python3 -I mock_theta_tower.py N ℓ … --appell=<mock_theta_appell output> --wheel=<flint wheel>)
import os, sys, time, json
WHEEL = next((a[8:] for a in sys.argv[1:] if a.startswith('--wheel=')), None)
if WHEEL: sys.path.insert(0, WHEEL)
APPELL = next((a[9:] for a in sys.argv[1:] if a.startswith('--appell=')), None)   # exact A(q) from mock_theta_appell.py
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(1 if WHEEL else 0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception:
    try: os.nice(19)
    except Exception: pass
try:
    import numpy as np
except ImportError:
    np = None                                                        # WSL: A(q) comes from --appell
import flint
from flint import nmod_poly
flint.ctx.threads = 8
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; N, PS = nums[0], nums[1:]; KEEP = '--keep' in sys.argv
tag = 'N%d_%s' % (N, '_'.join(map(str, PS)))
fh = open(os.path.join(HERE, 'mock_theta_tower_%s.txt' % tag), 'x', encoding='utf-8') if KEEP else None
out = {'N': N, 'primes': {}}


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


def appell_py(n_):
    """The same, in pure Python (small n_ only: the identity check)."""
    a = [0] * n_; a[0] = 1; n = 1
    while n * (3 * n + 1) // 2 < n_:
        s = n * (3 * n + 1) // 2
        for i, e in enumerate(range(s, n_, n)): a[e] += 4 * (-1) ** n * (1 if i % 2 == 0 else -1)
        n += 1
    return a


def appell(n_):
    """1 + 4·Σ_(n≥1) (−1)ⁿ q^(n(3n+1)/2)/(1 + qⁿ) to n_ terms, as int64."""
    a = np.zeros(n_, dtype=np.int64); a[0] = 1; n = 1
    while n * (3 * n + 1) // 2 < n_:
        s = n * (3 * n + 1) // 2; cnt = len(range(s, n_, n))
        pat = np.empty(cnt, dtype=np.int64); pat[0::2] = 1; pat[1::2] = -1
        a[s::n] += 4 * (-1) ** n * pat; n += 1
    return a


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


# the identity, checked exactly for n < 400
S = 400; defn = [0] * S; den = [1] + [0] * (S - 1)                      # Σ q^(n²)/(−q;q)_n², exact
n = 0
while n * n < S:
    if n:
        for _ in range(2):                                              # divide by (1 + qⁿ)
            for i in range(n, S): den[i] -= den[i - n]
    for i in range(S - n * n): defn[n * n + i] += den[i]
    n += 1
ap = appell_py(S); e1 = [0] * S
inv = [0] * S; inv[0] = 1                                               # 1/(q;q)_∞ = Σ p(n)qⁿ, exact
for i in range(1, S):
    s_ = 0; j = 1
    while j * (3 * j - 1) // 2 <= i:
        sg = 1 if j % 2 else -1
        s_ += sg * inv[i - j * (3 * j - 1) // 2]
        if j * (3 * j + 1) // 2 <= i: s_ += sg * inv[i - j * (3 * j + 1) // 2]
        j += 1
    inv[i] = s_
viaid = [sum(ap[j] * inv[i - j] for j in range(i + 1)) for i in range(S)]
say('f(q) by its definition: %s …; the Appell–Lerch identity agrees for n < %d: %s' % (defn[:11], S, defn == viaid))
if defn != viaid: sys.exit()
t0 = time.time()
for L in PS:
    t = time.time(); K = 1
    while L ** (K + 1) < 2 ** 30: K += 1
    M = L ** K
    if APPELL:
        from array import array
        raw = array('i'); raw.fromfile(open(APPELL, 'rb'), N)
        A = nmod_poly([x % M for x in raw], M); del raw
    else:
        a = appell(N) % M
        A = nmod_poly([int(x) for x in a], M); del a
    F = A.mul_low(euler(M, 1, N).inverse_series_trunc(N), N); del A
    lv = {}
    for b in range(0, 40):
        if L ** b >= N: break
        d0 = pow(24, -1, L ** b) if b else 0
        lv[b] = [int(F[L ** b * n + d0]) % M for n in range(12) if L ** b * n + d0 < N]
    del F
    out['primes'][str(L)] = {'K': K, 'levels': {str(b): v for b, v in lv.items()}}
    say(); say('ℓ = %d, lane ℓ^%d (%.0f s): levels 0..%d; (3|ℓ) = %+d, (−1|ℓ) = %+d' % (
        L, K, time.time() - t, max(lv), 1 if pow(3, (L - 1) // 2, L) == 1 else -1, 1 if L % 4 == 1 else -1))
    say('   the data alone: from level 1: %s;  from level 2: %s' % (data_law(lv, L, 1), data_law(lv, L, 2)))
if KEEP: json.dump(out, open(os.path.join(HERE, 'mock_theta_tower_%s.json' % tag), 'x', encoding='utf-8'))
say(); say('seconds: %.0f' % (time.time() - t0))
