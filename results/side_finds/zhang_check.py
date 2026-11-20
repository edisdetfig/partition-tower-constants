# A dig, 12 Nov 2026: ZHANG'S CONJECTURES FOR E₄/(j − j(C))^r, TESTED AND EXTENDED (README "Zhang's conjectures for
# E₄/(j − j(C))^r"; predicted there first). P. Zhang, arXiv 2510.23200v3: Theorem 1.3, Conjectures 1.4 and 1.5.
#   F^(r) = E₄·Δ^r/(E₄³ − j₀Δ)^r = E₄/(j − j₀)^r: weight 4, a pole of order r at the points where j = j₀ (default j₀ = 2: no
#   complex multiplication). Its coefficients a(n), n < N, on the lane p^K (p^K < 2^BITS) with FLINT. E₄ comes from Euler's
#   series as in level1_tower.py; the first 40 coefficients are checked here against exact integers.
#   a_p(C) comes from counting points on y² = x³ + Ax + B, A = 3j₀(1728 − j₀), B = 2j₀(1728 − j₀)²: a curve with j = j₀.
#   A twist changes only the sign of a_p; the laws below use only u² and a_p².
# The tower is T_l(n) = a(n·pˡ) for n prime to p, with a(non-integer) = 0. (n divisible by p repeats a higher level.)
# A law ∏(X − ρ) acts with X as one level up. Its depth at level l is the least p-adic valuation, over n, of the law's sum
# ending at level l.
#   ordinary p:      L_I = ∏_(i<I) (X − pⁱ·u^(2−2i)), u the unit root of X² − a_pX + p.
#                    L_3 = (X − u²)(X − p)(X − ū²) is Zhang's cubic. L_4, L_5, L_6 add our next roots p³u⁻⁴, p⁴u⁻⁶, p⁵u⁻⁸.
#                    Probes at slope 3: L_I·(X − p³) (E₄'s critical root), I = 3..6 (added after a first small run at
#                    N = 2·10⁵ showed L_4 alone falling back to 3 digits a level).
#   supersingular p: a_p = 0, so α² = β² = −p. Zhang's cubic (X + p)²(X − p), its pieces, and slope-3 probes after X² − p².
# COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/zhang_check.py N p … [--j0=2] [--bits=63] [--keep]
#          (WSL: python3 -I zhang_check.py N p … --wheel=<unpacked python-flint wheel>)
import os, sys, time, json
WHEEL = next((a[8:] for a in sys.argv[1:] if a.startswith('--wheel=')), None)
if WHEEL: sys.path.insert(0, WHEEL)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception:
    try: os.nice(19)
    except Exception: pass
import flint
from flint import nmod_poly
flint.ctx.threads = 8
HERE = os.path.dirname(os.path.abspath(__file__))
J0 = next((int(a[5:]) for a in sys.argv[1:] if a.startswith('--j0=')), 2)
BITS = next((int(a[7:]) for a in sys.argv[1:] if a.startswith('--bits=')), 63)
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; N, PS = nums[0], nums[1:]; KEEP = '--keep' in sys.argv
NMAX = 200                                                           # n ≤ NMAX at every level
tag = 'j0%d_N%d_%s' % (J0, N, '_'.join(map(str, PS)))
fh = open(os.path.join(HERE, 'zhang_check_%s.txt' % tag), 'x', encoding='utf-8') if KEEP else None
out = {'j0': J0, 'N': N, 'primes': {}}


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def euler(M, step, n):
    """∏(1 − q^(step·k)) to n terms: Euler's pentagonal series with exponents times step."""
    P = nmod_poly([], M); terms = [(0, 1)]; j = 1
    while step * j * (3 * j - 1) // 2 < n:
        terms.append((step * j * (3 * j - 1) // 2, (-1) ** j))
        if step * j * (3 * j + 1) // 2 < n: terms.append((step * j * (3 * j + 1) // 2, (-1) ** j))
        j += 1
    for e, c in sorted(terms, reverse=True): P[e] = c % M
    return P


def forms(M, n):
    """F^(1), F^(2), F^(3) to n terms on lane M, as a generator (one at a time, to keep memory low)."""
    E1 = euler(M, 1, n); E2 = euler(M, 2, n)
    A = E1.pow_trunc(16, n).mul_low(E2.pow_trunc(8, n).inverse_series_trunc(n), n)
    B = E2.pow_trunc(16, n).mul_low(E1.pow_trunc(8, n).inverse_series_trunc(n), n); del E2
    E = A + B.mul_low(nmod_poly([0, 256 % M], M), n); del A, B          # E₄ = η¹⁶(τ)/η⁸(2τ) + 256·η¹⁶(2τ)/η⁸(τ)
    D = E1.pow_trunc(24, n).mul_low(nmod_poly([0, 1], M), n); del E1    # Δ = q·∏(1 − qᵐ)²⁴
    G = D.mul_low((E.pow_trunc(3, n) - D * (J0 % M)).inverse_series_trunc(n), n); del D   # 1/(j − j₀) = Δ/(E₄³ − j₀Δ)
    F = E
    for r in (1, 2, 3):
        F = F.mul_low(G, n)
        if r == 1: del E
        yield r, F


def exact_first(n):
    """F^(1..3) to n terms in exact integers: E₄ = 1 + 240σ₃, Δ from the product."""
    def mul(a, b): return [sum(a[i] * b[k - i] for i in range(k + 1)) for k in range(n)]
    def inv(a):
        b = [1] + [0] * (n - 1)
        for k in range(1, n): b[k] = -sum(a[i] * b[k - i] for i in range(1, k + 1))
        return b
    E = [1] + [240 * sum(d ** 3 for d in range(1, k + 1) if k % d == 0) for k in range(1, n)]
    P = [1] + [0] * (n - 1)
    for m in range(1, n):
        for _ in range(24): P = [P[k] - (P[k - m] if k >= m else 0) for k in range(n)]
    D = [0] + P[:n - 1]
    E3 = mul(mul(E, E), E); G = mul(D, inv([E3[k] - J0 * D[k] for k in range(n)]))
    res = {}; F = E
    for r in (1, 2, 3): F = mul(F, G); res[r] = F
    return res


def val(x, p, K):
    x %= p ** K
    if x == 0: return K
    v = 0
    while x % p == 0: x //= p; v += 1
    return v


def pmul(A, B, M):
    o = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        for j, b in enumerate(B): o[i + j] = (o[i + j] + a * b) % M
    return o


def from_roots(roots, M):
    Q = [1]
    for r in roots: Q = pmul(Q, [(-r) % M, 1], M)
    return Q


def a_p_curve(p):
    A = 3 * J0 * (1728 - J0); B = 2 * J0 * (1728 - J0) ** 2
    if (4 * A ** 3 + 27 * B ** 2) % p == 0: return None
    s = 0
    for x in range(p):
        y = (x ** 3 + A * x + B) % p
        if y: s += 1 if pow(y, (p - 1) // 2, p) == 1 else -1
    return -s


t0 = time.time()
EX = exact_first(40)
say('F^(r) = E₄/(j − %d)^r, r = 1, 2, 3, weight 4; a(n) for n < %d on lane p^K < 2^%d (FLINT %s%s)' % (J0, N, BITS, flint.__version__, ', WSL' if WHEEL else ''))
say('   exact a(1..8): r = 1: %s;  r = 2: %s;  r = 3: %s' % (EX[1][1:9], EX[2][1:9], EX[3][1:9]))
say('   Zhang: the cubic law reaches depth 3l − 2 − r at level l (Conj. 1.5); a(p) of F^(1) ≡ a_p(C)² mod p (Thm 1.3), ≡ 0 mod p² when p is supersingular (Conj. 1.4)')
for p in PS:
    t = time.time(); K = 1
    while p ** (K + 1) < 2 ** BITS: K += 1
    M = p ** K; tw = {}
    for r, F in forms(M, N):
        ok = all(int(F[k]) == EX[r][k] % M for k in range(40))
        if not ok: say('p = %d, r = %d: the FLINT series disagrees with the exact one: stopped' % (p, r)); sys.exit()
        lv = {}
        for m in range(0, 64):
            if p ** m >= N: break
            lv[m] = [int(F[k * p ** m]) for k in range(1, min(NMAX, (N - 1) // p ** m) + 1)]
        tw[r] = lv
    del F
    a = a_p_curve(p)
    out['primes'][str(p)] = {'K': K, 'a_p': a, 'towers': {str(r): {str(m): v for m, v in lv.items()} for r, lv in tw.items()}}
    say(); say('p = %d, lane p^%d (%.0f s; the first 40 terms equal the exact ones): a_p(C) = %s%s' % (
        p, K, time.time() - t, a, '' if a is None else ' (%s)' % ('supersingular' if a % p == 0 else 'ordinary')))
    if a is None: say('   bad reduction: skipped'); continue
    if a % p == 0:
        laws = [('X + p', [p], False), ('X − p', [-p], False), ('(X + p)(X − p)', [p, -p], False), ('(X + p)²', [p, p], False),
                ('Zhang: (X + p)²(X − p)', [p, p, -p], True), ('(X + p)²(X − p)²', [p, p, -p, -p], False), ('(X + p)³(X − p)²', [p, p, p, -p, -p], False),
                ('(X² − p²)(X − p³)', [p, -p, -p ** 3], False), ('(X² − p²)(X + p³)', [p, -p, p ** 3], False), ('(X² − p²)(X² − p⁶)', [p, -p, p ** 3, -p ** 3], False)]
        laws = [(nm, from_roots([(-x) % M for x in negs], M), z) for nm, negs, z in laws]   # stored as −roots
    else:
        u = a % p
        for _ in range(100): u = (u - (u * u - a * u + p) * pow((2 * u - a) % M, -1, M)) % M
        ro = [p ** i * pow(u, 2 - 2 * i, M) % M for i in range(6)]
        laws = [('L_%d' % I, from_roots(ro[:I], M), I == 3) for I in range(1, 7)] + \
               [('L_%d·(X − p³)' % I, from_roots(ro[:I] + [p ** 3], M), False) for I in (3, 4, 5, 6)]
        zc = [(-p ** 3) % M, p * (a * a - p) % M, (-(a * a - p)) % M, 1]
        say('   u = the unit root of X² − a_pX + p, mod p^%d: %d;  L_3 equals Zhang\'s cubic (X − α²)(X − p)(X − β²): %s' % (
            min(K, 6), u % p ** min(K, 6), laws[2][1] == zc))
    top = max(m for m in tw[1] if any(k % p for k in range(1, len(tw[1][m]) + 1)))
    if a % p: say('   (depths at levels l = 1..%d; n prime to p, n ≤ %d; "·" = no such n; depth %d = the lane\'s full width)' % (top, NMAX, K))
    else: say('   (depths at levels l = 1..%d; n prime to p, n ≤ %d; depth %d = the lane\'s full width)' % (top, NMAX, K))
    for r in (1, 2, 3):
        lv = tw[r]; T = lambda m, k: lv[m][k - 1] if m >= 0 else 0
        if r == 1:
            ap1 = lv[1][0]
            say('   r = 1: a(p) = %d mod p, a_p(C)² = %d mod p (Thm 1.3: %s); v_p(a(p)) = %d; least v_p(a(np)) over n prime to p: %d%s' % (
                ap1 % p, a * a % p, ap1 % p == a * a % p, val(ap1, p, K), min(val(T(1, k), p, K) for k in range(1, len(lv[1]) + 1) if k % p),
                ' (Conj 1.4 asks ≥ 2)' if a % p == 0 else ''))
        say('   r = %d, Zhang\'s bound 3l − 2 − r by level: %s' % (r, [max(0, 3 * l - 2 - r) for l in range(1, top + 1)]))
        for nm, Q, z in laws:
            D = len(Q) - 1; ds = []; fail = []
            for l in range(1, top + 1):
                ks = [k for k in range(1, len(lv[l]) + 1) if k % p]
                if not ks: ds.append('·'); continue
                vs = [val(sum(Q[s] * T(l - D + s, k) for s in range(D + 1)), p, K) for k in ks]
                ds.append(min(vs))
                if z: fail += [(l, k) for k, v in zip(ks, vs) if v < min(K, 3 * l - 2 - r)]
            say('      %-24s %s%s' % (nm + ':', ds, ('   ZHANG\'S LAW: holds at every level and n' if not fail else '   ZHANG\'S LAW FAILS at (l, n) = %s' % fail[:6]) if z else ''))
if KEEP: json.dump(out, open(os.path.join(HERE, 'zhang_check_%s.json' % tag), 'x', encoding='utf-8'))
say(); say('seconds: %.0f' % (time.time() - t0))
