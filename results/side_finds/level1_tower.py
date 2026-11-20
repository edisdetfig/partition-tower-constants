# A dig, 12 Nov 2026: TOWERS OF TWO MORE LEVEL-1 FORMS, j AND 1/E₄, FROM EXACT VALUES. The tower of a form F = Σ c(n)qⁿ
# is F | U_ℓᵇ = Σ c(ℓᵇn) qⁿ, so H_b(n) = c(ℓᵇn) (n ≥ 0; H_b(0) = c(0) at every level).
#   j    = E₄³/Δ = q⁻¹ + 744 + 196884q + …   weight 0, a pole at the cusp (Lehner's congruences at 2, 3, 5, 7, 11);
#   1/E₄ = 1 − 240q + 55440q² − …            weight −4, holomorphic at the cusp but a pole at ρ = e^(2πi/3) (Hardy and
#                                             Ramanujan's coefficients);
#   1/E₆ = 1 + 504q + …                       weight −6, a pole at i (added 12 Nov 2026, after the 1/E₄ runs). E₆ comes from
#                                             Ramanujan's identity E₆ = E₂E₄ − 3θE₄, with E₂ = θΔ/Δ.
#   E₄/j = Δ/E₄²                              weight 4, a DOUBLE pole at ρ (Zhang's theorem: a_(npˡ) ≡ (−3|p)·p·a_(npˡ⁻¹)
#                                             mod p^(3l)); added 12 Nov as a check on a proven case;
#   1/E₄²                                     weight −8, a double pole at ρ; added 12 Nov.
#   1/(E₄³ − j₀Δ) = 1/(Δ·(j − j₀))            weight −12, ONE simple pole, at the point with j = j₀ (--j0=…, default 2:
#                                             no complex multiplication); added 12 Nov.
#   1/(Δ^d·H(j)) = 1/Σ cᵢE₄^(3i)Δ^(d−i)        weight −12d, simple poles at the roots of H (--jpoly=c₀,c₁,…,1, lowest first;
#                                             --name= for the file name). For H the class polynomial of discriminant −15:
#                                             two poles at two CM points that Frobenius may swap; added 12 Nov.
# E₄ is built from Euler's sparse series alone: E₄(τ) = η(τ)¹⁶/η(2τ)⁸ + 256·η(2τ)¹⁶/η(τ)⁸ (checked here against 1 + 240σ₃(n)
# for n < 2000). Everything on lane ℓ^K (ℓ^K < 2^BITS) to N terms with FLINT. Windows: N ≤ 1.6·10⁷ (the Windows build
# crashes on products beyond about 2·10⁷ terms); WSL with --wheel=<unpacked python-flint wheel> for more. No forms are used
# here. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/level1_tower.py j|invE4 N ℓ … [--bits=30] [--keep] [--wheel=…] [--out=…]
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
OUTD = next((a[6:] for a in sys.argv[1:] if a.startswith('--out=')), HERE)
FORM = sys.argv[1]; assert FORM in ('j', 'invE4', 'invE6', 'E4overj', 'invE4sq', 'poleJ0', 'poleH')
JPOLY = next(([int(c) for c in a[8:].split(',')] for a in sys.argv[1:] if a.startswith('--jpoly=')), None)   # H(j), lowest first, monic
NAME = next((a[7:] for a in sys.argv[1:] if a.startswith('--name=')), 'H')
J0 = next((int(a[5:]) for a in sys.argv[1:] if a.startswith('--j0=')), 2)
N = int(sys.argv[2]); PS = [int(a) for a in sys.argv[3:] if a.isdigit()]; KEEP = '--keep' in sys.argv
BITS = next((int(a[7:]) for a in sys.argv[1:] if a.startswith('--bits=')), 30)
tag = '%s_N%d_%s' % (FORM + ('%d' % J0 if FORM == 'poleJ0' else NAME if FORM == 'poleH' else ''), N, '_'.join(map(str, PS)))
fh = open(os.path.join(OUTD, 'level1_tower_%s.txt' % tag), 'x', encoding='utf-8') if KEEP else None
out = {'form': FORM, 'N': N, 'primes': {}, 'j0': J0, 'jpoly': JPOLY}


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


def e4(M, n):
    E1, E2 = euler(M, 1, n), euler(M, 2, n)
    A = E1.pow_trunc(16, n).mul_low(E2.pow_trunc(8, n).inverse_series_trunc(n), n)
    B = E2.pow_trunc(16, n).mul_low(E1.pow_trunc(8, n).inverse_series_trunc(n), n)
    return A + B.mul_low(nmod_poly([0, 256 % M], M), n)              # + 256·q·B


def den_H(M, n):
    """Σ cᵢ·E₄^(3i)·Δ^(d−i) to n terms modulo M (Δ = q·∏(1 − qᵐ)²⁴), for H(j) = Σ cᵢ jⁱ monic of degree d."""
    d = len(JPOLY) - 1; E3 = e4(M, n).pow_trunc(3, n); Dl = euler(M, 1, n).pow_trunc(24, n).mul_low(nmod_poly([0, 1], M), n)
    tot = nmod_poly([0], M)
    for i, c in enumerate(JPOLY):
        tot = tot + (E3.pow_trunc(i, n) if i else nmod_poly([1], M)).mul_low(Dl.pow_trunc(d - i, n) if d - i else nmod_poly([1], M), n) * (c % M)
    return tot


def e6(M, n):
    """E₆ = E₂·E₄ − 3·θE₄ (Ramanujan), E₂ = 1 + θD/D with Δ = q·D, D = ∏(1 − qⁿ)²⁴; θ = q·d/dq."""
    E = e4(M, n); D = euler(M, 1, n).pow_trunc(24, n)
    th = lambda f: f.derivative().mul_low(nmod_poly([0, 1], M), n)
    E2 = nmod_poly([1], M) + th(D).mul_low(D.inverse_series_trunc(n), n)
    return E2.mul_low(E, n) - th(E) * 3


t0 = time.time()
small = 2000; s3 = [0] * small                                       # the check of E₄ against 1 + 240σ₃(n)
for d in range(1, small):
    for m in range(d, small, d): s3[m] += d ** 3
chk = e4(2 ** 61 - 1, small)
ok = all(int(chk[n]) == (1 if n == 0 else 240 * s3[n]) % (2 ** 61 - 1) for n in range(small))
say('%s: c(n) on lane ℓ^K for n < %d (FLINT %s%s); E₄ from eta quotients equals 1 + 240σ₃(n) for n < %d: %s' % (
    {'j': 'j = E₄³/Δ', 'invE4': '1/E₄', 'invE6': '1/E₆', 'E4overj': 'E₄/j = Δ/E₄²', 'invE4sq': '1/E₄²', 'poleJ0': '1/(E₄³ − %dΔ)' % J0, 'poleH': '1/(Δ^d·H(j)), H = %s (lowest first)' % JPOLY}[FORM], N, flint.__version__, ', WSL' if WHEEL else '', small, ok))
if not ok: sys.exit()
BIG = 2 ** 61 - 1; Es = e4(BIG, 9)                                    # the first exact values, as a check against the known ones
Fs = (Es.pow_trunc(3, 9).mul_low(euler(BIG, 1, 9).pow_trunc(24, 9).inverse_series_trunc(9), 9) if FORM == 'j' else
      Es.inverse_series_trunc(9) if FORM == 'invE4' else
      Es.pow_trunc(2, 9).inverse_series_trunc(9) if FORM == 'invE4sq' else
      nmod_poly([0], BIG) + euler(BIG, 1, 9).pow_trunc(24, 9).mul_low(Es.pow_trunc(2, 9).inverse_series_trunc(9), 9) if FORM == 'E4overj' else
      (Es.pow_trunc(3, 9) - euler(BIG, 1, 9).pow_trunc(24, 9).mul_low(nmod_poly([0, J0 % BIG], BIG), 9)).inverse_series_trunc(9) if FORM == 'poleJ0' else
      den_H(BIG, 9).inverse_series_trunc(9) if FORM == 'poleH' else
      e6(BIG, 9).inverse_series_trunc(9))
sgn = lambda v: v - BIG if v > BIG // 2 else v
say('   the first values c(%d..): %s' % (-1 if FORM == 'j' else 0, [sgn(int(Fs[i])) for i in range(8)]))
for L in PS:
    t = time.time(); K = 1
    while L ** (K + 1) < 2 ** BITS: K += 1
    M = L ** K; NN = N + 1 if FORM == 'j' else N
    E = e4(M, NN) if FORM not in ('invE6', 'poleH') else (e6(M, NN) if FORM == 'invE6' else None)
    if FORM == 'E4overj':                                            # Δ/E₄² = q·∏(1−qⁿ)²⁴/E₄²: c(n) = coefficient n − 1
        F = euler(M, 1, NN).pow_trunc(24, NN).mul_low(E.pow_trunc(2, NN).inverse_series_trunc(NN), NN); sh = -1
    elif FORM == 'invE4sq':
        F = E.pow_trunc(2, NN).inverse_series_trunc(NN); sh = 0
    elif FORM == 'poleH':
        F = den_H(M, NN).inverse_series_trunc(NN); sh = 0
    elif FORM == 'poleJ0':                                           # 1/(E₄³ − j₀·q·∏(1 − qⁿ)²⁴)
        F = (E.pow_trunc(3, NN) - euler(M, 1, NN).pow_trunc(24, NN).mul_low(nmod_poly([0, J0 % M], M), NN)).inverse_series_trunc(NN); sh = 0
    elif FORM == 'j':                                                  # q·j = E₄³ / ∏(1 − qⁿ)²⁴; c(n) = coefficient n + 1
        F = E.pow_trunc(3, NN).mul_low(euler(M, 1, NN).pow_trunc(24, NN).inverse_series_trunc(NN), NN); sh = 1
    else:
        F = E.inverse_series_trunc(NN); sh = 0
    del E
    c = lambda n: (int(F[n + sh]) if n + sh >= 0 else 0)
    first = [c(n) for n in range(6)]
    lv = {}
    for b in range(0, 40):
        if L ** b >= N: break
        vals = [c(L ** b * n) for n in range(12) if L ** b * n < N]
        if not vals: break
        lv[b] = vals
    del F
    out['primes'][str(L)] = {'K': K, 'levels': {str(b): v for b, v in lv.items()}}
    say(); say('ℓ = %d, lane ℓ^%d (%.0f s): c(0..5) on that lane %s; levels %s' % (L, K, time.time() - t, first,
                                                                                   {b: len(v) for b, v in lv.items()}))
    def val(x):
        if x == 0: return K
        v = 0
        while x % L == 0: x //= L; v += 1
        return v
    for b, v in lv.items():
        say('   b = %d: ℓ-adic valuations %s;  on lane ℓ %s' % (b, [val(x) for x in v], [x % L for x in v]))
if KEEP:
    json.dump(out, open(os.path.join(OUTD, 'level1_tower_%s.json' % tag), 'x', encoding='utf-8'))
say(); say('seconds: %.0f' % (time.time() - t0))
