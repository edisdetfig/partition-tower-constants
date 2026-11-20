# A dig, 12 Nov 2026: DISTINCT PARTS, THE CONSTANTS FROM FORMS (README "Distinct parts: the constants from forms").
# Q(n) = the number of partitions of n into distinct parts, Σ Q(n)qⁿ = ∏(1 + qⁿ) = ∏(1 − q^(2n))/∏(1 − qⁿ).
# The tower is H_b(n) = Q(ℓᵇn + δ_b) with 24δ_b ≡ −1 (mod ℓᵇ).
# The forms: the composite C = T_ℓ∘T_ℓ on M_k(Γ0(2), v) (forms_eta2). Its ordinary part, the unit root X of
# X² − (C − 2c)X + c² = 0 (c = (2|ℓ)ℓ^(k−1)), X's characteristic polynomial at k = (ℓ − 1)(1 + x), x = 0..XW, and
# Newton's series at x = −1 (weight 0).
# Tested on the towers, in two-level steps: the law χ(y) applied to H_b, H_(b+2), H_(b+4), …, from every level b ≥ 1 and
# every n ≥ 1 (and with a factor y for the dying part where the lane-ℓ polynomial has one). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/distinct_parts_check.py N ℓ … [--XW=12] [--keep] [--wheel=…]
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
import forms_eta2 as FE
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; NS, PS = nums[0], nums[1:]
XW = next((int(a[5:]) for a in sys.argv[1:] if a.startswith('--XW=')), 12); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'distinct_parts_check_N%d_%s.txt' % (NS, '_'.join(map(str, PS)))), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def euler(M, step, n):
    P = nmod_poly([], M); terms = [(0, 1)]; j = 1
    while step * j * (3 * j - 1) // 2 < n:
        terms.append((step * j * (3 * j - 1) // 2, (-1) ** j))
        if step * j * (3 * j + 1) // 2 < n: terms.append((step * j * (3 * j + 1) // 2, (-1) ** j))
        j += 1
    for e, c in sorted(terms, reverse=True): P[e] = c % M
    return P


def family_poly(L, N):
    """χ(y) at weight 0 (lowest first, monic) modulo ℓ^prec, the lane-ℓ polynomial at weight ℓ − 1, and the details."""
    series = []; dims = []
    for x in range(XW + 1):
        k = (L - 1) * (1 + x); m = L ** N
        C, ok = FE.composite(k, L, m)
        d = len(C)
        Cm = F.mat_mod(d, d, [v for row in C for v in row], m)
        if x == 0:
            lane = [int(v) for v in nmod_mat(d, d, [v % L for row in C for v in row], L).charpoly().coeffs()]
        Y, R, r = F.ordinary_part(Cm, L, N)
        if x == 0:
            say('   weight %d: dimension %d, the images in the target space: %s; the composite on lane ℓ (in y): %s; ordinary part %d' % (
                k, d, ok, lane, r))
        if r == 0: series.append([1]); dims.append(0); continue
        Ml = F.restrict(Cm, Y, R, r)
        c = FE.kron2(L) * pow(L, k - 1, m) % m
        Xm = Ml
        for _ in range(N + 2):                                      # X = C − 2c + c²·X⁻¹ on the ordinary part
            Xinv = F.inverse_mod(Xm, L, N)
            Xm = Ml - F.identity(r, m) * (2 * c % m) + Xinv * (c * c % m)
        series.append(F.charpoly_small(Xm, r, m)); dims.append(r)
    if len(set(dims)) != 1:
        say('   the ordinary part changes size along the weights %s; stopped' % dims); return None, None, lane
    r = dims[0]
    poly = []; vals = []
    for j in range(r):
        v, vv = F1.newton_at([series[x][j] for x in range(XW + 1)], L, N, -1)
        poly.append(v); vals.append(vv)
    poly.append(1)
    prec = min(N, min((min(vv[1:]) if len(vv) > 1 else N) for vv in vals) + XW // 2) if vals else N
    return poly, vals, lane


t0 = time.time()
for L in PS:
    t = time.time(); K = 1
    while L ** (K + 1) < 2 ** 30: K += 1
    M = L ** K
    E1, E2 = euler(M, 1, NS), euler(M, 2, NS)
    Qs = E2.mul_low(E1.inverse_series_trunc(NS), NS); del E1, E2
    lv = {}
    for b in range(0, 40):
        if L ** b >= NS: break
        d_b = (-pow(24, -1, L ** b)) % (L ** b) if b else 0
        lv[b] = [int(Qs[L ** b * n + d_b]) % M for n in range(12) if L ** b * n + d_b < NS]
    del Qs
    top = max(lv)
    say(); say('ℓ = %d, lane ℓ^%d: tower levels 0..%d from n < %d (%.0f s)' % (L, K, top, NS, time.time() - t))
    poly, vals, lane = family_poly(L, K + 2)
    laws = []
    zero = next(i for i, c in enumerate(lane) if c % L)                 # y^zero divides the lane polynomial
    ordl = [c % L for c in lane[zero:]]
    laws.append(('lane ℓ, the ordinary part χ(y)', ordl, L, 1))
    if zero: laws.append(('lane ℓ, y^%d·χ(y)' % zero, [0] * zero + ordl, L, 1))
    if poly is not None:
        say('   the families\' polynomial χ(y) at weight 0 (lowest first), modulo ℓ^%d: %s;  differences\' valuations %s' % (K, [c % M for c in poly], vals))
        laws.append(('every digit: χ(y) at weight 0', [c % M for c in poly], M, K))
        if zero: laws.append(('every digit: y^%d·χ(y)' % zero, [0] * zero + [c % M for c in poly], M, K))
    for name, Q, mod, KK in laws:
        D = len(Q) - 1; rows = []; checks = 0
        for b in sorted(lv):
            if b < 1 or b + 2 * D > top: continue
            nn = min(len(lv[b + 2 * i]) for i in range(D + 1))
            ds = [val(sum(Q[i] * lv[b + 2 * i][n] for i in range(D + 1)), L, KK) for n in range(1, nn)]
            if ds: rows.append('%d: %d' % (b, min(ds))); checks += len(ds)
        say('   %-34s %s: smallest depth over n ≥ 1 (of %d), by level:  %s   (%d checks)' % (
            name, [c % L for c in Q], KK, '  '.join(rows) or 'no level span', checks))
say(); say('seconds: %.0f' % (time.time() - t0))
