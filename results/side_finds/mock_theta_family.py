# A dig, 12 Nov 2026: THE MOCK THETA FUNCTION f(q), THE CANDIDATE FAMILY (README "A mock theta function", "The unidentified
# piece"; predicted there first). The family through the level-6 newform of weight ℓ − 1 in the class with U₂, U₃ signs
# ((2|ℓ), (3|ℓ)), continued to the weights k = (ℓ − 1)(1 + x) with signs ((2|ℓ)^(x+1), (3|ℓ)^(x+1)) (the U₂, U₃ values must
# vary continuously, and 2^((ℓ−1)/2) ≡ (2|ℓ)).
# At each weight:
#   - the class exactly over the integers (class_basis_exact), saturated at ℓ;
#   - T_ℓ on it modulo ℓ^N, its ordinary part (one form, checked), and the unit root of x² − a·x + ℓ^(k−1);
#   - then Newton's series at x = −1, weight 0: u_f(0).
# Tested on the tower (mock_theta_tower json), in two-level steps y. By level: the smallest depth over n ≥ 1.
#   the shadow's polynomial S(y) = (y − (−3|ℓ))(y − χ₁₂(ℓ)ℓ)
#   S(y)·(y − (2|ℓ)u_f(0))                      PREDICTED
#   S(y)·(y − u_f(0)), S(y)·(y + u_f(0))        the other signs
#   (y − (−3|ℓ))·(y − (2|ℓ)u_f(0))             without the shadow's slope-1 root
# COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/mock_theta_family.py <tower json> ℓ … [--X=8] [--keep]
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
from math import gcd
from functools import reduce
from flint import nmod_mat
import forms6_flint as F
import forms1_flint as F1
SRC = next(a for a in sys.argv[1:] if a.endswith('.json')); PS = [int(a) for a in sys.argv[1:] if a.isdigit()]
X = next((int(a[4:]) for a in sys.argv[1:] if a.startswith('--X=')), 8); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'mock_theta_family_%s.txt' % '_'.join(map(str, PS))), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def saturate(cols, L):
    """Integer columns spanning a lattice; returns columns spanning its saturation at ℓ (full rank on lane ℓ)."""
    cols = [[x // reduce(gcd, [abs(y) for y in c if y]) for x in c] for c in cols]
    while True:
        d = len(cols); n = len(cols[0])
        Xl = nmod_mat(n, d, [cols[j][i] % L for i in range(n) for j in range(d)], L)
        if Xl.rank() == d: return cols
        ker, nul = Xl.nullspace()                                  # a combination ≡ 0 (mod ℓ)
        coef = [int(ker[j, 0]) for j in range(d)]
        jj = next(j for j in range(d) if coef[j] % L)
        new = [sum(coef[j] * cols[j][i] for j in range(d)) for i in range(n)]
        assert all(v % L == 0 for v in new)
        cols[jj] = [v // L for v in new]
        g = reduce(gcd, [abs(y) for y in cols[jj] if y]); cols[jj] = [v // g for v in cols[jj]]


def pmul(A, B, M):
    out = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        for j, b in enumerate(B): out[i + j] = (out[i + j] + a * b) % M
    return out


data = json.load(open(os.path.join(HERE, SRC), encoding='utf-8'))
t0 = time.time()
for L in PS:
    N = 12; m = L ** N; vals = []; dims = []
    say(); say('ℓ = %d: the class ((2|ℓ)^(x+1), (3|ℓ)^(x+1)) at weights (ℓ − 1)(1 + x), x = 0..%d, modulo ℓ^%d' % (L, X, N))
    for x in range(X + 1):
        k = (L - 1) * (1 + x); s2, s3 = leg(2, L) ** (x + 1), leg(3, L) ** (x + 1); t = time.time()
        K, d = F.class_basis_exact(k, s2, s3)
        if d == 0: say('   weight %d: the class is empty' % k); dims.append(0); continue
        cols = saturate([[int(K[i, j]) for i in range(k + 1)] for j in range(d)], L)
        Xm = nmod_mat(k + 1, d, [cols[j][i] % m for i in range(k + 1) for j in range(d)], m)
        Rn = F.pivot_rows(F.reduce_mat(Xm, L), L)
        H = F.hecke_matrices(k, m, [('T', L)])[('T', L)]
        TX = H * Xm
        XR = nmod_mat(d, d, [int(Xm[Rn[a], b]) for a in range(d) for b in range(d)], m)
        A = F.inverse_mod(XR, L, N) * nmod_mat(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m)
        Y, R2, r = F.ordinary_part(A, L, N)
        if r == 0: say('   weight %d: %d newforms in the class, none ordinary' % (k, d)); dims.append(0); continue
        Ml = F.restrict(A, Y, R2, r); U = F.unit_root_matrix(Ml, r, L, N, k); cp = F.charpoly_small(U, r, m)
        dims.append(r); vals.append(cp)
        say('   weight %d (signs %+d, %+d): %d newforms, %d ordinary; U_ℓ polynomial on lane ℓ %s  (%.0f s)' % (
            k, s2, s3, d, r, [c % L for c in cp], time.time() - t))
    if len(set(dims)) != 1 or dims[0] != 1:
        say('   the ordinary part is not one form at every weight (%s); stopped' % dims); continue
    u, vv = F1.newton_at([-vals[x][0] % m for x in range(X + 1)], L, N, -1)
    prec = min(N, X + 1); M = L ** prec; u %= M
    say('   u_f(0) (weight 0), digits lowest first: %s;  the differences\' valuations %s' % ([(u // L ** e) % L for e in range(prec)], vv))
    rec = data['primes'][str(L)]; K = min(rec['K'], prec); MK = L ** K
    lv = {int(b): [v % MK for v in vs] for b, vs in rec['levels'].items()}; top = max(lv)
    m3, c12, s2l = leg(-3, L), leg(3, L), leg(2, L)
    Sh = pmul([(-m3) % MK, 1], [(-c12 * L) % MK, 1], MK)
    laws = [('the shadow\'s polynomial S(y)', Sh),
            ('S(y)·(y − (2|ℓ)u_f(0))  PREDICTED', pmul(Sh, [(-s2l * u) % MK, 1], MK)),
            ('S(y)·(y − u_f(0))', pmul(Sh, [(-u) % MK, 1], MK)),
            ('S(y)·(y + u_f(0))', pmul(Sh, [u % MK, 1], MK)),
            ('(y − (−3|ℓ))·(y − (2|ℓ)u_f(0))', pmul([(-m3) % MK, 1], [(-s2l * u) % MK, 1], MK))]
    for name, Q in laws:
        D = len(Q) - 1; rows = []; checks = 0
        for b in sorted(lv):
            if b < 1 or b + 2 * D > top: continue
            nn = min(len(lv[b + 2 * i]) for i in range(D + 1))
            ds = [val(sum(Q[i] * lv[b + 2 * i][n] for i in range(D + 1)), L, K) for n in range(1, nn)]
            if ds: rows.append('%d: %d' % (b, min(ds))); checks += len(ds)
        say('   %-40s depth (of %d) by level  %s   (%d checks)' % (name, K, '  '.join(rows) or 'no level span', checks))
say(); say('seconds: %.0f' % (time.time() - t0))
