# A dig, 12 Nov 2026: r-COLORED PARTITIONS ON LEVEL 6, THE CONSTANTS FROM MODULAR FORMS (FLINT only; forms6_flint.py).
# For r prime to 3 the class is on Γ0(6) (colors5_census.txt for r = 5): the newforms of weight k₀ = ℓ − r − 2 with
# U₂ = ε₂·(2|ℓ)·2^((k−2)/2) and U₃ = ε₃·(3|ℓ)·3^((k−2)/2) (r = 1: ε₂ = ε₃ = +1, partitions_flint_family.py; r = 5: −1, −1).
# 1  the class at k₀, exactly over the integers; its Hecke numbers on lane ℓ;
# 2  at k = k₀ + (ℓ − 1)·x: the ordinary part of T_ℓ, the piece with the class's numbers on lane ℓ lifted to ℓ^N, the forms
#    in it with U₂, U₃ EXACTLY the class's values (signs carried by (2|ℓ)^x, (3|ℓ)^x), the unit root of T_ℓ there;
# 3  Newton's series at x = −1, i.e. weight −(r + 1): the polynomial of u(−(r+1)); for one form its digits (and of −u).
# PREDICTED before the run: the differences divisible by ℓⁿ; the constants (with the twist the first digits show) agree with
# every digit the exact r-colored partition numbers settle (colors_tower_r*.json). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors6_family.py ℓ X N --r=5 --e2=-1 --e3=-1 [--keep] [--threads=t]
import os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat
import forms6_flint as F
flint.ctx.threads = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--threads=')), 8)
opt = lambda name, dflt: next((int(a[len(name) + 3:]) for a in sys.argv[1:] if a.startswith('--%s=' % name)), dflt)
args = [int(a) for a in sys.argv[1:] if a.isdigit()]; L, X, N = args[0], args[1], args[2]
RC, E2, E3 = opt('r', 5), opt('e2', -1), opt('e3', -1); KEEP = '--keep' in sys.argv
OUT = os.path.join(HERE, 'colors6_family_r%d_%d.txt' % (RC, L))
if KEEP and os.path.exists(OUT): OUT = OUT[:-4] + '_N%d.txt' % N
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


t0 = time.time(); m = L ** N; k0 = L - RC - 2; s2b, s3b = E2 * leg(2, L), E3 * leg(3, L)
QS = [q for q in (5, 7, 11, 13) if q != L]; OPS = [('T', q) for q in QS] + [('T', L)]; OPSF = OPS + [('U', 2), ('U', 3)]
say('ℓ = %d, %d colours, k₀ = %d, level 6, U₂ sign %+d, U₃ sign %+d at k₀, %d weights, modulo ℓ^%d; FLINT %s' % (
    L, RC, k0, s2b, s3b, X + 1, N, flint.__version__))

# ---- 1. the class at the base weight, exactly ----
dz = F.class_basis_exact(k0, s2b, s3b)[1]                  # how many forms, exactly
Hb = F.hecke_matrices(k0, L ** 4, OPS + [('U', 2), ('U', 3)])
W_, RW_, r_, Wn, Rn, d, Np0 = F.class_padic(Hb, k0, L, 4, s2b, s3b)   # the class on lane ℓ, separated exactly
if d == 0 or d != dz or Wn is None: say('   the class: %d forms exactly, %d on lane ℓ; stopped' % (dz, d)); sys.exit()
Xc, Rc = F.reduce_mat(Wn, L), Rn
ref = {q: F.restrict(F.reduce_mat(Hb[('T', q)], L), Xc, Rc, d) for q in QS + [L]}
say('1  the class at weight %d: %d form%s;  T_ℓ polynomial on lane ℓ %s' % (k0, d, '' if d == 1 else 's', [int(c) for c in ref[L].charpoly().coeffs()]))

# ---- 2. the family ----
rng = random.Random(L); series = []; precs = []; coefs = None


def new_part(W, Wr, rt, k, x, S):
    c2 = s2b * leg(2, L) ** x * 2 ** ((k - 2) // 2); c3 = s3b * leg(3, L) ** x * 3 ** ((k - 2) // 2)
    U2W = F.restrict(S[('U', 2)], W, Wr, rt); U3W = F.restrict(S[('U', 3)], W, Wr, rt)
    Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(rt)] for i in range(rt)] + \
         [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(rt)] for i in range(rt)]
    ker, Np = F.padic_kernel(Bm, 2 * rt, rt, L, N)
    if len(ker) != d: return None, len(ker)
    mp_ = L ** Np
    Xn = F.mat_mod(rt, d, [ker[c][i] % mp_ for i in range(rt) for c in range(d)], mp_)
    Rn = F.pivot_rows(F.reduce_mat(Xn, L), L)
    XnRinv = F.inverse_mod(F.mat_mod(d, d, [int(Xn[Rn[a], b]) for a in range(d) for b in range(d)], mp_), L, Np)
    TW = F.restrict(S[('T', L)], W, Wr, rt); TWp = F.mat_mod(rt, rt, [int(v) % mp_ for v in TW.entries()], mp_)
    TX = TWp * Xn
    return XnRinv * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], mp_), Np


for x in range(X + 1):
    k = k0 + (L - 1) * x; t = time.time()
    H = F.hecke_matrices(k, m, OPSF)
    Y, Ry, r = F.ordinary_part(H[('T', L)], L, N)
    S = {op: F.restrict(H[op], Y, Ry, r) for op in OPSF}
    out, info = None, None
    for attempt in range(12):
        if coefs is None or attempt > 0: coefs = {op: rng.randrange(1, L) for op in OPS}
        Ar = sum((ref[p] * coefs[('T', p)] for p in QS + [L]), nmod_mat(d, d, [0] * (d * d), L))
        g = Ar.charpoly(); sq = g.factor_squarefree()[1]
        if len(sq) == 1 and sq[0][1] == 1:
            A = None
            for op in OPS: A = S[op] * coefs[op] if A is None else A + S[op] * coefs[op]
            out, info = F.lift_invariant(A, g, L, N)
            if out is not None: break
    if out is None: say('   weight %d: the piece was not found (%s); stopped' % (k, info)); sys.exit()
    Xk, Mk, Rk = out; rt = Xk.ncols(); Np = N; note = ''
    if rt == d: Ml = F.restrict(S[('T', L)], Xk, Rk, d)
    else:
        Ml, Np = new_part(Xk, Rk, rt, k, x, S)
        if Ml is None: say('   weight %d: the piece has %d forms and its new part %d, not %d; stopped' % (k, rt, Np, d)); sys.exit()
        note = ' (the piece had %d forms; kept the part new at 2 and 3, to ℓ^%d)' % (rt, Np)
    U = F.unit_root_matrix(Ml, d, L, Np, k); cp = F.charpoly_small(U, d, L ** Np)
    series.append(cp); precs.append(Np)
    say('2  weight %d (dimension %d, ordinary part %d): the class\'s piece %d%s; U_ℓ polynomial on lane ℓ %s  (%.0f s)' % (
        k, k + 1, r, d, note, [c % L for c in cp], time.time() - t))

# ---- 3. weight −(r + 1) ----
P = min(precs); mP = L ** P; poly = []; vals = []
for j in range(d):
    v, val, _ = F.newton_at_minus_one([series[x][j] % mP for x in range(X + 1)], L, P)
    poly.append(v); vals.append(val)
poly.append(1)
say('3  the differences of each coefficient are divisible by ℓ to %s' % vals)
K = min(P, X + 1)
say('   the polynomial of u(−%d) (lowest first) on lane ℓ: %s;  on lane ℓ^%d: %s' % (RC + 1, [c % L for c in poly], K, [c % L ** K for c in poly]))
dig = lambda v: [(v // L ** e) % L for e in range(K)]
if d == 1:
    u = (-poly[0]) % mP
    say('   u(−%d), digits lowest first: %s;  −u: %s' % (RC + 1, dig(u), dig((-u) % mP)))
say(); say('seconds: %.0f' % (time.time() - t0))
