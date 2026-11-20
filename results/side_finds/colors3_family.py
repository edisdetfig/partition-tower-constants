# A dig, 12 Nov 2026: 3-COLORED PARTITIONS, THE CONSTANTS FROM MODULAR FORMS (FLINT only; forms2_flint.py, level 2).
# The recipe, carried over (colors3_census.txt: the class is on Γ0(2), weight ℓ − 5, U₂ = (2|ℓ)·2^((ℓ−7)/2), exactly
# R_ℓ(3) forms): each constant is t(ℓ)·u(−4), u the U_ℓ-number of the ordinary (Hida) family through a class form,
# continued from the weights k = ℓ − 5 + (ℓ − 1)·x to x = −1, weight −4 (the Shimura weight of η^(−3)); t(ℓ) = ±1 is the
# twist the first digits show (colors3_compare.txt).
# 1  the class at weight ℓ − 5, on lane ℓ (exactly separated from forms of level 1);
# 2  at each weight: the ordinary part of T_ℓ, the piece with the class's Hecke numbers on lane ℓ (T₃, T₅, T₇, …),
#    lifted to ℓ^N; the forms new at 2 with U₂ = s₂(k)·2^((k−2)/2) exactly, s₂ carried along: s₂(k) = (2|ℓ)^(1+x);
#    T_ℓ on them, the unit root (the ordinary ℓ-stabilization), its characteristic polynomial;
# 3  Newton's series of each coefficient at x = −1; the constants' polynomial; its roots (on lane ℓ, or A ± B·√r).
# PREDICTED before the run: the first digits are those of colors3_compare.txt, and the digits agree with every digit the
# exact 3-colored partition numbers settle (colors3_tower_*.json) once the twist t(ℓ) is applied. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors3_family.py ℓ X N [--keep] [--threads=t]
import os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly, nmod_mat
import forms2_flint as F2
import forms6_flint as F
flint.ctx.threads = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--threads=')), 8)
args = [int(a) for a in sys.argv[1:] if a.isdigit()]; L, X, N = args[0], args[1], args[2]
KEEP = '--keep' in sys.argv
OUT = os.path.join(HERE, 'colors3_family_%d.txt' % L)
if KEEP and os.path.exists(OUT): OUT = OUT[:-4] + '_N%d.txt' % N
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


t0 = time.time(); m = L ** N; k0 = L - 5
QS = [q for q in (3, 5, 7, 11, 13) if q != L][:4]; OPS = [('T', q) for q in QS] + [('T', L)]; OPSF = OPS + [('U', 2)]
say('ℓ = %d (3 colours), k₀ = %d, level 2, %d weights, modulo ℓ^%d; FLINT %s' % (L, k0, X + 1, N, flint.__version__))

# ---- 1. the class at the base weight ----
Xc, Rc, d = F2.class_on_lane2(k0, L, leg(2, L))
if d == 0: say('   the class is empty'); sys.exit()
Hb = F2.hecke_matrices2(k0, L, OPS)
ref = {q: F.restrict(Hb[('T', q)], Xc, Rc, d) for q in QS + [L]}
say('1  the class at weight %d on lane ℓ: %d form%s;  T_ℓ polynomial on lane ℓ %s' % (
    k0, d, '' if d == 1 else 's', [int(c) for c in ref[L].charpoly().coeffs()]))

# ---- 2. the family ----
rng = random.Random(L); series = []; precs = []; coefs = None


def new_part(W, Wr, rt, k, x, S):
    c2 = leg(2, L) ** (1 + x) * 2 ** ((k - 2) // 2)
    U2W = F.restrict(S[('U', 2)], W, Wr, rt)
    Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(rt)] for i in range(rt)]
    ker, Np = F.padic_kernel(Bm, rt, rt, L, N)
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
    H = F2.hecke_matrices2(k, m, OPSF)
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
        note = ' (the piece had %d forms: a form of level 1 shares the numbers on lane ℓ; kept the new part, to ℓ^%d)' % (rt, Np)
    U = F.unit_root_matrix(Ml, d, L, Np, k); cp = F.charpoly_small(U, d, L ** Np)
    series.append(cp); precs.append(Np)
    say('2  weight %d (dimension %d, ordinary part %d): the class\'s piece %d%s; U_ℓ polynomial on lane ℓ %s  (%.0f s)' % (
        k, k // 4 + 1, r, d, note, [c % L for c in cp], time.time() - t))

# ---- 3. weight −4 ----
P = min(precs); mP = L ** P; poly = []; vals = []
for j in range(d):
    v, val, _ = F.newton_at_minus_one([series[x][j] % mP for x in range(X + 1)], L, P)
    poly.append(v); vals.append(val)
poly.append(1)
say('3  the differences of each coefficient are divisible by ℓ to %s' % vals)
K = min(P, X + 1)
say('   the polynomial of u(−4) (lowest first) on lane ℓ: %s;  on lane ℓ^%d: %s' % ([c % L for c in poly], K, [c % L ** K for c in poly]))
dig = lambda v: [(v // L ** e) % L for e in range(K)]
if d == 1:
    u = (-poly[0]) % mP
    say('   u(−4), digits lowest first: %s;  −u(−4): %s' % (dig(u), dig((-u) % mP)))
else:
    say('   (several forms: the polynomial above; its factors as in partitions_flint_family.py)')
say(); say('seconds: %.0f' % (time.time() - t0))
