# A dig, 12 Nov 2026: THE TOWER'S CONSTANTS FROM MODULAR FORMS, ON FLINT ALONE (no PARI; engine: forms6_flint.py).
# 1  At the base weight k₀ = ℓ − 3: the tower's class in M_k₀(Γ0(6)), found over the field with ℓ elements as the joint
#    generalized eigenspace of U₂ and U₃ for the eigenvalues (2|ℓ)·2^((k₀−2)/2) and (3|ℓ)·3^((k₀−2)/2) (the class's
#    Atkin–Lehner signs, read off 13 to 31 earlier). PREDICTED: its dimension is ⌊(ℓ + 12)/24⌋ and the polynomials of
#    (3|Q)·a_Q on lane ℓ are those PARI gave (partitions_split.txt, partitions_three_forms.txt).
# 2  At the weights k = k₀ + (ℓ − 1)·x, x = 0..X: the piece of M_k(Γ0(6)) with the class's Hecke numbers on lane ℓ (for a
#    combination A of T₅, T₇, T₁₁, T_ℓ), lifted to modulo ℓ^N; for ℓ ≥ 7 no lower-level or Eisenstein form shares them (a
#    level-lowering congruence would need a_p ≡ ±(p + 1)p^((k−2)/2), impossible at p = 2 and, for ℓ ≠ 5, at p = 3).
#    On that piece: the U_ℓ-action of the ordinary stabilizations and its characteristic polynomial.
# 3  Newton's series of each coefficient at x = −1 (weight −2), the twist by (3|ℓ): the polynomial of the constants, its
#    factors over the ℓ-adic integers (FLINT mod ℓ, Hensel), and the constants themselves (one on lane ℓ, or a pair
#    A ± B·√r with r a non-square mod ℓ). PREDICTED: they agree with every earlier result (13: partitions_hida_digits_x12;
#    37, 41, 43: partitions_twins; 47, 61: on lane ℓ, partitions_split and partitions_three_forms).
# 4  (added 13 Nov 2026 in this copy) The same for the step numbers at Q = 5: the polynomial of T₅ on the family's piece at
#    each weight, Newton's series at x = −1, the twist by (3|5): λ(5) = (3|5)·a₅(F at weight −2) to N digits, split over the
#    ℓ-adic integers as the constants are (README "97's tower is two towers", the second digit). PREDICTED there first.
# COUNTED.  usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_family_steps.py ℓ X N [--keep]
import os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly, nmod_mat
import forms6_flint as F
flint.ctx.threads = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--threads=')), 8)
args = [int(a) for a in sys.argv[1:] if a.isdigit()]; L, X, N = args[0], args[1], args[2]
KEEP = '--keep' in sys.argv
OUT = os.path.join(HERE, 'partitions_flint_family_steps_%d.txt' % L)
if KEEP and os.path.exists(OUT): OUT = OUT[:-4] + '_N%d.txt' % N      # a later, longer run keeps its own file
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


t0 = time.time(); m = L ** N; k0 = L - 3; dexp = (L + 12) // 24
QS = [q for q in (5, 7, 11, 13) if q != L]; OPS = [('T', q) for q in QS] + [('T', L)]
say('ℓ = %d, k₀ = %d, %d weights, modulo ℓ^%d; FLINT %s' % (L, k0, X + 1, N, flint.__version__))

# ---- 1. the class at the base weight ----
t = time.time(); H = F.hecke_matrices(k0, m, OPS + [('U', 2), ('U', 3)]); n0 = k0 + 1
c2 = leg(2, L) * pow(2, (k0 - 2) // 2, L) % L; c3 = leg(3, L) * pow(3, (k0 - 2) // 2, L) % L
if k0 <= 150:                                              # exactly, over the integers (a form of lower level can share
    Kz, r0 = F.class_basis_exact(k0, leg(2, L), leg(3, L))  # the class's U₂, U₃ values on lane ℓ: it did for 17 and 29)
    Xc0 = nmod_mat(n0, r0, [int(Kz[i, j]) % L for i in range(n0) for j in range(r0)], L)
    if Xc0.rank() != r0: say('   the exact class does not stay of full rank on lane ℓ; stopped'); sys.exit()
    how = 'exactly, over the integers'
else:
    U2 = F.reduce_mat(H[('U', 2)], L); U3 = F.reduce_mat(H[('U', 3)], L); I0 = F.identity(n0, L)
    K1, r1 = F.generalized_kernel(U2 - I0 * c2); X1, R1 = F.normalized_basis(K1, r1, L)
    M3 = F.restrict(U3, X1, R1, r1); K2, r0 = F.generalized_kernel(M3 - F.identity(r1, L) * c3)
    Xc0 = X1 * nmod_mat(r1, r0, [int(K2[i, j]) for i in range(r1) for j in range(r0)], L)
    how = 'on lane ℓ'
Xc, Rc = F.normalized_basis(Xc0, r0, L)
say('1  the class at weight %d (found %s): dimension %d (predicted ⌊(ℓ + 12)/24⌋ = %d)  →  %s  (%.0f s)' % (
    k0, how, r0, dexp, 'as predicted' if r0 == dexp else 'NOT as predicted', time.time() - t))
if r0 != dexp: sys.exit()
d = r0; ref = {q: F.restrict(F.reduce_mat(H[('T', q)], L), Xc, Rc, d) for q in QS + [L]}
for q in QS + [L]:
    cp = F.charpoly_small(ref[q], d, L); s = leg(3, q) if q != 3 else 1
    say('   (3|Q)·a_Q on lane ℓ, Q = %d: polynomial %s' % (q, [cp[j] * s ** (d - j) % L for j in range(d + 1)]))

# ---- 2. the family ----
rng = random.Random(L); series = []; precs = []; series5 = []
coefs = None; OPSF = OPS + [('U', 2), ('U', 3)]; s2b, s3b = leg(2, L), leg(3, L)


def new_part(W, Wr, rt, k, x):
    """Inside the lifted piece (rank rt > d): the forms new at 2 and 3, whose U₂, U₃ act EXACTLY as the class's values
    s·2^((k−2)/2), s′·3^((k−2)/2) (signs carried along the family: s₂(k) = s₂·(2|ℓ)^x); an old form that shares the
    numbers on lane ℓ (a level-raising congruence) differs from them deeper. Returns {op: d×d matrix mod ℓ^N′}, N′."""
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
    out = {}
    for op in (('T', QS[0]), ('T', L)):
        TW = F.restrict(S[op], W, Wr, rt); TWp = F.mat_mod(rt, rt, [int(v) % mp_ for v in TW.entries()], mp_)
        TX = TWp * Xn
        out[op] = XnRinv * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], mp_)
    return out, Np


for x in range(X + 1):
    k = k0 + (L - 1) * x; t = time.time()
    H = F.hecke_matrices(k, m, OPSF) if x > 0 else H
    Y, Ry, r = F.ordinary_part(H[('T', L)], L, N)           # the ordinary part: small, about constant along the family
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
    if rt == d:
        Ml = F.restrict(S[('T', L)], Xk, Rk, d); M5 = F.restrict(S[('T', QS[0])], Xk, Rk, d)
    else:
        nw, Np = new_part(Xk, Rk, rt, k, x)
        if nw is None: say('   weight %d: the piece has dimension %d and its new part %d, not %d; stopped' % (k, rt, Np, d)); sys.exit()
        Ml, M5 = nw[('T', L)], nw[('T', QS[0])]
        note = ' (the piece had dimension %d: an old form shares the numbers on lane ℓ; kept the new part, to ℓ^%d)' % (rt, Np)
    ok5 = F.charpoly_small(F.reduce_mat(M5, L), d, L) == F.charpoly_small(ref[QS[0]], d, L)
    U = F.unit_root_matrix(Ml, d, L, Np, k); cp = F.charpoly_small(U, d, L ** Np)
    series.append(cp); precs.append(Np); series5.append(F.charpoly_small(M5, d, L ** Np))
    say('2  weight %d (dimension %d, ordinary part %d): the class\'s piece has dimension %d (T_%d on lane ℓ as at the base: %s)%s; U_ℓ polynomial on lane ℓ: %s  (%.0f s)' % (
        k, k + 1, r, d, QS[0], 'yes' if ok5 else 'NO', note, [c % L for c in cp], time.time() - t))

# ---- 3. weight −2 ----
s3 = leg(3, L); poly = []; vals = []; N = min(precs); m = L ** N
for j in range(d):
    v, val, _ = F.newton_at_minus_one([series[x][j] % m for x in range(X + 1)], L, N)
    poly.append(v * s3 ** (d - j) % m); vals.append(val)
poly.append(1)
say('3  the differences of each coefficient are divisible by ℓ to %s' % vals)
K = min(N, X + 1); MK = L ** K
say('   the constants\' polynomial (lowest first) on lane ℓ: %s;  on lane ℓ^%d: %s' % ([c % L for c in poly], K, [c % MK for c in poly]))
fac = nmod_poly([c % L for c in poly], L).factor()[1]
rest = [c % m for c in poly]; pieces = []
for fp, e in fac:
    gl = [int(c) for c in fp.coeffs()]
    if e != 1: say('   a repeated factor mod ℓ; not split further'); break
    if len(fac) == 1: pieces.append(rest); break
    hbar = nmod_poly([c % L for c in rest], L) // fp
    G, Hh = F.hensel_factor(rest, gl, [int(c) for c in hbar.coeffs()], L, N)
    pieces.append(G); rest = Hh
    if len(pieces) == len(fac) - 1: pieces.append(rest); break
dig = lambda v: [(v // L ** e) % L for e in range(K)]
for pc in pieces:
    deg = len(pc) - 1
    if deg == 1: say('   a constant on lane ℓ, digits lowest first: %s' % dig((-pc[0]) % m))
    elif deg == 2:
        c0, c1 = pc[0] % m, pc[1] % m; r = next(r for r in range(2, L) if leg(r, L) == -1)
        disc = (c1 * c1 - 4 * c0) % m
        if leg(disc, L) == 1:
            sq = F.padic_sqrt(disc, L, N); inv2 = pow(2, -1, m)
            say('   two constants on lane ℓ: %s and %s' % (dig((-c1 + sq) * inv2 % m), dig((-c1 - sq) * inv2 % m)))
        else:
            w = F.padic_sqrt(disc * pow(r, -1, m) % m, L, N); inv2 = pow(2, -1, m)
            say('   two constants A ± B·√%d (no square root of %d on lane ℓ): A %s, B %s' % (r, r, dig((-c1) * inv2 % m), dig(w * inv2 % m)))
    else: say('   a factor of degree %d: %s' % (deg, [c % MK for c in pc]))
# ---- 4. the step numbers at Q = 5, at weight −2 ----
s5 = leg(3, QS[0]); poly5 = []; vals5 = []
for j in range(d):
    v, val, _ = F.newton_at_minus_one([series5[x][j] % m for x in range(X + 1)], L, N)
    poly5.append(v * s5 ** (d - j) % m); vals5.append(val)
poly5.append(1)
say('4  the step numbers at Q = %d, λ = (3|Q)·a_Q(F at weight −2): the differences of each coefficient are divisible by ℓ to %s' % (QS[0], vals5))
say('   their polynomial (lowest first) on lane ℓ: %s;  on lane ℓ^%d: %s' % ([c % L for c in poly5], K, [c % MK for c in poly5]))
for fp, e in nmod_poly([c % L for c in poly5], L).factor()[1]:
    if fp.degree() == 1 and e == 1:
        hb = nmod_poly([c % L for c in poly5], L) // fp
        Gp, Hp = F.hensel_factor([c % m for c in poly5], [int(c) for c in fp.coeffs()], [int(c) for c in hb.coeffs()], L, N)
        say('   a root on lane ℓ, digits lowest first: %s;  the rest (lowest first) on lane ℓ^%d: %s' % (dig((-Gp[0]) % m), K, [c % MK for c in Hp]))
say(); say('seconds: %.0f' % (time.time() - t0))
