# ghost_series.py, 14 Nov 2026 (README "The full ghost series from the seeds, against every slope"; a relayed review:
# "a plain ghost-series computation ... compared with your 30 counted slopes").
#   1. The seeds, as serre_weights.py finds them (its code, repeated): the class forms of weights 4 … ℓ + 1 at level 6,
#      by twisted system. From them, each system's class counts at every step z: d_z (level 6, a sum over the factors of
#      Sym^(k−2)) and D_z (level 6ℓ, over Ind_B^G Sym^(k−2)), with serre_weights.py's correction for det^m′ factors.
#   2. Bergdall and Pollack's ghost series, as the review states it: g_i has a zero at step z of multiplicity
#      min(i − d_z, D_z − d_z − i) when d_z < i < D_z − d_z. At step x, v_ℓ(g_i) = Σ_z mult·(1 + v_ℓ(x − z)). The slopes of
#      the Newton polygon of (i, v_ℓ(g_i)), i = 0 … d_x, are the predicted slopes of the system's d_x forms at level 6.
#   3. Counted beside it, at level 6 directly: the class (forms6_flint, as fading_forms5.py), split by system, and the
#      slopes of T_ℓ from the characteristic polynomial on each system (exact segments only).
# DERIVED (the prediction, from the seeds alone) against COUNTED (the slopes).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/ghost_series.py ℓ X [--N=45] [--keep]
#              [--threads=t]
import os, sys, time, json
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat
import forms6_flint as F


def opt(name, default):
    return next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--' + name + '=')), default)


flint.ctx.threads = int(opt('threads', 2))
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; L, X = nums[0], nums[1]
N6 = int(opt('N', 45)); NPW = 10; KEEP = '--keep' in sys.argv
QS = [q for q in (5, 7, 11, 13) if q != L]; LM = L - 1; ZMAX = X + 6 * L
base = os.path.join(HERE, 'ghost_series_%d_x0-%d' % (L, X)); OUT, JS = base + '.txt', base + '.jsonl'
if KEEP and (os.path.exists(OUT) or os.path.exists(JS)): sys.exit('stopped before computing anything: %s exists' % OUT)
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None
jh = open(JS, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def sign(p, k): return leg(p, L) ** (1 + (k - (L - 3)) // (L - 1))      # the class's sign of U_p (fading_forms5.py)


def vl(x):
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


# ---- GL2(F_ℓ) (serre_weights.py's), with the check of Sym^n by Brauer characters up to this run's largest n ----
def jh_ind(s, t):
    s, t = s % LM, t % LM
    if s == t: return [(0, t), (LM, t)]
    r = (s - t) % LM
    return [(r, t), (LM - r, s)]


def jh_sym(n):
    out = []
    j = 0
    while True:                                          # Glover's recursion, unrolled: det^j ⊗ (the j-th layer)
        if n <= LM: out.append((n, j % LM)); return out
        out += [(r, (m + j) % LM) for r, m in jh_ind(n, 0)]
        if n < L + 1: return out
        n -= L + 1; j += 1


LL = L * L - 1
UU, VV = [a.ravel() for a in np.meshgrid(np.arange(LM), np.arange(LM), indexing='ij')]
WW = np.array([w for w in range(LL) if w % (L + 1)])
z_ = lambda e, N: np.exp(2j * np.pi * (np.asarray(e) % N) / N)
def ch_serre(r, m):
    return np.concatenate([sum(z_((m + i) * UU + (m + r - i) * VV, LM) for i in range(r + 1)),
                           sum(z_(WW * ((L + 1) * m + i + L * (r - i)), LL) for i in range(r + 1))])
def ch_sym(n):
    return np.concatenate([sum(z_(i * UU + (n - i) * VV, LM) for i in range(n + 1)),
                           sum(z_(WW * (i + L * (n - i)), LL) for i in range(n + 1))])
kX = -2 + LM * (X + 1)
bad = [n for n in range(0, kX - 1) if float(np.max(np.abs(ch_sym(n) - sum(ch_serre(*g) for g in jh_sym(n))))) > 1e-6]
say('ℓ = %d, steps x = 0 … %d (weights up to %d); the factors of Sym^n for n ≤ %d by Brauer characters: %s' % (
    L, X, kX, kX - 2, 'right' if not bad else 'WRONG at %s' % bad[:5]))
if bad: sys.exit(1)


# ---- the class at level 6, by system (fading_forms5.py's computation) ----
def split_pieces(Mq, d, Nn):
    m1 = L ** Nn; pieces = [((), F.identity(d, m1), list(range(d)), d)]
    for q in QS:
        newp = []
        for lab, Yp, Rp, rp in pieces:
            Mqr = F.restrict(Mq[q], Yp, Rp, rp)
            for f, ee in F.reduce_mat(Mqr, L).charpoly().factor()[1]:
                co = [int(c) for c in f.coeffs()]
                Bq = Mqr - F.identity(rp, m1) * ((-co[0]) % L) if len(co) == 2 else F.mat_poly(co, Mqr)
                Zm, Rz, r2 = F.kernel_part(Bq, L, Nn)
                if r2: newp.append((lab + (((-co[0]) % L) if len(co) == 2 else tuple(co),), Yp * Zm, [Rp[j] for j in Rz], r2))
        pieces = newp
    return pieces


def class_at(w, s2, s3, N, with_tl):
    """The forms of weight w on Γ0(6), new at 2 and 3 with U_2 = s2·2^((w−2)/2), U_3 = s3·3^((w−2)/2) exactly, by system:
    {label: (count, T_ℓ's characteristic polynomial mod ℓ^N1 or None)}, N1."""
    n = w + 1; m = L ** N; ops = [('T', q) for q in QS] + ([('T', L)] if with_tl else []) + [('U', 2), ('U', 3)]
    H = F.hecke_matrices(w, m, ops); c2 = s2 * 2 ** ((w - 2) // 2) % m; c3 = s3 * 3 ** ((w - 2) // 2) % m
    Y, R, r = F.kernel_part(H[('U', 2)] - F.identity(n, m) * c2, L, N)
    if r == 0: return {}, N
    Zm, Rz, r2 = F.kernel_part(F.restrict(H[('U', 3)], Y, R, r) - F.identity(r, m) * c3, L, N)
    if r2 == 0: return {}, N
    Y, R = Y * Zm, [R[j] for j in Rz]; r = r2
    U2W = F.restrict(H[('U', 2)], Y, R, r); U3W = F.restrict(H[('U', 3)], Y, R, r)
    Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)] + \
         [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(r)] for i in range(r)]
    ker, N1 = F.padic_kernel(Bm, 2 * r, r, L, N); d = len(ker)
    if d == 0: return {}, N1
    assert N1 >= 1
    m1 = L ** N1
    Xm = F.mat_mod(r, d, [ker[c][i] % m1 for i in range(r) for c in range(d)], m1)
    Rn = F.pivot_rows(nmod_mat(r, d, [ker[c][i] % L for i in range(r) for c in range(d)], L), L)
    XR = F.inverse_mod(F.mat_mod(d, d, [ker[b][Rn[a]] % m1 for a in range(d) for b in range(d)], m1), L, N1)
    def on_new(T):
        TX = F.mat_mod(r, r, [int(v) % m1 for v in T.entries()], m1) * Xm
        return XR * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m1)
    pieces = split_pieces({q: on_new(F.restrict(H[('T', q)], Y, R, r)) for q in QS}, d, N1)
    Mt = on_new(F.restrict(H[('T', L)], Y, R, r)) if with_tl else None
    out = {}
    for lab, Yp, Rp, rp in pieces:
        chi = None
        if with_tl:
            Mtp = F.restrict(Mt, Yp, Rp, rp)
            chi = F.berkowitz([[int(Mtp[i, j]) for j in range(rp)] for i in range(rp)], m1)
        out[lab] = (rp, chi)
    return out, N1


def twist(lab, m):
    out = []
    for q, v in zip(QS, lab):
        f = pow(q, m, L)
        out.append(v * f % L if isinstance(v, int) else tuple(v[i] * pow(f, len(v) - 1 - i, L) % L for i in range(len(v))))
    return tuple(out)


# ---- 1. the seeds, and the counts d_z, D_z from them ----
c2 = sign(2, L - 3) * pow(2, (L - 5) // 2, L) % L; c3 = sign(3, L - 3) * pow(3, (L - 5) // 2, L) % L
mu = {}; t0 = time.time()
for r in range(2, L, 2):
    for m in [m for m in range(LM) if (r + 2 * m - (L - 5)) % LM == 0]:
        u2 = pow(2, -m, L) * c2 % L; u3 = pow(3, -m, L) * c3 % L
        s2 = u2 * pow(2, -(r // 2), L) % L; s3 = u3 * pow(3, -(r // 2), L) % L
        ex, _ = class_at(r + 2, 1 if s2 == 1 else -1, 1 if s3 == 1 else -1, NPW, False)
        for lab, (c, _) in ex.items(): mu.setdefault(twist(lab, m), {})[(r, m)] = c
E = {lab: {mp for mp in range(LM) if all(isinstance(v, int) and v == pow(q, mp, L) * (1 + q) % L for q, v in zip(QS, lab))
           and pow(2, mp, L) == c2 and pow(3, mp, L) == c3} for lab in mu}
val = lambda lab, g: mu[lab].get(g, 0) - (1 if g[0] == 0 and g[1] in E[lab] else 0)
def counts(lab, z):
    k = -2 + LM * (z + 1)
    return (sum(val(lab, g) for g in jh_sym(k - 2)), sum(val(lab, g) for i in range(k - 1) for g in jh_ind(k - 2 - i, i)))
dD = {lab: [counts(lab, z) for z in range(ZMAX + 1)] for lab in mu}
say('the seeds (%.0f s): %s' % (time.time() - t0, '; '.join('%r: %s' % (lab, ', '.join('weight %d·ω^%d' % (r + 2, m)
    for (r, m), c in sorted(s.items()) for _ in range(c))) for lab, s in sorted(mu.items(), key=repr))))
for lab in sorted(mu, key=repr):
    say('   %r: d = %s …; D = %s …' % (lab, [dD[lab][z][0] for z in range(X + 1)], [dD[lab][z][1] for z in range(X + 1)]))


# ---- 2. the ghost slopes ----
def hull_slopes(pts):
    """Lower convex hull of the points (i, v_i), i = 0 … n: the slopes with their lengths, smallest first."""
    h = [0]
    for j in range(1, len(pts)):
        while len(h) >= 2:
            a, b = h[-2], h[-1]
            if (pts[b] - pts[a]) * (j - a) >= (pts[j] - pts[a]) * (b - a): h.pop()
            else: break
        h.append(j)
    return [(Fraction(pts[b] - pts[a], b - a), b - a) for a, b in zip(h, h[1:])]


def ghost(lab, x):
    d = dD[lab][x][0]
    v = [0] * (d + 1)
    for i in range(1, d + 1):
        for z in range(ZMAX + 1):
            dz, Dz = dD[lab][z]
            if dz < i < Dz - dz and z != x: v[i] += min(i - dz, Dz - dz - i) * (1 + vl(abs(x - z)))
    return hull_slopes(v)


def expand(sl): return sorted(s for s, n in sl for _ in range(n))


# ---- 3. counted at level 6, and the comparison ----
say(); say('at each step: per system, the ghost slopes from the seeds / the slopes of T_ℓ counted at level 6')
agree = total = 0; notexact = 0; t0 = time.time()
for x in range(0, X + 1):
    k = -2 + LM * (x + 1)
    if k < 4: say('   x = %d (weight %d): no cusp forms at level 6' % (x, k)); continue
    got, N1 = class_at(k, sign(2, k), sign(3, k), N6, True)
    labs = sorted(set(got) | {lab for lab in mu if dD[lab][x][0] > 0}, key=repr)
    parts = []; row = {'l': L, 'x': x, 'k': k, 'systems': {}}
    for lab in labs:
        cnt, chi = got.get(lab, (0, None))
        dpred = dD[lab][x][0] if lab in mu else 0
        gs = expand(ghost(lab, x)) if lab in mu and dpred > 0 else []
        cs = None
        if chi is not None:                                  # the Newton polygon of T_ℓ's characteristic polynomial
            m1 = L ** N1
            vv = [vl(c % m1) if c % m1 else N1 for c in chi]
            hh = [0]
            for j in range(1, len(vv)):
                while len(hh) >= 2:
                    a, b = hh[-2], hh[-1]
                    if (vv[b] - vv[a]) * (j - a) >= (vv[j] - vv[a]) * (b - a): hh.pop()
                    else: break
                hh.append(j)
            if all(vv[j] < N1 for j in hh):
                cs = sorted(Fraction(vv[a] - vv[b], b - a) for a, b in zip(hh, hh[1:]) for _ in range(b - a))
        if cnt == 0 and dpred == 0: continue
        if cnt != dpred: ok, note = False, '  DIFFERENT (forms predicted %d, counted %d)' % (dpred, cnt)
        elif cs is None: ok, note = None, '  (not exact at ℓ^%d)' % N1
        else: ok, note = (cs == gs), ('' if cs == gs else '  DIFFERENT')
        if ok is None: notexact += 1
        else: total += 1; agree += ok
        parts.append('%r: %s / %s%s' % (lab, [str(s) for s in gs] or '—', [str(s) for s in cs] if cs else '—', note))
        row['systems'][repr(lab)] = {'d_pred': dpred, 'd': cnt, 'ghost': [str(s) for s in gs],
                                     'counted': None if cs is None else [str(s) for s in cs]}
    say('   x = %d (weight %d): %s' % (x, k, '; '.join(parts) if parts else 'no forms'))
    if jh: jh.write(json.dumps(row) + '\n'); jh.flush(); os.fsync(jh.fileno())
say(); say('agreeing (system, step) pairs, every slope of every form: %d of %d%s  (%.0f s)' % (
    agree, total, '' if not notexact else '; %d not exact at ℓ^%d' % (notexact, N6), time.time() - t0))
