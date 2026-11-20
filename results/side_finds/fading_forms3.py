# A dig, 13 Nov 2026 (version 3 of fading_forms: several roots of slope S at once; README "What fades at 19 and 23: two
# quiet notes at once", predicted there first). As fading_forms2.py (whose header and method are kept below), with:
#   --any: at each weight, all forms new at 2 and 3 with the class's signs; the roots of slope S counted, their residues
#          u/ℓ^S on lane ℓ (the roots of the segment's polynomial), and for each simple one its Hecke numbers on lane ℓ;
#   --split (with --any): the forms with the class's signs split by their T_q numbers on lane ℓ (q = 5, 7, 11, 13 in turn;
#          rational numbers, or an irreducible factor for numbers of a bigger field), and T_ℓ's Newton polygon read in each
#          piece; a polygon's left end is marked NOT determined when its last coefficient is 0 to the working precision
#          (then the precision, not the forms, decides it: the "crowds" of fading_forms/fading_forms2 were this). A slope-S
#          segment whose two ends have valuations below the precision is exact all the same (the capped points lie left of
#          it, on steeper segments); only such segments are used. Each piece's slope-S roots are carried to weight −2 by
#          (system, residue), as below;
#   --sys=q:a,q:a/q:a,…: one or more systems ('/' between them); at each weight the piece of each system, its slope-S
#          root followed and carried to weight −2 (Lagrange in x at x = −1 through the weights known to --minp digits or
#          more, default 8, as fading_carry.py), compared with the fading notes measured from partition numbers
#          (tower_fading.txt for 13 and 17; fading_two.txt for 19 and 23: β/ℓ² = (3|ℓ)·u(−2)/ℓ², by first digit).
# The Hecke matrices are made once per weight for all systems.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/fading_forms3.py ℓ X N [--extra=E] [--slope=S]
#          [--any | --sys=…] [--skip] [--minp=P] [--keep] [--threads=t]
# ---- fading_forms2.py's method (abridged from its header): the system on lane ℓ is the weight-(ℓ − 3) class form; at each
# weight k = −2 + (ℓ − 1)·ℓ·(x + 1), modulo ℓ^(N + E): the piece with the system's T_q numbers (q = 5, 7, 11, 13, not ℓ),
# inside it the forms new at 2 and 3 with the class's U₂, U₃ values, T_ℓ on them, its characteristic polynomial
# (Berkowitz), its Newton polygon, the root of slope S (Hensel), u from it (u = ℓ^S·w, w² − (λ/ℓ^S)·w + ℓ^(k−1−2S) = 0).
# No partition number is used. COUNTED.
import os, sys, time
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat, nmod_poly
import forms6_flint as F
flint.ctx.threads = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--threads=')), 8)
args = [int(a) for a in sys.argv[1:] if a.isdigit()]; L, X, N = args[0], args[1], args[2]
EXTRA = next((int(a[8:]) for a in sys.argv[1:] if a.startswith('--extra=')), 30); NP = N + EXTRA
KEEP = '--keep' in sys.argv
S = next((int(a[8:]) for a in sys.argv[1:] if a.startswith('--slope=')), 2)
ANY = '--any' in sys.argv
SPLIT = ANY and '--split' in sys.argv
SKIP = '--skip' in sys.argv
MINP = next((int(a[7:]) for a in sys.argv[1:] if a.startswith('--minp=')), 8)
X0 = next((int(a[7:]) for a in sys.argv[1:] if a.startswith('--from=')), 0)   # --from=x: start at weight x (resume a stopped run)
SMALL = '--small' in sys.argv          # --small: weights −2 + (ℓ − 1)·(x + 1), only ℓ − 1 apart (outside the family's guaranteed disc)
WSTEP = (lambda L_: L_ - 1) if SMALL else (lambda L_: (L_ - 1) * L_)
SYSS = next(([{int(kv.split(':')[0]): int(kv.split(':')[1]) for kv in part.split(',')} for part in a[6:].split('/')]
             for a in sys.argv[1:] if a.startswith('--sys=')), [])
OUT = os.path.join(HERE, 'fading_forms3_%d_s%d%s%s.txt' % (L, S, '_split' if SPLIT else ('_any' if ANY else ('_sys' if SYSS else '')), '_small' if SMALL else ''))
if KEEP and os.path.exists(OUT): OUT = OUT[:-4] + '_N%d_X%d%s.txt' % (N, X, '_from%d' % X0 if X0 else '')
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None
SETTLED = {13: {10: [10, 11, 6, 2, 10, 10]}, 17: {3: [3, 9, 1, 13, 4, 3]},                  # tower_fading.txt
           19: {2: [2, 1, 0, 17], 13: [13, 14, 2, 15]}, 23: {8: [8, 20, 21, 8], 17: [17, 20, 19, 21]},   # fading_two.txt
           29: {23: [23, 3, 10]}, 31: {2: [2, 28, 13]}}                                       # tower_fading.txt (slope 1)


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def segment(chi, ell, s, Nn):
    """The Newton polygon of chi (lowest first, modulo ℓ^Nn): its slopes, and the segment of slope s as the residues on
    lane ℓ of root/ℓ^s for its roots (the roots of h(y) = chi(ℓ^s·y)/ℓ^c between the segment's ends), with multiplicity."""
    d = len(chi) - 1; m = ell ** Nn
    def val(x):
        x %= m
        if x == 0: return Nn
        v = 0
        while x % ell == 0: x //= ell; v += 1
        return v
    vs = [val(c) for c in chi]; hull = [0]
    for j in range(1, d + 1):
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            if (vs[b] - vs[a]) * (j - a) >= (vs[j] - vs[a]) * (b - a): hull.pop()
            else: break
        hull.append(j)
    slopes = [(Fraction(vs[a] - vs[b], b - a), b - a) for a, b in zip(hull, hull[1:])]
    seg = [(a, b) for a, b in zip(hull, hull[1:]) if Fraction(vs[a] - vs[b], b - a) == s]
    if not seg: return (slopes, 0, []) if not EXACT else (slopes, 0, [], True)
    j1, j2 = seg[0]; c = vs[j1] + s * j1
    res = [((chi[j] % m) // ell ** (c - s * j)) % ell if c - s * j >= 0 else 0 for j in range(j1, j2 + 1)]
    rts = sorted((int(r), e) for r, e in nmod_poly(res, ell).roots())
    # the segment is exact when its two ends have valuations below the precision: points capped at the precision lie
    # left of it, on steeper segments, and their true (higher) valuations cannot change it
    return (slopes, j2 - j1, rts) if not EXACT else (slopes, j2 - j1, rts, vs[j1] < Nn and vs[j2] < Nn)


EXACT = False                                                          # segment() also says whether its slope-S part is exact


def sign(p, k): return leg(p, L) ** (1 + (k - (L - 3)) // (L - 1))      # the class's sign of U_p at weight k


def compose(Y, R, Z, Rz): return Y * Z, [R[j] for j in Rz]


QS = [q for q in (5, 7, 11, 13) if q != L]; t0 = time.time()
say('ℓ = %d, slope %d, weights −2 + %d·(x + 1) for x = 0..%d, modulo ℓ^%d (N = %d and %d spare); FLINT %s' % (
    L, S, WSTEP(L), X, NP, N, EXTRA, flint.__version__))
kr = L - 3
if SPLIT and L > 100:                  # split mode never uses the class's own system; a class of many forms (607) has none
    class_ref = None
    say('the class at weight %d is not read here (split mode; the class has many forms)' % kr)
else:
    Kz, dr = F.class_basis_exact(kr, sign(2, kr), sign(3, kr))
    Xr, Rr = F.normalized_basis(nmod_mat(kr + 1, dr, [int(Kz[i, j]) % L for i in range(kr + 1) for j in range(dr)], L), dr, L)
    Hr = F.hecke_matrices(kr, L, [('T', q) for q in QS])
    if dr != 1: say('   the class at weight %d has %d forms; stopped' % (kr, dr)); sys.exit()
    class_ref = {q: int(F.restrict(Hr[('T', q)], Xr, Rr, 1)[0, 0]) for q in QS}
    say('the system on lane ℓ: the weight-%d form, the only one of its class; T_q numbers %s' % (kr, class_ref))
if SYSS: say('   --sys: the pieces with the Hecke numbers %s' % ' and '.join(str(sy) for sy in SYSS))
carry = {i: [] for i in range(len(SYSS))}

for x in range(X0, X + 1):
    k = -2 + WSTEP(L) * (x + 1); t = time.time(); m = L ** NP; n = k + 1
    H = F.hecke_matrices(k, m, [('T', q) for q in QS] + [('T', L), ('U', 2), ('U', 3)])
    c2 = sign(2, k) * 2 ** ((k - 2) // 2) % m; c3 = sign(3, k) * 3 ** ((k - 2) // 2) % m
    for si, ref in ([(None, None)] if ANY else list(enumerate(SYSS))):
        Y, R, r, dims = None, None, n, []
        for q in ([] if ref is None else QS):
            Tq = H[('T', q)] if Y is None else F.restrict(H[('T', q)], Y, R, r)
            Z, Rz, r2 = F.kernel_part(Tq - F.identity(r, m) * ref[q], L, NP)
            if r2 == 0: r = 0; break
            if Y is None: Y, R = Z, Rz
            else: Y, R = compose(Y, R, Z, Rz)
            r = r2; dims.append(r)
        tag = 'every system' if ref is None else 'the system %s' % ref
        if r == 0: say('   weight %d, %s: no form with these numbers' % (k, tag)); continue
        if Y is None: Z, Rz, r2 = F.kernel_part(H[('U', 2)] - F.identity(r, m) * c2, L, NP); Y, R = Z, Rz
        else: Z, Rz, r2 = F.kernel_part(F.restrict(H[('U', 2)], Y, R, r) - F.identity(r, m) * c2, L, NP); Y, R = (compose(Y, R, Z, Rz) if r2 else (None, None))
        if r2 == 0: say('   weight %d, %s: no form with the class\'s U₂ value' % (k, tag)); continue
        r = r2
        Z, Rz, r2 = F.kernel_part(F.restrict(H[('U', 3)], Y, R, r) - F.identity(r, m) * c3, L, NP)
        if r2 == 0: say('   weight %d, %s: no form with the class\'s U₃ value' % (k, tag)); continue
        Y, R = compose(Y, R, Z, Rz); r = r2
        U2W = F.restrict(H[('U', 2)], Y, R, r); U3W = F.restrict(H[('U', 3)], Y, R, r); TLW = F.restrict(H[('T', L)], Y, R, r)
        Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)] + \
             [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(r)] for i in range(r)]
        ker, N1 = F.padic_kernel(Bm, 2 * r, r, L, NP); d = len(ker); m1 = L ** N1
        if d == 0: say('   weight %d, %s: no form new at 2 and 3' % (k, tag)); continue
        Xn = [[ker[c][i] % m1 for c in range(d)] for i in range(r)]
        Xm = F.mat_mod(r, d, [Xn[i][c] for i in range(r) for c in range(d)], m1)
        Rn = F.pivot_rows(nmod_mat(r, d, [Xn[i][c] % L for i in range(r) for c in range(d)], L), L)
        XR = F.inverse_mod(F.mat_mod(d, d, [Xn[Rn[a]][b] for a in range(d) for b in range(d)], m1), L, N1)
        def on_new(T):
            TX = F.mat_mod(r, r, [int(v) % m1 for v in T.entries()], m1) * Xm
            return XR * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m1)
        Mt = on_new(TLW)
        chi = F.berkowitz([[int(Mt[i, j]) for j in range(d)] for i in range(d)], m1)
        slopes, nS, rts = segment(chi, L, S, N1)
        poly = ', '.join('%s×%d' % (sl, ln) for sl, ln in slopes)
        say('   weight %d (dimension %d), %s: the piece %s, new at 2 and 3 with the class\'s signs: %d (to ℓ^%d); slopes of T_ℓ: %s%s; '
            '%d root%s of slope %d, u/ℓ^%d on lane ℓ: %s  (%.0f s)' % (k, n, tag, dims + [r], d, N1, poly,
                                                                 '' if (not SPLIT or chi[0] % m1) else ' (NOT determined: the last coefficient is 0 to ℓ^%d)' % N1,
                                                                 nS, '' if nS == 1 else 's', S, S,
                                                                 ', '.join('%d%s' % (y, '' if e == 1 else ' (×%d)' % e) for y, e in rts) or '—', time.time() - t))
        if SPLIT:
            Mq = {q: on_new(F.restrict(H[('T', q)], Y, R, r)) for q in QS}
            pieces = [((), F.identity(d, m1), list(range(d)), d)]
            for q in QS:
                newp = []
                for lab, Yp, Rp, rp in pieces:
                    Mqr = F.restrict(Mq[q], Yp, Rp, rp)
                    for f, ee in F.reduce_mat(Mqr, L).charpoly().factor()[1]:
                        co = [int(c) for c in f.coeffs()]
                        B = Mqr - F.identity(rp, m1) * ((-co[0]) % L) if len(co) == 2 else F.mat_poly(co, Mqr)
                        Z, Rz, r2 = F.kernel_part(B, L, N1)
                        if r2 == 0: continue
                        Y2, R2 = compose(Yp, Rp, Z, Rz)
                        newp.append((lab + (((-co[0]) % L) if len(co) == 2 else tuple(co),), Y2, R2, r2))
                pieces = newp
            for lab, Yp, Rp, rp in pieces:
                Mtp = F.restrict(Mt, Yp, Rp, rp)
                chip = F.berkowitz([[int(Mtp[i, j]) for j in range(rp)] for i in range(rp)], m1)
                EXACT = True; sl_p, nSp, rtsp, exact = segment(chip, L, S, N1); EXACT = False
                dtm = chip[0] % m1 != 0
                say('      the system %s: %d form%s; slopes of T_ℓ: %s%s; %d root%s of slope %d%s%s' % (
                    dict(zip(QS, lab)), rp, '' if rp == 1 else 's', ', '.join('%s×%d' % (a_, b_) for a_, b_ in sl_p),
                    '' if dtm else ' (its left end NOT determined: the last coefficient is 0 to ℓ^%d)' % N1, nSp, '' if nSp == 1 else 's', S,
                    (', u/ℓ^%d on lane ℓ: %s' % (S, ', '.join('%d%s' % (y, '' if e == 1 else ' (×%d)' % e) for y, e in rtsp))) if nSp else '',
                    ('' if exact else ' — NOT exact, skipped') if nSp else ''))
                if nSp and not exact: continue
                for y0, e in rtsp:
                    if e != 1: continue
                    lam, pr, _ = F.slope_root(chip, L, S, N1, pick=y0 if nSp > 1 else None)
                    if lam is None: continue
                    mp = L ** (pr - S); w = (lam // L ** S) % mp; ee = k - 1 - 2 * S
                    for _ in range(pr + 2): w = ((lam // L ** S) - (L ** ee % mp) * pow(w, -1, mp)) % mp if ee < pr - S else (lam // L ** S) % mp
                    dg = [((leg(3, L) * w % L ** (pr - S)) // L ** i) % L for i in range(min(8, pr - S))]
                    say('         the root with u/ℓ^%d ≡ %d: known to ℓ^%d; (3|ℓ)·u/ℓ^%d, digits lowest first: %s' % (S, y0, pr - S, S, dg))
                    carry.setdefault((lab, 0 if SMALL else y0), []).append((x, w, pr - S))   # --small: the residue moves with the weight
            continue
        if nS == 0: continue
        Mq = {q: on_new(F.restrict(H[('T', q)], Y, R, r)) for q in QS}
        for y0, e in rts:
            if e != 1: continue
            lam, pr, _ = F.slope_root(chi, L, S, N1, pick=y0 if nS > 1 else None)
            if lam is None: say('      the root ≡ %d: not found alone' % y0); continue
            mp = L ** (pr - S); w = (lam // L ** S) % mp; ee = k - 1 - 2 * S
            for _ in range(pr + 2): w = ((lam // L ** S) - (L ** ee % mp) * pow(w, -1, mp)) % mp if ee < pr - S else (lam // L ** S) % mp
            mpr = L ** min(pr, N1); Bk = [[(int(Mt[i, j]) - (lam if i == j else 0)) % mpr for j in range(d)] for i in range(d)]
            kv, Nk = F.padic_kernel(Bk, d, d, L, min(pr, N1))
            sysq = {}
            if len(kv) == 1:
                wv = kv[0]; i0 = next(i for i in range(d) if wv[i] % L)
                for q in QS: sysq[q] = sum(int(Mq[q][i0, j]) * wv[j] for j in range(d)) % L * pow(wv[i0] % L, -1, L) % L
            twists = [j for j in range(L - 1) if sysq and all(sysq[q] == pow(q, j, L) * class_ref[q] % L for q in QS)]
            dg = [((leg(3, L) * w % L ** (pr - S)) // L ** i) % L for i in range(min(8, pr - S))]
            say('      the root with u/ℓ^%d ≡ %d: known to ℓ^%d; (3|ℓ)·u/ℓ^%d, digits lowest first: %s; its Hecke numbers on lane ℓ: %s; '
                'θ-twists of the class: %s' % (S, y0, pr - S, S, dg, sysq or '(eigenvector not alone)', twists or 'none'))
            if ref is not None and nS == 1: carry[si].append((x, w, pr - S))

# ---- weight −2 ----
for si, rows in carry.items():
    name = SYSS[si] if not isinstance(si, tuple) else ('%s, its slope-%d root' % (dict(zip(QS, si[0])), S) if SMALL else '%s, root ≡ %d' % (dict(zip(QS, si[0])), si[1]))
    say(); say('the system %s: weights carried (x, digits known): %s' % (name, [(xx, pp) for xx, _, pp in rows]))
    good = []
    for row in rows:
        if row[2] >= MINP: good.append(row)
        else: break
    for nn in range(2, len(good) + 1):
        P = min(pp for _, _, pp in good[:nn]); mP = L ** P; v = 0; xs = [xx for xx, _, _ in good[:nn]]
        for i, xi in enumerate(xs):
            num = den = 1
            for j, xj in enumerate(xs):
                if j != i: num *= (-1 - xj); den *= (xi - xj)
            v = (v + good[i][1] % mP * num * pow(den % mP, -1, mP)) % mP
        dg = [((leg(3, L) * v % mP) // L ** j) % L for j in range(P)]
        target = SETTLED.get(L, {}).get(dg[0], [])
        same = next((i for i in range(min(len(dg), len(target))) if dg[i] != target[i]), min(len(dg), len(target)))
        say('   through %d weights (each known to ℓ^%d or more): (3|ℓ)·u(−2)/ℓ^%d = %s …;  measured on partition numbers: %s  →  %d of %d digits agree' % (
            nn, MINP, S, dg[:6], target or '(no note with this first digit)', same, len(target)))
say(); say('seconds: %.0f' % (time.time() - t0))
