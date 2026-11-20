# A dig, 13 Nov 2026: FOLLOWING A WHOLE GROUP OF NOTES (README "607's quiet braids"; predicted there first). As
# fading_forms3.py --any --split (its method kept), but for every piece of forms with one system on lane ℓ the whole set of
# its roots of slope S is followed: the factor G(y) of T_ℓ's characteristic polynomial whose roots are root/ℓ^S for the
# roots of valuation exactly S (Hensel's lemma on the Newton segment; exact when the segment's two ends are below the
# precision). Its coefficients are carried from weight to weight by the piece's system and taken to weight −2 by Lagrange
# at x = −1. Printed for weight −2: G mod ℓ^P, its factors on lane ℓ, and the polynomial of the constants (3|ℓ)·u/ℓ^S.
# Weights −2 + (ℓ − 1)·(x + 1) (--small, the default here) or −2 + ℓ(ℓ − 1)·(x + 1) (--disc). No partition number is
# used. COUNTED (the forms), DERIVED (the value at −2).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/fading_forms4.py ℓ X N [--slope=S]
#          [--from=x] [--extra=E] [--minp=P] [--disc] [--keep] [--threads=t]
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
EXTRA = next((int(a[8:]) for a in sys.argv[1:] if a.startswith('--extra=')), 20); NP = N + EXTRA
S = next((int(a[8:]) for a in sys.argv[1:] if a.startswith('--slope=')), 1)
X0 = next((int(a[7:]) for a in sys.argv[1:] if a.startswith('--from=')), 1)
MINP = next((int(a[7:]) for a in sys.argv[1:] if a.startswith('--minp=')), 4)
DISC = '--disc' in sys.argv; KEEP = '--keep' in sys.argv
WSTEP = (L - 1) * L if DISC else L - 1
OUT = os.path.join(HERE, 'fading_forms4_%d_s%d%s.txt' % (L, S, '_disc' if DISC else ''))
if KEEP and os.path.exists(OUT): OUT = OUT[:-4] + '_N%d_X%d.txt' % (N, X)
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def sign(p, k): return leg(p, L) ** (1 + (k - (L - 3)) // (L - 1))      # the class's sign of U_p at weight k


def compose(Y, R, Z, Rz): return Y * Z, [R[j] for j in Rz]


def segment_factor(chi, s, Nn):
    """The Newton polygon of chi (lowest first, mod ℓ^Nn); for its segment of slope s: (G monic lowest first, whose roots
    are root/ℓ^s for the roots of valuation s, mod ℓ^prec; prec; exact?; the polygon's slopes)."""
    d = len(chi) - 1; m = L ** Nn
    def val(x):
        x %= m
        if x == 0: return Nn
        v = 0
        while x % L == 0: x //= L; v += 1
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
    if not seg: return None, 0, False, slopes
    j1, j2 = seg[0]; exact = vs[j1] < Nn and vs[j2] < Nn; c = vs[j1] + s * j1; prec = Nn - c
    if not exact or prec <= 0: return None, 0, exact, slopes
    mp = L ** prec
    h = [((chi[j] % m) * L ** (s * j - c)) % mp if s * j - c >= 0 else ((chi[j] % m) // L ** (c - s * j)) % mp for j in range(d + 1)]
    lc = h[j2] % L; inv = pow(lc, -1, L)
    g = [h[j] % L * inv % L for j in range(j1, j2 + 1)]
    G, _ = F.hensel_factor(h, g, [0] * j1 + [lc], L, prec)
    return [x % mp for x in G], prec, True, slopes


QS = [q for q in (5, 7, 11, 13) if q != L]; t0 = time.time(); carry = {}
say('ℓ = %d, slope %d, weights −2 + %d·(x + 1) for x = %d..%d, modulo ℓ^%d; FLINT %s' % (L, S, WSTEP, X0, X, NP, flint.__version__))
for x in range(X0, X + 1):
    k = -2 + WSTEP * (x + 1); t = time.time(); m = L ** NP; n = k + 1
    H = F.hecke_matrices(k, m, [('T', q) for q in QS] + [('T', L), ('U', 2), ('U', 3)])
    c2 = sign(2, k) * 2 ** ((k - 2) // 2) % m; c3 = sign(3, k) * 3 ** ((k - 2) // 2) % m
    Y, R, r = F.kernel_part(H[('U', 2)] - F.identity(n, m) * c2, L, NP)
    Z, Rz, r2 = F.kernel_part(F.restrict(H[('U', 3)], Y, R, r) - F.identity(r, m) * c3, L, NP); Y, R = compose(Y, R, Z, Rz); r = r2
    U2W = F.restrict(H[('U', 2)], Y, R, r); U3W = F.restrict(H[('U', 3)], Y, R, r)
    Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)] + \
         [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(r)] for i in range(r)]
    ker, N1 = F.padic_kernel(Bm, 2 * r, r, L, NP); d = len(ker); m1 = L ** N1
    Xn = [[ker[c][i] % m1 for c in range(d)] for i in range(r)]
    Xm = F.mat_mod(r, d, [Xn[i][c] for i in range(r) for c in range(d)], m1)
    Rn = F.pivot_rows(nmod_mat(r, d, [Xn[i][c] % L for i in range(r) for c in range(d)], L), L)
    XR = F.inverse_mod(F.mat_mod(d, d, [Xn[Rn[a]][b] for a in range(d) for b in range(d)], m1), L, N1)
    def on_new(T):
        TX = F.mat_mod(r, r, [int(v) % m1 for v in T.entries()], m1) * Xm
        return XR * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m1)
    Mt = on_new(F.restrict(H[('T', L)], Y, R, r)); Mq = {q: on_new(F.restrict(H[('T', q)], Y, R, r)) for q in QS}
    pieces = [((), F.identity(d, m1), list(range(d)), d)]
    for q in QS:
        newp = []
        for lab, Yp, Rp, rp in pieces:
            Mqr = F.restrict(Mq[q], Yp, Rp, rp)
            for f, ee in F.reduce_mat(Mqr, L).charpoly().factor()[1]:
                co = [int(c) for c in f.coeffs()]
                B = Mqr - F.identity(rp, m1) * ((-co[0]) % L) if len(co) == 2 else F.mat_poly(co, Mqr)
                Z, Rz, r2 = F.kernel_part(B, L, N1)
                if r2: Y2, R2 = compose(Yp, Rp, Z, Rz); newp.append((lab + (((-co[0]) % L) if len(co) == 2 else tuple(co),), Y2, R2, r2))
        pieces = newp
    say('   weight %d (dimension %d): %d forms new at 2 and 3 with the class\'s signs (to ℓ^%d), in %d systems  (%.0f s)' % (k, n, d, N1, len(pieces), time.time() - t))
    for lab, Yp, Rp, rp in pieces:
        Mtp = F.restrict(Mt, Yp, Rp, rp)
        chip = F.berkowitz([[int(Mtp[i, j]) for j in range(rp)] for i in range(rp)], m1)
        G, prec, exact, sl = segment_factor(chip, S, N1)
        degs = [1 if not isinstance(v, tuple) else len(v) - 1 for v in lab]
        line = '      a system of degree %s on lane ℓ: %d form%s; slopes %s' % (degs[0], rp, '' if rp == 1 else 's', ', '.join('%s×%d' % (a_, b_) for a_, b_ in sl))
        if G is not None:
            Gb = nmod_poly([c % L for c in G], L); fdeg = sorted(ff.degree() for ff, e in Gb.factor()[1] for _ in range(e))
            line += ';  its %d root%s of slope %d: G known to ℓ^%d, factors on lane ℓ of degrees %s' % (len(G) - 1, '' if len(G) == 2 else 's', S, prec, fdeg)
            carry.setdefault(lab, []).append((x, G, prec))
        elif exact is False and any(a_ == S for a_, b_ in sl): line += ';  its slope-%d segment NOT exact at this precision' % S
        say(line)
say()
s3 = leg(3, L)
for lab, rows in carry.items():
    good = [rw for rw in rows if rw[2] >= MINP]
    if len(good) < 2 or len({len(rw[1]) for rw in good}) != 1: continue
    deg = len(good[0][1]) - 1; P = min(rw[2] for rw in good); mP = L ** P; xs = [rw[0] for rw in good]
    Gm2 = []
    for i in range(deg + 1):
        v = 0
        for a, xa in enumerate(xs):
            num = den = 1
            for b, xb in enumerate(xs):
                if b != a: num *= (-1 - xb); den *= (xa - xb)
            v = (v + good[a][1][i] * num * pow(den % mP, -1, mP)) % mP
        Gm2.append(v)
    Gs = [Gm2[i] * pow(s3, deg - i, mP) % mP for i in range(deg + 1)]      # the polynomial of (3|ℓ)·y: roots (3|ℓ)·u/ℓ^S
    fb = nmod_poly([c % L for c in Gs], L); fdeg = sorted(ff.degree() for ff, e in fb.factor()[1] for _ in range(e))
    dig = lambda v: [(v // L ** i) % L for i in range(min(P, 4))]
    say('a group of %d root%s of slope %d (system of degree %s), carried through weights x = %s (each to ℓ^%d or more):' % (
        deg, '' if deg == 1 else 's', S, 1 if not isinstance(lab[0], tuple) else len(lab[0]) - 1, xs, MINP))
    say('   at weight −2, the polynomial of the constants (3|ℓ)·u/ℓ^%d: coefficients lowest first, each as digits lowest first: %s' % (
        S, [dig(c) for c in Gs]))
    say('   its factors on lane ℓ: degrees %s%s' % (fdeg, ';  the root: %s' % dig((-Gs[0]) % mP) if deg == 1 else ''))
say(); say('seconds: %.0f' % (time.time() - t0))
