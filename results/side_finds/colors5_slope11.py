# A dig, 12 Nov 2026: 5-COLORED PARTITIONS, THE LAW AT 11 FROM MODULAR FORMS (STATEMENT G3: "its family not yet computed").
# The exact 5-colour tower at 11 has the law 11^(⌊b/2⌋+1) (colors_tower_r5_*.json): two levels add one factor of 11, so
# by the reading that worked for 5, 7, 11 (partitions_flint_slopes.py) the tower follows a family of slope s = 1, with
# constant α = (3|ℓ)·u(−6)/ℓ (weight −(r + 1) = −6, twist (3|ℓ) for 5 colours).
# The class (5 colours: U_p = −(p|ℓ)·p^((k−2)/2) at weight ℓ − 7, the sign carried by (p|ℓ) per step of ℓ − 1) is empty at
# weights 4 and 14; at 24, 34, 44 it has exactly one form, of slope 1, with T_q numbers {5: 10, 7: 5, 13: 10} on lane 11;
# at 54 a slope-3 form joins (found 12 Nov). Method as partitions_flint_slopes.py: at k = −6 + 110·(x + 1) the piece with
# those T_q numbers, inside it the class (exact U₂, U₃), T_ℓ, its Newton polygon, the one root of slope 1, u from it,
# Newton's series at x = −1.
# PREDICTED before the run: exactly one root of slope 1 at every weight; the digits of (3|ℓ)·u(−6)/ℓ agree with the ratios
# H′_(b+2)/H′_b of the exact 5-colour tower (H′ = H/11^(⌊b/2⌋+1)) as far as those settle. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors5_slope11.py X N [--extra=E] [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat
import forms6_flint as F
flint.ctx.threads = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--threads=')), 8)
args = [int(a) for a in sys.argv[1:] if a.isdigit()]; L = 11; X, N = args[0], args[1]
EXTRA = next((int(a[8:]) for a in sys.argv[1:] if a.startswith('--extra=')), 30); NP = N + EXTRA
KEEP = '--keep' in sys.argv
OUT = os.path.join(HERE, 'colors5_slope11.txt')
if KEEP and os.path.exists(OUT): OUT = OUT[:-4] + '_N%d.txt' % N
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None
SETTLED = {11: []}                                             # filled from the exact tower below


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


S = 1; QS = [5, 7, 13]; t0 = time.time(); W0 = -6
sign = lambda p, k: -leg(p, L) ** (1 + (k - (L - 7)) // (L - 1))           # 5 colours: −(p|ℓ) at weight ℓ − 7, carried
say('5 colours, ℓ = %d, slope %d, weights −6 + %d·(x + 1) for x = 0..%d, modulo ℓ^%d (N = %d and %d spare); FLINT %s' % (
    L, S, (L - 1) * L, X, NP, N, EXTRA, flint.__version__))

# ---- the system on lane ℓ (found at weights 24, 34, 44) ----
ref = {5: 10, 7: 5, 13: 10}; how = 'the slope-1 form of the class at weights 24, 34, 44'
say('the system on lane ℓ: %s; T_q numbers %s' % (how, ref))


def compose(Y, R, Z, Rz):
    return Y * Z, [R[j] for j in Rz]


# ---- the members ----
ws = []; precs = []
for x in range(X + 1):
    k = W0 + (L - 1) * L * (x + 1); t = time.time(); m = L ** NP; n = k + 1
    ops = [('T', q) for q in QS] + [('T', L), ('U', 2), ('U', 3)]
    H = F.hecke_matrices(k, m, ops)
    Y, R, r = None, None, n                                    # the piece: the system's numbers, one T_q after another
    dims = []
    for q in QS:
        Tq = H[('T', q)] if Y is None else F.restrict(H[('T', q)], Y, R, r)
        Z, Rz, r2 = F.kernel_part(Tq - F.identity(r, m) * ref[q], L, NP)
        if r2 == 0: Y = None; r = 0; break
        if Y is None: Y, R = Z, Rz
        else: Y, R = compose(Y, R, Z, Rz)
        r = r2; dims.append(r)
    if r == 0: say('   weight %d: no form with the system\'s numbers; stopped' % k); break
    c2 = sign(2, k) * 2 ** ((k - 2) // 2) % m; c3 = sign(3, k) * 3 ** ((k - 2) // 2) % m
    U2 = F.restrict(H[('U', 2)], Y, R, r); Z, Rz, r2 = F.kernel_part(U2 - F.identity(r, m) * c2, L, NP)
    if r2 == 0: say('   weight %d: no form with the class\'s U₂ value; stopped' % k); break
    Y, R = compose(Y, R, Z, Rz); r = r2
    U3 = F.restrict(H[('U', 3)], Y, R, r); Z, Rz, r2 = F.kernel_part(U3 - F.identity(r, m) * c3, L, NP)
    if r2 == 0: say('   weight %d: no form with the class\'s U₃ value; stopped' % k); break
    Y, R = compose(Y, R, Z, Rz); r = r2
    U2W = F.restrict(H[('U', 2)], Y, R, r); U3W = F.restrict(H[('U', 3)], Y, R, r); TLW = F.restrict(H[('T', L)], Y, R, r)
    Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)] + \
         [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(r)] for i in range(r)]
    ker, N1 = F.padic_kernel(Bm, 2 * r, r, L, NP); d = len(ker); m1 = L ** N1
    Xn = [[ker[c][i] % m1 for c in range(d)] for i in range(r)]  # the new forms, in the piece's coordinates (r × d)
    Rn = F.pivot_rows(nmod_mat(r, d, [Xn[i][c] % L for i in range(r) for c in range(d)], L), L)
    XR = F.inverse_mod(F.mat_mod(d, d, [Xn[Rn[a]][b] for a in range(d) for b in range(d)], m1), L, N1)
    TX = F.mat_mod(r, r, [int(v) % m1 for v in TLW.entries()], m1) * F.mat_mod(r, d, [Xn[i][c] for i in range(r) for c in range(d)], m1)
    Mt = XR * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m1)
    chi = F.berkowitz([[int(Mt[i, j]) for j in range(d)] for i in range(d)], m1)
    lam, pr, slopes = F.slope_root(chi, L, S, N1)
    poly = ', '.join('%s×%d' % (sl, ln) for sl, ln in slopes)
    if lam is None:
        say('   weight %d: the piece %s, new at 2 and 3: %d (to ℓ^%d); slopes %s: no single root of slope %d; stopped' % (k, dims + [r], d, N1, poly, S)); break
    mp = L ** (pr - S); w = (lam // L ** S) % mp; e = k - 1 - 2 * S          # u = ℓ^S·w, w² − (λ/ℓ^S)·w + ℓ^e = 0
    for _ in range(pr + 2): w = ((lam // L ** S) - (L ** e % mp) * pow(w, -1, mp)) % mp if e < pr - S else (lam // L ** S) % mp
    ws.append(w); precs.append(pr - S)
    say('   weight %d (dimension %d): the piece %s, new at 2 and 3 with the class\'s signs: %d (to ℓ^%d); slopes of T_ℓ there: %s; '
        'u/ℓ^%d on lane ℓ: %d, known to ℓ^%d  (%.0f s)' % (k, n, dims + [r], d, N1, poly, S, w % L, pr - S, time.time() - t))

# ---- weight −6, and the exact tower ----
import json, glob
src = sorted(glob.glob(os.path.join(HERE, 'colors_tower_r5_N*_11_*.json')) + glob.glob(os.path.join(HERE, 'colors_tower_r5_N*11*.json')))
tower = None
for f_ in src:
    dj = json.load(open(f_, encoding='utf-8'))
    if '11' in dj['primes']: tower = dj; break
ratios = []
if tower:
    rec = tower['primes']['11']; Kt = rec['K']; Mt = L ** Kt; lv = {int(b): v for b, v in rec['levels'].items()}
    def val(z):
        z %= Mt
        if z == 0: return Kt
        v = 0
        while z % L == 0: z //= L; v += 1
        return v
    for b in sorted(lv):
        if b + 2 not in lv: continue
        e1, e2 = b // 2 + 1, (b + 2) // 2 + 1
        for n_ in range(min(len(lv[b]), len(lv[b + 2]))):
            h1, h2 = lv[b][n_] % Mt, lv[b + 2][n_] % Mt
            if val(h1) != e1: continue
            prec = Kt - e2
            if prec <= 0: continue
            a1 = (h1 // L ** e1) % L ** prec; a2 = (h2 // L ** e2) % L ** prec
            ratios.append((b, n_, a2 * pow(a1, -1, L ** prec) % L ** prec, prec))
if len(ws) == X + 1:
    P = min(precs); v, val_, _ = F.newton_at_minus_one([w % L ** P for w in ws], L, P)
    a = leg(3, L) * v % L ** P; K = min(P, X + 1)
    say('the differences are divisible by ℓ to %s' % val_)
    dg = [(a // L ** j) % L for j in range(P)]
    say("(3|ℓ)·u(−6)/ℓ, digits lowest first (the first %d are the series' reach, x + 1 = %d weights): %s" % (K, X + 1, dg[:K + 2]))
    if ratios:
        for b_ in sorted(set(r_[0] for r_ in ratios)):
            rs = [r_ for r_ in ratios if r_[0] == b_]
            agree = [next((j for j in range(min(K, r_[3])) if (r_[2] // L ** j) % L != dg[j]), min(K, r_[3])) for r_ in rs]
            say("   exact tower, H′_%d/H′_%d over n = %s: digits agreeing with the forms' constant %s (of at most %s)" % (
                b_ + 2, b_, [r_[1] for r_ in rs], agree, [min(K, r_[3]) for r_ in rs]))
    else:
        say('   no exact 5-colour tower at 11 found to compare')
say(); say('seconds: %.0f' % (time.time() - t0))
