# A dig, 12 Nov 2026: 3-COLORED PARTITIONS, THE LAWS AT 11 AND 17 FROM MODULAR FORMS. The exact towers
# (colors3_tower_N150000000_*.txt) show p₃(11ᵇn + δ_b) divisible by 11ᵇ and p₃(17ᵇn + δ_b) by 17^⌈b/2⌉ (R_ℓ(3) = 0 there):
# two levels add ℓ² and ℓ, so by the reading that worked for 5, 7, 11 (partitions_flint_slopes.py) the towers should follow
# families of slope s = 2 (11) and 1 (17), with constant α = (−1|ℓ)·u(−4)/ℓ^s (the twist of 3 colours).
# Method (FLINT only; forms2_flint.py, level 2): at the weights k = −4 + (ℓ − 1)·ℓ·(x + 1), the class (U₂ = s₂(k)·2^((k−2)/2),
# s₂(k) = (2|ℓ)^(1 + (k − ℓ + 5)/(ℓ − 1)) as the sign rule carries it) lifted over the ℓ-adic integers and separated exactly,
# T_ℓ on it, its characteristic polynomial (Berkowitz), its Newton polygon, the one root of slope s (Hensel), u from it,
# Newton's series at x = −1.
# PREDICTED before the run: exactly one root of slope s at every weight; the digits of (−1|ℓ)·u(−4)/ℓ^s agree with the
# ratios H′_(b+2)/H′_b of the exact towers (H′ = H divided by the law's power) as far as those settle. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors3_slopes.py ℓ X N [--extra=E] [--pick=a] [--keep] [--threads=t]
# For 11 two forms of slope 2 share the class at weight 106: --pick=a follows the one with u/ℓ² ≡ a on lane ℓ (both are
# followed, and the exact towers decide).
import os, sys, time, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat
import forms2_flint as F2
import forms6_flint as F
flint.ctx.threads = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--threads=')), 8)
args = [int(a) for a in sys.argv[1:] if a.isdigit()]; L, X, N = args[0], args[1], args[2]
EXTRA = next((int(a[8:]) for a in sys.argv[1:] if a.startswith('--extra=')), 40); NP = N + EXTRA
KEEP = '--keep' in sys.argv; S = {11: 2, 17: 1}[L]
PICK = next((int(a[7:]) for a in sys.argv[1:] if a.startswith('--pick=')), None)   # u/ℓ^s on lane ℓ, when two roots share slope s
OUT = os.path.join(HERE, 'colors3_slopes_%d%s.txt' % (L, '_pick%d' % PICK if PICK is not None else ''))
if KEEP and os.path.exists(OUT): OUT = OUT[:-4] + '_N%d.txt' % N
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


t0 = time.time(); m = L ** NP
say('ℓ = %d (3 colours), slope %d, weights −4 + %d·(x + 1), x = 0..%d, modulo ℓ^%d (N = %d and %d spare); FLINT %s' % (
    L, S, (L - 1) * L, X, NP, N, EXTRA, flint.__version__))
ws = []; precs = []; pairs = []
for x in range(X + 1):
    k = -4 + (L - 1) * L * (x + 1); t = time.time(); n = k // 4 + 1
    s2 = leg(2, L) ** (1 + (k - (L - 5)) // (L - 1)); c2 = s2 * 2 ** ((k - 2) // 2) % m
    H = F2.hecke_matrices2(k, m, [('U', 2), ('T', L)])
    Y, R, r = F.kernel_part(H[('U', 2)] - F.identity(n, m) * c2, L, NP)
    if r == 0: say('   weight %d: the class is empty; stopped' % k); break
    U2W = F.restrict(H[('U', 2)], Y, R, r); TLW = F.restrict(H[('T', L)], Y, R, r)
    ker, N1 = F.padic_kernel([[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)], r, r, L, NP)
    d = len(ker); m1 = L ** N1
    Xn = [[ker[c][i] % m1 for c in range(d)] for i in range(r)]
    Rn = F.pivot_rows(nmod_mat(r, d, [Xn[i][c] % L for i in range(r) for c in range(d)], L), L)
    XR = F.inverse_mod(F.mat_mod(d, d, [Xn[Rn[a]][b] for a in range(d) for b in range(d)], m1), L, N1)
    TX = F.mat_mod(r, r, [int(v) % m1 for v in TLW.entries()], m1) * F.mat_mod(r, d, [Xn[i][c] for i in range(r) for c in range(d)], m1)
    Mt = XR * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m1)
    chi = F.berkowitz([[int(Mt[i, j]) for j in range(d)] for i in range(d)], m1)
    lam, pr, slopes = F.slope_root(chi, L, S, N1, pick=PICK)
    poly = ', '.join('%s×%d' % (sl, ln) for sl, ln in slopes)
    if lam is None and PICK is None:                            # several roots of slope s: follow their polynomial
        G, prg = F.slope_factor(chi, L, S, N1)
        if G is None: say('   weight %d: the class %d forms (to ℓ^%d); slopes %s: no factor of slope %d; stopped' % (k, d, N1, poly, S)); break
        pairs.append(G); precs.append(prg)
        say('   weight %d (dimension %d): the class %d forms (to ℓ^%d); slopes of T_ℓ %s; the %d roots of slope %d: the polynomial of '
            'root/ℓ^%d on lane ℓ %s, known to ℓ^%d  (%.0f s)' % (k, n, d, N1, poly, len(G) - 1, S, S, [c % L for c in G], prg, time.time() - t))
        continue
    if lam is None: say('   weight %d: the class %d forms (to ℓ^%d); slopes %s: no single root of slope %d; stopped' % (k, d, N1, poly, S)); break
    mp = L ** (pr - S); e = k - 1 - 2 * S; w = (lam // L ** S) % mp
    for _ in range(pr + 2): w = ((lam // L ** S) - (L ** e % mp) * pow(w, -1, mp)) % mp if e < pr - S else (lam // L ** S) % mp
    ws.append(w); precs.append(pr - S)
    say('   weight %d (dimension %d): the class %d forms (to ℓ^%d); slopes of T_ℓ %s; u/ℓ^%d on lane ℓ %d, known to ℓ^%d  (%.0f s)' % (
        k, n, d, N1, poly, S, w % L, pr - S, time.time() - t))
if len(pairs) == X + 1:                                       # the polynomial of the roots at weight −4, twisted by (−1|ℓ)
    P = min(precs); mP = L ** P; deg = len(pairs[0]) - 1; t_ = leg(-1, L); co = []
    for j in range(deg):
        v, val, _ = F.newton_at_minus_one([G[j] % mP for G in pairs], L, P)
        co.append(v * t_ ** (deg - j) % mP); say('the differences of coefficient %d are divisible by ℓ to %s' % (j, val))
    co.append(1); K = min(P, X + 1)
    say("the constants' polynomial at weight −4, twisted by (−1|ℓ) = %+d (lowest first), on lane ℓ^%d: %s;  on lane ℓ: %s" % (
        t_, K, [c % L ** K for c in co], [c % L for c in co]))
    from flint import nmod_poly
    say('   on lane ℓ it factors as %s' % [([int(x) for x in f.coeffs()], e) for f, e in nmod_poly([c % L for c in co], L).factor()[1]])
if len(ws) == X + 1:
    P = min(precs); v, val, _ = F.newton_at_minus_one([w % L ** P for w in ws], L, P)
    a = leg(-1, L) * v % L ** P; K = min(P, X + 1)
    say('the differences are divisible by ℓ to %s' % val)
    say('(−1|ℓ)·u(−4)/ℓ^%d, digits lowest first (the series reaches %d): %s' % (S, K, [(a // L ** j) % L for j in range(min(P, K + 3))]))
say(); say('seconds: %.0f' % (time.time() - t0))
