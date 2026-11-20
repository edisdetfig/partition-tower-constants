# A dig, 12 Nov 2026: TOWERS OF LEVEL-1 FORMS, THE CONSTANTS FROM HIDA FAMILIES (FLINT only; forms1_flint.py).
# A weakly holomorphic form F of weight w and level 1 (1/Δ: w = −12, the 24-colored partitions; j: w = 0; 1/E₄: w = −4)
# has the tower F | U_ℓᵇ = Σ c(ℓᵇn) qⁿ. By the classical route of PROOF.md (no Shimura step at integral weight), its ℓ-adic
# constants, one level apart, should be the U_ℓ-values at weight w of the ordinary families of level 1 through the weights
# k ≡ w (mod ℓ − 1): EXACTLY 1 for the Eisenstein family (U_ℓ = 1 at every weight), and A_f(w) for each ordinary cusp form.
# Method: at k = k₀ + (ℓ − 1)·x (k₀ the least such weight ≥ 4), the ordinary part of T_ℓ on M_k(SL₂(ℤ)) modulo ℓ^N, the
# unit roots (the ordinary ℓ-stabilizations), their characteristic polynomial; Newton's series at x₀ = (w − k₀)/(ℓ − 1).
# PREDICTED before the run: the polynomial has the factor x − 1 exactly at every weight and at w; its other roots are the
# cusp forms' families. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/forms1_family.py ℓ w X N [--keep] [--threads=t]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly
import forms1_flint as F1
import forms6_flint as F
flint.ctx.threads = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--threads=')), 8)
nums = [a for a in sys.argv[1:] if a.lstrip('-').isdigit()]; L, W, X, N = (int(a) for a in nums[:4])
KEEP = '--keep' in sys.argv
OUT = os.path.join(HERE, 'forms1_family_w%d_%d.txt' % (W, L))
if KEEP and os.path.exists(OUT): OUT = OUT[:-4] + '_N%d.txt' % N
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); m = L ** N
t_ = 0
while W + (L - 1) * t_ < 4: t_ += 1
k0 = W + (L - 1) * t_; x0 = -t_
say('ℓ = %d, a level-1 form of weight %d: families through k₀ = %d, read at x = %d (weight %d); %d weights; modulo ℓ^%d; FLINT %s' % (
    L, W, k0, x0, W, X + 1, N, flint.__version__))
series = []; dords = []
for x in range(X + 1):
    k = k0 + (L - 1) * x; t = time.time()
    T = F1.hecke_matrices1(k, m, [L])[L]
    Y, R, r = F.ordinary_part(T, L, N)
    Ml = F.restrict(T, Y, R, r)
    U = F.unit_root_matrix(Ml, r, L, N, k); cp = F.charpoly_small(U, r, m)
    series.append(cp); dords.append(r)
    say('   weight %d (dimension %d): the ordinary part has %d forms; U_ℓ polynomial on lane ℓ %s  (%.0f s)' % (
        k, F1.dim_mk(k), r, [c % L for c in cp], time.time() - t))
if len(set(dords)) != 1: say('   the ordinary part changes size along the weights; stopped'); sys.exit()
d = dords[0]; poly = []; vals = []
for j in range(d):
    v, val = F1.newton_at([series[x][j] for x in range(X + 1)], L, N, x0)
    poly.append(v); vals.append(val)
poly.append(1)
K = min(N, X + 1); MK = L ** K
say('the differences of each coefficient are divisible by ℓ to %s' % vals)
say('the constants\' polynomial at weight %d (lowest first), on lane ℓ^%d: %s;  on lane ℓ: %s' % (W, K, [c % MK for c in poly], [c % L for c in poly]))
one = sum(c for c in poly) % MK
say('   x = 1 is a root to ℓ^%d (the Eisenstein family): %s' % (K, one == 0))
fac = nmod_poly([c % L for c in poly], L).factor()[1]
say('   on lane ℓ it factors as %s' % [([int(c) for c in f.coeffs()], e) for f, e in fac])
if d >= 2 and one == 0:                                                    # divide out x − 1: the cusp forms' polynomial
    q = [0] * d; acc = 0
    for i in range(d, 0, -1):                                             # synthetic division by (x − 1), highest first
        acc = (poly[i] + acc) % MK; q[i - 1] = acc
    say('   the cusp forms\' polynomial (lowest first), on lane ℓ^%d: %s' % (K, q))
    if d == 2:
        a = (-q[0]) % MK
        say('   the cusp form\'s constant, digits lowest first: %s' % [(a // L ** e) % L for e in range(K)])
say(); say('seconds: %.0f' % (time.time() - t0))
