# A dig, 12 Nov 2026: 3-COLORED PARTITIONS, THE TOWERS FROM EXACT VALUES. p₃(n) is the coefficient of qⁿ in ∏(1 − qⁿ)^(−3)
# = 1/J(q), J(q) = ∏(1 − qⁿ)³ = Σ_j (−1)^j (2j + 1) q^(j(j+1)/2) (Jacobi). Here 1/J is computed modulo ℓ^K to N terms
# (FLINT's power-series inverse; exact on each lane), and the tower H_b(n) = p₃(ℓᵇn + δ_b), 24δ_b ≡ 3 (mod ℓᵇ) — Belmont,
# Lee, Musat and Trebat-Leder's p₃((ℓᵇN + 3)/24) — is read off for every level b and n that fit below N.
# WHERE IT RUNS: in WSL (Linux) with python-flint 0.9.0 unpacked from its manylinux wheel. The Windows build of the same
# FLINT crashes on products of more than about 2·10⁷ terms (a 32-bit size overflow), counted 12 Nov 2026.
# No forms are used here; colors3_compare.py compares. COUNTED.
#   usage (WSL): python3 -I results/side_finds/colors3_tower.py <path of the unpacked wheel> N ℓ … [--bits=62] [--keep] [--out=<dir>]
# A first run at N = 3·10⁸ (60-bit lanes) died without output: FLINT's products at that length need more than WSL's 31 GB
# (counted: 10⁸ terms peak at 13 GB with 30-bit lanes, 20 GB with 60-bit). Its file colors3_tower_13_…_47.txt holds only
# the header.
import os, sys, time, json
sys.path.insert(0, sys.argv[1])
import flint
from flint import nmod_poly
flint.ctx.threads = 8
HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = next((a[6:] for a in sys.argv[1:] if a.startswith('--out=')), HERE)   # a trial writes elsewhere
N = int(sys.argv[2]); PS = [int(a) for a in sys.argv[3:] if a.isdigit()]; KEEP = '--keep' in sys.argv
BITS = next((int(a[7:]) for a in sys.argv[1:] if a.startswith('--bits=')), 62)   # the lane ℓ^K stays below 2^BITS
tag = 'N%d_%s' % (N, '_'.join(map(str, PS)))
fh = open(os.path.join(OUTD, 'colors3_tower_%s.txt' % tag), 'x', encoding='utf-8') if KEEP else None
out = {'N': N, 'primes': {}}


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); say('p₃(n) modulo ℓ^K for n < %d (FLINT %s, WSL)' % (N, flint.__version__))
for L in PS:
    t = time.time(); K = 1
    while L ** (K + 1) < 2 ** BITS: K += 1
    M = L ** K; J = nmod_poly([], M); terms = []; j = 0
    while j * (j + 1) // 2 < N:
        terms.append((j * (j + 1) // 2, (-1) ** j * (2 * j + 1))); j += 1
    for e, c in reversed(terms): J[e] = c % M
    P = J.inverse_series_trunc(N); del J
    lv = {}
    for b in range(1, 40):
        if L ** b >= N: break
        db = 3 * pow(24, -1, L ** b) % L ** b
        vals = [int(P[L ** b * n + db]) for n in range(12) if L ** b * n + db < N]
        if not vals: break
        lv[b] = vals
    del P
    out['primes'][str(L)] = {'K': K, 'levels': {str(b): v for b, v in lv.items()}}
    say(); say('ℓ = %d, modulo ℓ^%d (%.0f s): levels %s' % (L, K, time.time() - t, {b: len(v) for b, v in lv.items()}))
    def val(x):
        if x == 0: return K
        v = 0
        while x % L == 0: x //= L; v += 1
        return v
    for b, v in lv.items():
        say('   b = %d: ℓ-adic valuations %s;  on lane ℓ %s' % (b, [val(x) for x in v], [x % L for x in v]))
if KEEP:
    json.dump(out, open(os.path.join(OUTD, 'colors3_tower_%s.json' % tag), 'x', encoding='utf-8'))
say(); say('seconds: %.0f' % (time.time() - t0))
