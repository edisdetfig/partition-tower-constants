# A dig, 12 Nov 2026: r-COLORED PARTITIONS, THE TOWERS FROM EXACT VALUES (any odd r; colors3_tower.py was r = 3).
# p_r(n) is the coefficient of qⁿ in ∏(1 − qⁿ)^(−r) = 1/E_r(q); E_r = J·E₁^(r−3) for r ≥ 3, with J = ∏(1 − qⁿ)³ (Jacobi's
# sparse series) and E₁ = ∏(1 − qⁿ) (Euler's sparse pentagonal series), else E₁ itself. 1/E_r is computed on lane ℓ^K
# (ℓ^K below 2^BITS) to N terms with FLINT, and the tower H_b(n) = p_r(ℓᵇn + δ_b), 24δ_b ≡ r (mod ℓᵇ), is read off.
# WHERE IT RUNS: in WSL with python-flint 0.9.0 from its manylinux wheel (the Windows build crashes on products of more
# than about 2·10⁷ terms; colors3_tower.py). The first values are printed as a check against the known p_r(0), p_r(1), ….
# No forms are used here. COUNTED.
#   usage (WSL): python3 -I results/side_finds/colors_tower.py <path of the unpacked wheel> r N ℓ … [--bits=30] [--keep] [--out=<dir>]
import os, sys, time, json
sys.path.insert(0, sys.argv[1])
import flint
from flint import nmod_poly
flint.ctx.threads = 8
HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = next((a[6:] for a in sys.argv[1:] if a.startswith('--out=')), HERE)
r = int(sys.argv[2]); N = int(sys.argv[3]); PS = [int(a) for a in sys.argv[4:] if a.isdigit()]; KEEP = '--keep' in sys.argv
BITS = next((int(a[7:]) for a in sys.argv[1:] if a.startswith('--bits=')), 30)
tag = 'r%d_N%d_%s' % (r, N, '_'.join(map(str, PS)))
fh = open(os.path.join(OUTD, 'colors_tower_%s.txt' % tag), 'x', encoding='utf-8') if KEEP else None
out = {'r': r, 'N': N, 'primes': {}}


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def sparse(M, kind):
    P = nmod_poly([], M); terms = []
    if kind == 'J':                                            # Σ (−1)^j (2j + 1) q^(j(j+1)/2)
        j = 0
        while j * (j + 1) // 2 < N: terms.append((j * (j + 1) // 2, (-1) ** j * (2 * j + 1))); j += 1
    else:                                                      # Σ (−1)^j q^(j(3j−1)/2), j ∈ Z
        terms.append((0, 1)); j = 1
        while j * (3 * j - 1) // 2 < N:
            terms.append((j * (3 * j - 1) // 2, (-1) ** j))
            if j * (3 * j + 1) // 2 < N: terms.append((j * (3 * j + 1) // 2, (-1) ** j))
            j += 1
        terms.sort()
    for e, c in reversed(terms): P[e] = c % M
    return P


t0 = time.time(); say('p_%d(n) on lane ℓ^K for n < %d (FLINT %s, WSL)' % (r, N, flint.__version__))
for L in PS:
    t = time.time(); K = 1
    while L ** (K + 1) < 2 ** BITS: K += 1
    M = L ** K
    if r >= 3:
        E = sparse(M, 'J')
        if r > 3:
            E1 = sparse(M, 'E'); E = E.mul_low(E1.pow_trunc(r - 3, N), N); del E1
    else:
        E = sparse(M, 'E').pow_trunc(r, N)
    P = E.inverse_series_trunc(N); del E
    first = [int(P[i]) for i in range(8)]
    lv = {}
    for b in range(1, 40):
        if L ** b >= N: break
        db = r * pow(24, -1, L ** b) % L ** b
        vals = [int(P[L ** b * n + db]) for n in range(12) if L ** b * n + db < N]
        if not vals: break
        lv[b] = vals
    del P
    out['primes'][str(L)] = {'K': K, 'levels': {str(b): v for b, v in lv.items()}}
    say(); say('ℓ = %d, lane ℓ^%d (%.0f s): p_%d(0..7) on that lane %s; levels %s' % (L, K, time.time() - t, r, first, {b: len(v) for b, v in lv.items()}))
    def val(x):
        if x == 0: return K
        v = 0
        while x % L == 0: x //= L; v += 1
        return v
    for b, v in lv.items():
        say('   b = %d: ℓ-adic valuations %s;  on lane ℓ %s' % (b, [val(x) for x in v], [x % L for x in v]))
if KEEP:
    json.dump(out, open(os.path.join(OUTD, 'colors_tower_%s.json' % tag), 'x', encoding='utf-8'))
say(); say('seconds: %.0f' % (time.time() - t0))
