# A dig, 12 Nov 2026: ANDREWS' SMALLEST-PARTS FUNCTION spt(n), ITS TOWERS. spt(n) counts the smallest parts over all
# partitions of n (1, 3, 5, 10, 14, 26, …). Andrews' identity (checked against brute force for n ≤ 20):
#     Σ spt(n) qⁿ = Σ n·p(n) qⁿ + P(q)·B(q),   P = 1/∏(1 − qⁿ),   B = Σ_{n≥1} (−1)ⁿ q^(n(3n+1)/2) (1 + qⁿ)/(1 − qⁿ)².
# spt's modular object has weight 3/2 (12·spt(n) + (24n − 1)·p(n) are the coefficients of a mock modular form whose shadow
# is η(24z)), so by the recipe its towers H_b(n) = spt(ℓᵇn + δ_b), 24δ_b ≡ 1 (mod ℓᵇ), follow families read at the
# Shimura weight 2: a classical weight, where an ordinary family of tame level 6 meets a newform of level 6ℓ that is new at
# ℓ, whose U_ℓ-number is exactly +1 or −1 (−w_ℓ: split or non-split reduction for an elliptic curve).
# PREDICTED before the run (derived):
#   - every constant is ±1, so H_(b+4) ≡ H_b on ever deeper lanes ℓᵉ;
#   - on lane ℓ, H_(b+2) ≡ (3|ℓ)·H_b (the twist (3|ℓ) times +1; Ono's congruence for spt has this first digit, from memory
#     of PNAS 2011, to be checked), and deeper too where every family has U_ℓ = +1 at weight 2;
#   - laws at 5, 7 and 13 (Garvan, Trans. AMS 2012, from memory: spt(ℓᵃn + δ_a) ≡ 0 mod ℓ^⌊(a+1)/2⌋).
# Exact values on lane ℓ^K (ℓ^K < 2^30): B built with numpy (strided sums), the products with FLINT (Windows; n < N ≤ 1.6·10⁷).
# COUNTED.  usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/spt_tower.py N ℓ … [--keep]
import os, sys, time, json
import numpy as np
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly
flint.ctx.threads = 8
HERE = os.path.dirname(os.path.abspath(__file__))
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; N, PS = nums[0], nums[1:]; KEEP = '--keep' in sys.argv
tag = 'N%d_%s' % (N, '_'.join(map(str, PS)))
fh = open(os.path.join(HERE, 'spt_tower_%s.txt' % tag), 'x', encoding='utf-8') if KEEP else None
out = {'N': N, 'primes': {}}


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


t0 = time.time()
Bz = np.zeros(N, dtype=np.int64); n = 1                        # B over the integers (the sums stay far below 2^63)
while n * (3 * n + 1) // 2 < N:
    st = n * (3 * n + 1) // 2; cnt = (N - 1 - st) // n + 1
    Bz[st::n][:cnt] += (1 if n % 2 == 0 else -1) * (2 * np.arange(cnt, dtype=np.int64) + 1)
    n += 1
say('spt(n) on lane ℓ^K for n < %d; B built in %.0f s' % (N, time.time() - t0))
for L in PS:
    t = time.time(); K = 1
    while L ** (K + 1) < 2 ** 30: K += 1
    M = L ** K
    E = nmod_poly([], M); terms = [(0, 1)]; j = 1                  # Euler's pentagonal series, inverted: P
    while j * (3 * j - 1) // 2 < N:
        terms.append((j * (3 * j - 1) // 2, (-1) ** j))
        if j * (3 * j + 1) // 2 < N: terms.append((j * (3 * j + 1) // 2, (-1) ** j))
        j += 1
    for e, c in sorted(terms, reverse=True): E[e] = c % M
    P = E.inverse_series_trunc(N); del E
    B = nmod_poly((Bz % M).tolist(), M)
    PB = P.mul_low(B, N); del B
    spt = lambda m: (m * int(P[m]) + int(PB[m])) % M
    first = [spt(m) for m in range(1, 11)]
    lv = {}
    for b in range(1, 40):
        if L ** b >= N: break
        db = pow(24, -1, L ** b)
        vals = [spt(L ** b * nn + db) for nn in range(12) if L ** b * nn + db < N]
        if not vals: break
        lv[b] = vals
    del P, PB
    out['primes'][str(L)] = {'K': K, 'levels': {str(b): v for b, v in lv.items()}}
    def val(x):
        if x == 0: return K
        v = 0
        while x % L == 0: x //= L; v += 1
        return v
    say(); say('ℓ = %d, lane ℓ^%d (%.0f s): spt(1..10) on that lane %s; levels %s; (3|ℓ) = %+d' % (
        L, K, time.time() - t, first, {b: len(v) for b, v in lv.items()}, leg(3, L)))
    for b, v in lv.items():
        say('   b = %d: ℓ-adic valuations %s;  on lane ℓ %s' % (b, [val(x) for x in v], [x % L for x in v]))
    s3 = leg(3, L)
    for b in lv:                                                   # the predicted laws, measured
        if b + 2 in lv:
            nn = min(len(lv[b]), len(lv[b + 2]))
            d2 = [val((lv[b + 2][i] - s3 * lv[b][i]) % M) for i in range(nn)]
            say('   H_%d − (3|ℓ)·H_%d: ℓ-adic depth %s' % (b + 2, b, d2))
        if b + 4 in lv:
            nn = min(len(lv[b]), len(lv[b + 4]))
            d4 = [val((lv[b + 4][i] - lv[b][i]) % M) for i in range(nn)]
            say('   H_%d − H_%d: ℓ-adic depth %s' % (b + 4, b, d4))
if KEEP:
    json.dump(out, open(os.path.join(HERE, 'spt_tower_%s.json' % tag), 'x', encoding='utf-8'))
say(); say('seconds: %.0f' % (time.time() - t0))
