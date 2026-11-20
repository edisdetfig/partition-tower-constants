# A dig, 12 Nov 2026: 1/E₄, THE CONSTANT OF THE POLE. 1/E₄ (weight −4) has a pole at ρ = e^(2πi/3), the point of the curve
# y² = x³ + 1 (j = 0, complex multiplication by ℤ[ζ₃]). At ℓ ≡ 2 mod 3 that curve is supersingular, and the tower follows
# the forms' constants alone. At ℓ ≡ 1 mod 3 it is ordinary: its Frobenius has a unit root u (the root of x² − aₗx + ℓ that
# is an ℓ-adic unit, 4ℓ = aₗ² + 3b²; u⁶ is the same for all six twists of the curve), and a small first run (n < 3·10⁵,
# at 7 and 13) showed the tower needs one constant more than the forms give, equal there to u⁻⁶ (9 digits at 7, 3 at 13;
# a search over u^m, m = −60…60, both signs).
# PREDICTED before this run (derived from that first run, and from the weight: u^(w−2) with w = −4): at every ℓ ≡ 1 mod 3
# the extra constant is u⁻⁶, and the law P(x)·(x − u⁻⁶) holds, deeper level by level. Tested at 19, 31, 37, 43 (not used
# to find it) and again at 7 and 13 with deeper data.
# Per prime: (1) the extra constant read off the data alone: the remainders G_b(n) = Σ_i [xⁱ]P·H_(b+i)(n) of the forms'
# law, one level apart, G_(b+1)/G_b at the last level, agreed over n; (2) the power u^m (m = −60…60, ±) that equals it to
# the most digits; (3) the depth of the law P(x)·(x − u^m) from every level, for n = 0…11. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/level1_pole_check.py <1/E₄ tower json> [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRC = next(a for a in sys.argv[1:] if a.endswith('.json'))
fh = open(os.path.join(HERE, 'level1_pole_check_%s.txt' % os.path.basename(SRC)[13:-5]), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def unit_root(L, K):
    """the unit root u of x² − a x + L on lane L^K, 4L = a² + 3b² (the curve y² = x³ + 1 up to twist), and a."""
    M = L ** K
    a = next(a for a in range(1, 2 * L) if 4 * L - a * a > 0 and (4 * L - a * a) % 3 == 0
             and round(((4 * L - a * a) // 3) ** 0.5) ** 2 == (4 * L - a * a) // 3)
    u = a % L
    for _ in range(4 * K + 8): u = (u - (u * u - a * u + L) * pow((2 * u - a) % M, -1, M)) % M
    return u, a


data = json.load(open(os.path.join(HERE, SRC), encoding='utf-8')); assert data['form'] == 'invE4'
say('1/E₄ towers: %s (n < %d)' % (SRC, data['N']))
for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
    L = int(Ls); f = os.path.join(HERE, 'forms1_family_w-4_%d.txt' % L)
    mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(f, encoding='utf-8').read())
    KF = int(mt.group(1)); P = [int(v) for v in mt.group(2).split(',')]; D = len(P) - 1
    K = min(KF, rec['K']); M = L ** K; lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}; top = max(lv)
    say(); say('ℓ = %d (%s), lane ℓ^%d; the forms\' polynomial P on lane ℓ %s' % (
        L, '≡ 1 mod 3: ρ ordinary' if L % 3 == 1 else '≡ 2 mod 3: ρ supersingular', K, [c % L for c in P]))
    G = {}
    for b in sorted(lv):
        if b + D <= top:
            nn = min(len(lv[b + i]) for i in range(D + 1))
            G[b] = [sum(P[i] * lv[b + i][n] for i in range(D + 1)) % M for n in range(nn)]
    say('   P\'s law alone, depth from levels %s: %s' % (sorted(G), [min(val(g, L, K) for g in G[b][1:]) for b in sorted(G)]))
    if L % 3 != 1: continue
    pairs = [b for b in G if b + 1 in G]                                  # (1) the extra constant from the last level
    bl = max(pairs) if pairs else None
    ests = []
    for n in (range(1, min(len(G[bl]), len(G[bl + 1]))) if bl is not None else []):
        v = val(G[bl][n], L, K)
        if v >= K - 1 or G[bl + 1][n] % L ** v: continue
        mk = L ** (K - v); ests.append(((G[bl + 1][n] // L ** v) * pow((G[bl][n] // L ** v) % mk, -1, mk) % mk, K - v))
    u, a = unit_root(L, K)
    if not ests: say('   too few levels to read the constant from the remainders')
    else:
      kk = min(k for _, k in ests); lam = ests[0][0] % L ** kk
      agree = next((e for e in range(kk) if any((x - lam) % L ** (e + 1) for x, _ in ests)), kk)
      say('   the extra constant from the data, G_%d/G_%d over %d values of n: digits lowest first %s (they agree to %d digits)' % (
        bl + 1, bl, len(ests), [(lam // L ** e) % L for e in range(agree)], agree))
      hits = []                                                           # (2) which power of u it is
      for m in range(-60, 61):
          for sg in (1, -1):
              w = sg * pow(u, m, M) % M
              hits.append((next((e for e in range(agree) if (w - lam) % L ** (e + 1)), agree), sg, m))
      hits.sort(key=lambda h: (-h[0], abs(h[2])))
      say('   the best powers ±u^m (digits agreeing, sign, m): %s' % hits[:4])
    say('   y² = x³ + 1 at ℓ: aₗ = ±%d up to twist; u⁻⁶ digits %s' % (a, [(pow(u, -6, M) // L ** e) % L for e in range(K)]))
    lamM = pow(u, -6, M)                                                  # (3) the law P(x)·(x − u⁻⁶)
    Q = [0] * (D + 2)
    for i, c in enumerate(P): Q[i + 1] = (Q[i + 1] + c) % M; Q[i] = (Q[i] - lamM * c) % M
    for b in sorted(lv):
        if b + D + 1 > top: continue
        nn = min(len(lv[b + i]) for i in range(D + 2))
        deps = [val(sum(Q[i] * lv[b + i][n] for i in range(D + 2)), L, K) for n in range(nn)]
        say('   P(x)·(x − u⁻⁶) from level %d: depth for n = 0, 1, …  %s' % (b, deps))
