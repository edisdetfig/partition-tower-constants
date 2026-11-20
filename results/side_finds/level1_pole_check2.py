# A dig, 12 Nov 2026: THE POLE RULE ON A SECOND CURVE. level1_pole_check.py found, for 1/E₄ (weight w = −4, pole at ρ, the
# curve y² = x³ + 1), that at primes where that curve is ordinary the tower needs one constant beyond the forms' P: u^(w−2) =
# u⁻⁶, u the unit root of the curve's Frobenius. The reading (derived): T_ℓ multiplies the residue at a CM pole by
# π^(w−2) + π̄^(w−2), and only the unit root can stay in a tower of integers.
# Here the same rule is tested where it was not found: 1/E₆ (weight −6), whose pole is at i, the point of y² = x³ − x
# (complex multiplication by ℤ[i]). PREDICTED before the 1/E₆ run (derived from the rule):
#   ℓ ≡ 3 mod 4 (the curve supersingular): the forms' P alone, deeper level by level;
#   ℓ ≡ 1 mod 4 (ordinary): the law P(x)·(x − u⁻⁸), u the unit root of x² − aₗx + ℓ with aₗ = 2a, ℓ = a² + b²; u⁴ is the same
#   for all four twists of the curve, so u⁻⁸ is canonical (and 4 divides w − 2 = −8, as a residue at i needs).
# Per prime: P's law alone; at the ordinary primes the extra constant read off the data alone (the remainders of P's law one
# level apart), the power ±u^m (m = −60…60) it equals to the most digits, and the depth of P(x)·(x − u^(w−2)). The same
# script reads 1/E₄ towers (curve y² = x³ + 1, u^(w−2) = u⁻⁶). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/level1_pole_check2.py <tower json> [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRC = next(a for a in sys.argv[1:] if a.endswith('.json'))
fh = open(os.path.join(HERE, 'level1_pole_check2_%s.txt' % os.path.basename(SRC)[13:-5]), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def unit_root(L, K, form):
    """the unit root u on lane L^K of the CM curve at the pole: y² = x³ + 1 (1/E₄) or y² = x³ − x (1/E₆)."""
    M = L ** K
    if form == 'invE4':
        a = next(a for a in range(1, 2 * L) if 4 * L - a * a > 0 and (4 * L - a * a) % 3 == 0
                 and round(((4 * L - a * a) // 3) ** 0.5) ** 2 == (4 * L - a * a) // 3)
    else:
        a = 2 * next(a for a in range(1, L) if round((L - a * a) ** 0.5) ** 2 == L - a * a and L - a * a > 0)
    u = a % L
    for _ in range(4 * K + 8): u = (u - (u * u - a * u + L) * pow((2 * u - a) % M, -1, M)) % M
    return u, a


data = json.load(open(os.path.join(HERE, SRC), encoding='utf-8')); FORM = data['form']
W, curve, split = {'invE4': (-4, 'y² = x³ + 1', lambda L: L % 3 == 1), 'invE6': (-6, 'y² = x³ − x', lambda L: L % 4 == 1)}[FORM]
E = W - 2
say('%s towers: %s (n < %d); the pole\'s curve %s; the rule\'s constant u^(w−2) = u^%d' % (FORM, SRC, data['N'], curve, E))
for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
    L = int(Ls); f = os.path.join(HERE, 'forms1_family_w%d_%d.txt' % (W, L))
    mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(f, encoding='utf-8').read())
    KF = int(mt.group(1)); P = [int(v) for v in mt.group(2).split(',')]; D = len(P) - 1
    K = min(KF, rec['K']); M = L ** K; lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}; top = max(lv)
    say(); say('ℓ = %d (%s), lane ℓ^%d; the forms\' polynomial P on lane ℓ %s' % (
        L, 'the curve ordinary' if split(L) else 'the curve supersingular', K, [c % L for c in P]))
    G = {}
    for b in sorted(lv):
        if b + D <= top:
            nn = min(len(lv[b + i]) for i in range(D + 1))
            G[b] = [sum(P[i] * lv[b + i][n] for i in range(D + 1)) % M for n in range(nn)]
    say('   P\'s law alone, the smallest depth over n ≥ 1 from levels %s: %s' % (
        sorted(G), [min(val(g, L, K) for g in G[b][1:]) if len(G[b]) > 1 else '-' for b in sorted(G)]))
    if not split(L): continue
    u, a = unit_root(L, K, FORM)
    pairs = [b for b in G if b + 1 in G]; bl = max(pairs) if pairs else None
    ests = []
    for n in (range(1, min(len(G[bl]), len(G[bl + 1]))) if bl is not None else []):
        v = val(G[bl][n], L, K)
        if v >= K - 1 or G[bl + 1][n] % L ** v: continue
        mk = L ** (K - v); ests.append(((G[bl + 1][n] // L ** v) * pow((G[bl][n] // L ** v) % mk, -1, mk) % mk, K - v))
    if ests:
        kk = min(k for _, k in ests); lam = ests[0][0] % L ** kk
        agree = next((e for e in range(kk) if any((x - lam) % L ** (e + 1) for x, _ in ests)), kk)
        hits = []
        for m in range(-60, 61):
            for sg in (1, -1):
                w = sg * pow(u, m, M) % M
                hits.append((next((e for e in range(agree) if (w - lam) % L ** (e + 1)), agree), sg, m))
        hits.sort(key=lambda h: (-h[0], abs(h[2])))
        say('   the extra constant from the data, G_%d/G_%d over %d values of n: digits %s (agreeing to %d);  the best ±u^m '
            '(digits, sign, m): %s' % (bl + 1, bl, len(ests), [(lam // L ** e) % L for e in range(agree)], agree, hits[:4]))
    else:
        say('   too few levels, or no remainder, to read the extra constant from the data')
    lamM = pow(u, E, M)
    say('   %s at ℓ: aₗ = ±%d up to twist; u^%d digits %s' % (curve, a, E, [(lamM // L ** e) % L for e in range(K)]))
    Q = [0] * (D + 2)
    for i, c in enumerate(P): Q[i + 1] = (Q[i + 1] + c) % M; Q[i] = (Q[i] - lamM * c) % M
    for b in sorted(lv):
        if b + D + 1 > top: continue
        nn = min(len(lv[b + i]) for i in range(D + 2))
        deps = [val(sum(Q[i] * lv[b + i][n] for i in range(D + 2)), L, K) for n in range(nn)]
        say('   P(x)·(x − u^%d) from level %d: depth for n = 0, 1, …  %s' % (E, b, deps))
