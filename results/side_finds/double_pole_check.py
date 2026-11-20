# A dig, 12 Nov 2026: DOUBLE POLES. Two forms with a double pole at ρ (the curve y² = x³ + 1):
#   E₄/j = Δ/E₄² (weight 4). Zhang (arXiv 2510.23200, Thm 1.2) PROVED a_(npˡ) ≡ (−3|p)·p·a_(npˡ⁻¹) (mod p^(3l)) for every p ≥ 5:
#     a single constant (−3|p)·p. This is a check of our tools on a proven case.
#   1/E₄² (weight −8): new.
# The residue reasoning of the simple pole (derived; README "1/E₄"), extended to a pole of order m:
#   - T_ℓ multiplies the coefficient of (τ − ρ)^(−j) by π^(j−1)·π̄^(w−1−j) + π̄^(j−1)·π^(w−1−j).
#   - The terms allowed at ρ are those with −j ≡ w (mod 3), by the stabilizer of ρ: for both forms only j = 2.
#   - At ℓ ≡ 1 mod 3 the integral one of the pair is π·π̄^(w−3) = ℓ·u^(w−4), u the unit root of the curve's Frobenius.
#   - For E₄/j that is ℓ·u⁰ = ℓ = (−3|ℓ)·ℓ, Zhang's constant (a check of the reasoning).
# PREDICTED before the run, for 1/E₄²:
#   - ℓ ≡ 1 mod 3: the law is P(x)·(x − ℓ·u⁻¹²), with P the weight −8 families' polynomial (forms1_family_w-8_ℓ.txt);
#   - ℓ ≡ 2 mod 3: P alone, or P with one constant of slope ≥ 1 that the data must name (Zhang's −p at w = 4 suggests one).
# Per prime: the extra constant read from the data (ratio of remainders of P's law one level apart, with its ℓ-adic
# valuation), and the depth of the predicted laws. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/double_pole_check.py <E4overj json> <invE4sq json> [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRCS = [a for a in sys.argv[1:] if a.endswith('.json')]
fh = open(os.path.join(HERE, 'double_pole_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def unit_root(L, K):
    M = L ** K
    a = next(a for a in range(1, 2 * L) if 4 * L - a * a > 0 and (4 * L - a * a) % 3 == 0
             and round(((4 * L - a * a) // 3) ** 0.5) ** 2 == (4 * L - a * a) // 3)
    u = a % L
    for _ in range(4 * K + 8): u = (u - (u * u - a * u + L) * pow((2 * u - a) % M, -1, M)) % M
    return u


def ratios(G, L, K):
    """G_(b+1)(n)/G_b(n) as (valuation, unit part mod ℓ^prec, prec) over b, n where defined."""
    out = []
    for b in sorted(G):
        if b + 1 not in G: continue
        for n in range(min(len(G[b]), len(G[b + 1]))):
            v0, v1 = val(G[b][n], L, K), val(G[b + 1][n], L, K)
            if v0 >= K - 1 or v1 >= K: continue
            prec = K - max(v0, v1)
            if prec < 1: continue
            u0 = (G[b][n] // L ** v0) % L ** prec; u1 = (G[b + 1][n] // L ** v1) % L ** prec
            out.append((b, n, v1 - v0, u1 * pow(u0, -1, L ** prec) % L ** prec, prec))
    return out


for src in SRCS:
    data = json.load(open(os.path.join(HERE, src), encoding='utf-8')); FORM = data['form']
    W = {'E4overj': 4, 'invE4sq': -8}[FORM]
    say(); say('%s towers: %s (n < %d), weight %d' % (FORM, src, data['N'], W))
    for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
        L = int(Ls); K = rec['K']; M = L ** K; lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}
        if FORM == 'E4overj':
            P = [M - 1, 1]; Pname = 'x − 1 (the Eisenstein family; weight 4 has no cusp form)'
        else:
            f = os.path.join(HERE, 'forms1_family_w-8_%d.txt' % L)
            mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(f, encoding='utf-8').read())
            KF = int(mt.group(1)); K = min(K, KF); M = L ** K
            P = [int(v) % M for v in mt.group(2).split(',')]; lv = {b: [x % M for x in v] for b, v in lv.items()}
            Pname = 'the weight −8 families, on lane ℓ %s' % [c % L for c in P]
        D = len(P) - 1; top = max(lv)
        G = {}
        for b in sorted(lv):
            if b + D <= top:
                nn = min(len(lv[b + i]) for i in range(D + 1))
                G[b] = [sum(P[i] * lv[b + i][n] for i in range(D + 1)) % M for n in range(nn)]
        split = L % 3 == 1
        say('ℓ = %d (%s), lane ℓ^%d; P = %s' % (L, 'ρ ordinary' if split else 'ρ supersingular', K, Pname))
        say('   P\'s law alone, smallest depth over n ≥ 1 from levels %s: %s' % (sorted(G), [min((val(g, L, K) for g in G[b][1:]), default='-') for b in sorted(G)]))
        rs = [r for r in ratios(G, L, K) if r[0] >= 1 and r[1] >= 1]
        if rs:
            vs = sorted(set(r[2] for r in rs))
            best = max(rs, key=lambda r: (r[0], r[4]))
            c_unit, prec = best[3], best[4]
            say('   the extra constant from the data: valuations of G_(b+1)/G_b seen %s; at the last level, unit part digits %s' % (
                vs, [(c_unit // L ** e) % L for e in range(min(prec, 6))]))
            if split:
                u = unit_root(L, K); cand = pow(u, W - 4, M)
                say('   predicted ℓ·u^(w−4) = ℓ·u^%d: valuation 1, unit part digits %s' % (W - 4, [(cand // L ** e) % L for e in range(min(prec, 6))]))
            else:
                say('   compare (−3|ℓ)·ℓ: valuation 1, unit part digits %s' % [((-1) % L ** prec // L ** e) % L for e in range(min(prec, 6))])
        # the law P(x)·(x − c) with c = the predicted constant (split) or (−3|ℓ)·ℓ (inert)
        c = (L * pow(unit_root(L, K), W - 4, M)) % M if split else (leg(-3, L) * L) % M
        Q = [0] * (D + 2)
        for i, ci in enumerate(P): Q[i + 1] = (Q[i + 1] + ci) % M; Q[i] = (Q[i] - c * ci) % M
        rows = []
        for b in sorted(lv):
            if b < 1 or b + D + 1 > top: continue
            nn = min(len(lv[b + i]) for i in range(D + 2))
            ds = [val(sum(Q[i] * lv[b + i][n] for i in range(D + 2)), L, K) for n in range(1, nn)]
            if ds: rows.append('%d: %d' % (b, min(ds)))
        say('   the law P(x)·(x − %s) from levels b ≥ 1, smallest depth over n ≥ 1: %s' % (
            'ℓ·u^%d' % (W - 4) if split else '(−3|ℓ)·ℓ', ';  '.join(rows) or 'no level span'))
