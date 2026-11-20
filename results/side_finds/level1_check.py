# A dig, 12 Nov 2026: j AND 1/E₄, THE FORMS' CONSTANTS AGAINST THE EXACT TOWERS (level1_tower.py), and the law the data
# themselves ask for. forms1_family_w<W>_ℓ.txt gives, from level-1 forms alone, the polynomial P(x) of the constants at
# weight W (j: W = 0; 1/E₄: W = −4): the Eisenstein family's 1 times the ordinary cusp forms' families.
# PREDICTED before the run (derived):
#   j, every ℓ: Σ_i [xⁱ]P · H_(b+i)(n) ≡ 0 on lane ℓ, deeper as b grows (j has its only pole at the cusp, as 1/Δ did); at
#     5, 7 and 11 (only the Eisenstein family is ordinary) H_(b+1) ≡ H_b: Lehner's congruences say the same.
#   1/E₄, ℓ ≡ 2 mod 3: the same, since ρ (j = 0) is then supersingular and 1/E₄ is an ℓ-adic modular form (E₄ divides a
#     power of the Hasse invariant on lane ℓ).
#   1/E₄, ℓ ≡ 1 mod 3: ρ is ordinary and 1/E₄ is not an ℓ-adic modular form; the pole may bring constants of its own, so P
#     alone may not be the whole law.
# For each prime: (1) the depth of P's law from each level; (2) where it does not deepen, the ratio of its remainders one
# level apart, G_(b+1)(n)/G_b(n) with G_b(n) = Σ_i [xⁱ]P·H_(b+i)(n): one more constant, if a single one is missing;
# (3) the least monic law on lane ℓ that the data alone allow (unknown coefficients, solved over F_ℓ from every level
# b ≥ 1 and every n), against P on lane ℓ. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/level1_check.py <tower json> … [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRCS = [a for a in sys.argv[1:] if a.endswith('.json')]
fh = None
if KEEP:
    name = 'level1_check_%s.txt' % '_'.join(os.path.basename(s)[13:-5] for s in SRCS)
    fh = open(os.path.join(HERE, name), 'x', encoding='utf-8')


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def solve_mod(rows, rhs, p):
    """Gaussian elimination over F_p: (consistent, solution or None, rank)."""
    m = len(rows[0]) if rows else 0
    A = [list(r) + [b] for r, b in zip(rows, rhs)]; piv = []; r = 0
    for c in range(m):
        k = next((i for i in range(r, len(A)) if A[i][c] % p), None)
        if k is None: continue
        A[r], A[k] = A[k], A[r]; inv = pow(A[r][c], -1, p)
        A[r] = [x * inv % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] % p:
                f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if any(all(x % p == 0 for x in row[:m]) and row[m] % p for row in A): return False, None, r
    sol = [0] * m
    for i, c in enumerate(piv): sol[c] = A[i][m]
    return True, sol, r


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


for src in SRCS:
    data = json.load(open(os.path.join(HERE, src), encoding='utf-8'))
    W = {'j': 0, 'invE4': -4, 'invE6': -6}[data['form']]
    say('towers: %s (n < %d), weight %d' % (src, data['N'], W))
    for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
        L = int(Ls); f = os.path.join(HERE, 'forms1_family_w%d_%d.txt' % (W, L))
        if not os.path.exists(f): say('   ℓ = %d: no family file' % L); continue
        mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(f, encoding='utf-8').read())
        KF = int(mt.group(1)); P = [int(v) for v in mt.group(2).split(',')]; D = len(P) - 1
        K = min(KF, rec['K']); M = L ** K; lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}
        top = max(lv)
        say('   ℓ = %d (%s), lane ℓ^%d; %d constant%s from the forms, P on lane ℓ %s' % (
            L, '≡ 1 mod 3' if L % 3 == 1 else '≡ 2 mod 3', K, D, '' if D == 1 else 's', [c % L for c in P]))
        G = {}
        for b in sorted(lv):                                         # (1) P's law, from every level b ≥ 0
            if b + D > top: continue
            nn = min(len(lv[b + i]) for i in range(D + 1))
            G[b] = [sum(P[i] * lv[b + i][n] for i in range(D + 1)) % M for n in range(nn)]
            say('      P\'s law from level %d: ℓ-adic depth for n = 0, 1, …  %s' % (b, [val(g, L, K) for g in G[b]]))
        rat = []                                                     # (2) one more constant, if a single one is missing
        for b in sorted(G):
            if b + 1 not in G or b < 1: continue
            for n in range(min(len(G[b]), len(G[b + 1]))):
                v = val(G[b][n], L, K)
                if v >= K - 1: continue
                u = G[b][n] // L ** v; w = G[b + 1][n]
                if w % L ** v: rat.append((b, n, None)); continue
                mk = L ** (K - v)
                rat.append((b, n, (w // L ** v) * pow(u % mk, -1, mk) % mk, K - v))
        if rat:
            good = [r for r in rat if r[2] is not None]
            if good:
                kmin = min(r[3] for r in good); ref = good[0][2] % L ** kmin
                same = all(r[2] % L ** kmin == ref for r in good)
                say('      the remainders one level apart, G_(b+1)/G_b, over %d pairs (b ≥ 1): on lane ℓ %s; '
                    'the same number to ℓ^%d for all: %s%s' % (
                        len(good), sorted(set(r[2] % L for r in good)), kmin, same,
                        (';  that number, digits lowest first: %s' % [(ref // L ** e) % L for e in range(kmin)]) if same else ''))
            bad = [r for r in rat if r[2] is None]
            if bad: say('      remainders that do not divide: %d' % len(bad))
        best = None                                                  # (3) the least law the data alone allow, on lane ℓ
        for d in range(1, 7):
            rows, rhs = [], []
            for b in sorted(lv):
                if b < 1 or b + d > top: continue
                nn = min(len(lv[b + i]) for i in range(d + 1))
                for n in range(nn):
                    rows.append([lv[b + i][n] % L for i in range(d)]); rhs.append((-lv[b + d][n]) % L)
            if len(rows) <= d: break
            ok, sol, rk = solve_mod(rows, rhs, L)
            if ok:
                best = (d, sol + [1], len(rows), rk); break
        if best:
            d, poly, ne, rk = best
            say('      the least law on lane ℓ from the data alone: degree %d, %s (lowest first), from %d equations, '
                'unique: %s;  the forms\' P on lane ℓ: %s' % (d, poly, ne, rk == d, [c % L for c in P]))
        else:
            say('      the data allow no law of degree ≤ 6 on lane ℓ (or too few levels)')
