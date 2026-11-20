# A dig, 12 Nov 2026: A POLE AT A POINT WITH NO COMPLEX MULTIPLICATION. F = 1/(E₄³ − j₀Δ) = 1/(Δ·(j − j₀)), weight w = −12,
# one simple pole, at the point where j = j₀ (level1_tower.py poleJ0, default j₀ = 2). The curve there:
#   C: y² = x³ + 3j₀(1728 − j₀)·x + 2j₀(1728 − j₀)²   (j-invariant j₀; for j₀ = 2: y² = x³ + 10356x + 11916304).
# The reasoning (derived, README): U_ℓ moves the pole by the lift of Frobenius φ (the quotient by the canonical subgroup).
# At a CM point with ℓ split φ fixes it, and the residue gains u^(w−2) per level. At j₀ = 2 the pole moves, φᵇ(C) → the
# canonical lift of C mod ℓ (a CM curve with the same point count), so the factor converges to u_C^(w−2) = u_C⁻¹⁴, u_C the
# unit root of x² − aₗx + ℓ, aₗ = ℓ + 1 − #C(𝔽ₗ) (the same for every twist).
# PREDICTED before the run:
#   - C ordinary at ℓ: the remainders of the weight −12 families' law P carry one more constant, a unit, = u_C⁻¹⁴; its
#     agreement may grow with the level rather than hold to every digit at once;
#   - C supersingular at ℓ (29 and 31 for j₀ = 2): P alone, deeper level by level.
# Then a BLIND reading, which does not use j₀: from the data alone, per prime, the traces a with u_a⁻¹⁴ equal to the
# constant read (or aₗ ≡ 0 where the remainders die), the residues j mod ℓ of curves with such a trace, and the integers j
# with |j| ≤ R that fit every prime. COUNTED.
# Fixed after the first run (12 Nov), before any reading of j: a constant is read only from unit ratios (the first run took
# the ratio 0 of dying remainders at 31 as a constant), and the blind reading uses the deepest read that at least 3 values of
# n confirm (the first run used one value of n at 5, whose last 2 of 11 digits were not yet settled).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/pole_j0_check.py <poleJ0 tower json> [--R=…] [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRC = next(a for a in sys.argv[1:] if a.endswith('.json'))
R = next((int(a[4:]) for a in sys.argv[1:] if a.startswith('--R=')), 10 ** 6)
fh = open(os.path.join(HERE, 'pole_j0_check_%s.txt' % os.path.basename(SRC)[13:-5]), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def trace(A, B, L):
    """aₗ of y² = x³ + Ax + B over 𝔽ₗ, by counting points."""
    sq = [0] * L
    for y in range(L): sq[y * y % L] += 1
    return L + 1 - (1 + sum(sq[(x ** 3 + A * x + B) % L] for x in range(L)))


def traces_of_j(j, L):
    """the traces of every curve over 𝔽ₗ with j-invariant j (all twists)."""
    j %= L
    if j == 0: return {trace(0, B, L) for B in range(1, L)}
    if j == 1728 % L: return {trace(A, 0, L) for A in range(1, L)}
    A, B = 3 * j * (1728 - j) % L, 2 * j * (1728 - j) ** 2 % L
    a = trace(A, B, L)
    return {a, -a}


def unit_root(a, L, K):
    M = L ** K; u = a % L
    for _ in range(4 * K + 8): u = (u - (u * u - a * u + L) * pow((2 * u - a) % M, -1, M)) % M
    return u


def digits(x, L, k): return [(x // L ** e) % L for e in range(k)]


data = json.load(open(os.path.join(HERE, SRC), encoding='utf-8')); assert data['form'] == 'poleJ0'
J0 = data['j0']; W = -12; EX = W - 2
A0, B0 = 3 * J0 * (1728 - J0), 2 * J0 * (1728 - J0) ** 2
say('1/(E₄³ − %dΔ) towers: %s (n < %d); the curve at the pole y² = x³ + %d·x + %d (j = %d); the predicted constant u^(w−2) = u^%d'
    % (J0, SRC, data['N'], A0, B0, J0, EX))
blind = {}
for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
    L = int(Ls); a = trace(A0, B0, L); ordinary = a % L != 0
    mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(os.path.join(HERE, 'forms1_family_w%d_%d.txt' % (W, L)), encoding='utf-8').read())
    K = min(int(mt.group(1)), rec['K']); M = L ** K
    P = [int(v) % M for v in mt.group(2).split(',')]; D = len(P) - 1
    lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}; top = max(lv)
    G = {b: [sum(P[i] * lv[b + i][n] for i in range(D + 1)) % M for n in range(min(len(lv[b + i]) for i in range(D + 1)))]
         for b in lv if b + D <= top}
    say(); say('ℓ = %d, lane ℓ^%d: aₗ(C) = %d (%s); P (the weight −12 families) on lane ℓ %s' % (
        L, K, a, 'ordinary' if ordinary else 'supersingular', [c % L for c in P]))
    mins = {b: min((val(g, L, K) for g in G[b][1:]), default=None) for b in sorted(G)}
    say('   P\'s law alone, the smallest depth over n ≥ 1, by level: %s' % mins)
    # the constant read from the data, level pair by level pair: G_(b+1)(n)/G_b(n) for n ≥ 1 where G_b(n) is a unit times ℓ^v
    reads = {}
    for b in sorted(G):
        if b < 1 or b + 1 not in G: continue
        ests = []
        for n in range(1, min(len(G[b]), len(G[b + 1]))):
            v0 = val(G[b][n], L, K)
            if v0 >= K - 1 or val(G[b + 1][n], L, K) != v0: continue        # a unit ratio only
            mk = L ** (K - v0); ests.append((G[b + 1][n] // L ** v0) * pow((G[b][n] // L ** v0) % mk, -1, mk) % mk)
            ests[-1] = (ests[-1], K - v0)
        if not ests: continue
        kk = min(k for _, k in ests); lam = ests[0][0] % L ** kk
        agree = next((e for e in range(kk) if any((x - lam) % L ** (e + 1) for x, _ in ests)), kk)
        if agree == 0: continue
        reads[b] = (lam % L ** agree, agree, len(ests))
    if ordinary:
        u = unit_root(a, L, K); c14 = pow(u, EX, M)
        say('   u_C^%d digits %s' % (EX, digits(c14, L, K)))
        for b, (lam, ag, ne) in sorted(reads.items()):
            hit = next((e for e in range(ag) if (lam - c14) % L ** (e + 1)), ag)
            say('   levels %d → %d: the data\'s constant (from %d values of n, agreeing to %d digits) %s; equal to u_C^%d to %d of %d digits'
                % (b, b + 1, ne, ag, digits(lam, L, ag), EX, hit, ag))
        if reads:
            b, (lam, ag, ne) = max(reads.items(), key=lambda kv: (kv[1][1], kv[0]))
            best = []
            for mm in range(-60, 61):
                for sg in (1, -1):
                    w_ = sg * pow(u, mm, M) % M
                    best.append((next((e for e in range(ag) if (w_ - lam) % L ** (e + 1)), ag), sg, mm))
            best.sort(key=lambda h: (-h[0], abs(h[2])))
            say('   the best ±u_C^m for the read at levels %d → %d (digits, sign, m): %s' % (b, b + 1, best[:4]))
        Q = [0] * (D + 2)
        for i, c in enumerate(P): Q[i + 1] = (Q[i + 1] + c) % M; Q[i] = (Q[i] - c14 * c) % M
        rows = []
        for b in sorted(lv):
            if b < 1 or b + D + 1 > top: continue
            nn = min(len(lv[b + i]) for i in range(D + 2))
            rows.append('%d: %d' % (b, min(val(sum(Q[i] * lv[b + i][n] for i in range(D + 2)), L, K) for n in range(1, nn))))
        say('   the law P(x)·(x − u_C^%d), the smallest depth over n ≥ 1 from level: %s' % (EX, ';  '.join(rows) or 'no level span'))
    # the blind reading: what the data alone say about the curve at the pole
    bs = [b for b in sorted(mins) if b >= 1 and mins[b] is not None]
    dying = len(bs) >= 2 and all(mins[bs[i + 1]] > mins[bs[i]] for i in range(len(bs) - 1))
    if dying:
        tr = {0}; why = 'the remainders die level by level: aₗ ≡ 0'
    elif [r for r in reads.values() if r[2] >= 3]:                    # the deepest read that at least 3 values of n confirm
        b, (lam, ag, ne) = max(((b, r) for b, r in reads.items() if r[2] >= 3), key=lambda kv: (kv[1][1], kv[0]))
        tr = set()
        for t in range(-int(2 * L ** 0.5), int(2 * L ** 0.5) + 1):
            if t % L == 0: continue
            if (pow(unit_root(t, L, ag), EX, L ** ag) - lam) % L ** ag == 0: tr.add(t)
        why = 'the constant at levels %d → %d, %d digits' % (b, b + 1, ag)
    else:
        say('   blind reading: nothing to read'); continue
    js = sorted(j for j in range(L) if traces_of_j(j, L) & tr)
    blind[L] = js
    say('   blind reading (%s): traces %s; j mod ℓ in %s (%d of %d); the true j₀ = %d mod ℓ = %d is %s' % (
        why, sorted(tr), js, len(js), L, J0, J0 % L, 'among them' if J0 % L in js else 'NOT among them'))
say()
if blind:
    mod = 1; sols = [0]
    for L, js in sorted(blind.items()):
        sols = [s + mod * t for s in sols for t in range(L) if (s + mod * t) % L in js]; mod *= L
        if len(sols) > 2 * 10 ** 6: say('too many residues to list; stopped at ℓ = %d' % L); break
    fits = sorted({x for s in sols for x in range(s - (R // mod + 1) * mod, R + 1, mod) if -R <= x <= R} if mod <= 4 * R else
                  {x for s in sols for x in (s, s - mod) if -R <= x <= R}, key=abs)
    exp = 1.0
    for L, js in blind.items(): exp *= len(js) / L
    say('blind reading over %d primes: %d residues mod %d; the integers j with |j| ≤ %d that fit every prime: %d (by chance about %.1f)'
        % (len(blind), len(sols), mod, R, len(fits), exp * (2 * R + 1)))
    say('   the smallest: %s' % fits[:20])
    say('   the true j₀ = %d is %s' % (J0, 'among them' if J0 in fits else 'NOT among them'))
