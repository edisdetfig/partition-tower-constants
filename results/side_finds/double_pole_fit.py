# A dig, 12 Nov 2026: 1/E₄², THE TWO SLOPE-1 CONSTANTS SEPARATED. double_pole_check.py showed that the remainders of the
# ordinary families' law, G_b(n) = Σ_i [xⁱ]P·H_(b+i)(n), carry two constants of slope 1:
#   - the pole's: ℓ·u^(w−4) = ℓ·u⁻¹² (derived; u the unit root of Frobenius of y² = x³ + 1 at ℓ ≡ 1 mod 3);
#   - one of a slope-1 holomorphic family.
# Here the two-term law G_(b+2) = s·G_(b+1) − p·G_b is fitted over every b ≥ 1 and n ≥ 1, modulo ℓ^K. Then
# x = ℓ·y turns x² − s·x + p into y² − (s/ℓ)·y + p/ℓ², whose two roots y₁, y₂ (Hensel, when distinct on lane ℓ) are the
# unit parts of the two constants.
# PREDICTED before this run (from the derivation): one root y equals u⁻¹² to every digit the data hold. COUNTED.
# Also (added before the deeper run): a d-term law for d = 2, 3 solved over the ℓ-adic integers (padic_kernel). Its slope-1
# factor (slope_factor, roots/ℓ) is evaluated at u⁻¹², in case more than one slope-1 holomorphic family is present.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/double_pole_fit.py <invE4sq tower json> … [--keep]
import os, sys, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forms6_flint as F
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRCS = [a for a in sys.argv[1:] if a.endswith('.json')]
fh = open(os.path.join(HERE, 'double_pole_fit_%s.txt' % '_'.join(os.path.basename(s)[13:-5] for s in SRCS)), 'x', encoding='utf-8') if KEEP else None


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
    M = L ** K
    a = next(a for a in range(1, 2 * L) if 4 * L - a * a > 0 and (4 * L - a * a) % 3 == 0
             and round(((4 * L - a * a) // 3) ** 0.5) ** 2 == (4 * L - a * a) // 3)
    u = a % L
    for _ in range(4 * K + 8): u = (u - (u * u - a * u + L) * pow((2 * u - a) % M, -1, M)) % M
    return u


def digits(x, L, k): return [(x // L ** e) % L for e in range(k)]


for src in SRCS:
    data = json.load(open(os.path.join(HERE, src), encoding='utf-8')); assert data['form'] == 'invE4sq'
    say('1/E₄² towers: %s (n < %d)' % (src, data['N']))
    for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
        L = int(Ls)
        if L % 3 != 1: continue
        K = rec['K']
        mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(os.path.join(HERE, 'forms1_family_w-8_%d.txt' % L), encoding='utf-8').read())
        K = min(K, int(mt.group(1))); M = L ** K
        P = [int(v) % M for v in mt.group(2).split(',')]; D = len(P) - 1
        lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}; top = max(lv)
        u12 = pow(unit_root(L, K), -12, M)
        G = {b: [sum(P[i] * lv[b + i][n] for i in range(D + 1)) % M for n in range(min(len(lv[b + i]) for i in range(D + 1)))]
             for b in lv if b + D <= top}
        eqs = [(G[b + 1][n], (-G[b][n]) % M, G[b + 2][n], b, n) for b in sorted(G) if b >= 1 and b + 2 in G
               for n in range(1, min(len(G[b]), len(G[b + 1]), len(G[b + 2])))]
        best = None
        for i in range(len(eqs)):
            for j in range(i + 1, len(eqs)):
                det = (eqs[i][0] * eqs[j][1] - eqs[j][0] * eqs[i][1]) % M; vd = val(det, L, K)
                if best is None or vd < best[0]: best = (vd, i, j)
        if best is None: say('ℓ = %d: too few equations' % L); continue
        vd, i, j = best; (a1, b1, c1), (a2, b2, c2) = eqs[i][:3], eqs[j][:3]
        prec = K - vd; Mp = L ** prec; det = (a1 * b2 - a2 * b1) % M; inv = pow((det // L ** vd) % Mp, -1, Mp)
        s = (((c1 * b2 - c2 * b1) % M) // L ** vd) * inv % Mp; p = (((a1 * c2 - a2 * c1) % M) // L ** vd) * inv % Mp
        cons = [all((a * s + bb * p - c) % L ** e == 0 for (a, bb, c, _, _) in eqs) for e in range(1, prec + 1)]
        okdepth = max([e + 1 for e in range(prec) if cons[e]], default=0)
        say('ℓ = %d, lane ℓ^%d: the two-term law x² − s·x + p fitted from %d equations; s and p known to ℓ^%d; every equation '
            'holds to ℓ^%d; v(s) = %d, v(p) = %d' % (L, K, len(eqs), prec, okdepth, val(s, L, prec), val(p, L, prec)))
        if val(s, L, prec) < 1 or val(p, L, prec) < 2: say('   not two slope-1 roots; stopped'); continue
        q = prec - 2; Mq = L ** q                                    # y² − σ·y + π0 with σ = s/ℓ, π0 = p/ℓ² (mod ℓ^q)
        sig, pi0 = (s // L) % Mq, (p // L ** 2) % Mq
        roots = []
        for y0 in range(L):
            if (y0 * y0 - sig * y0 + pi0) % L == 0 and (2 * y0 - sig) % L: roots.append(y0)
        lifted = []
        for y in roots:
            for _ in range(2 * q + 4): y = (y - (y * y - sig * y + pi0) * pow((2 * y - sig) % Mq, -1, Mq)) % Mq
            lifted.append(y)
        u12q = u12 % Mq
        say('   the two unit parts y (digits lowest first, to ℓ^%d): %s' % (q, [digits(y, L, q) for y in lifted] if lifted else 'not distinct on lane ℓ'))
        say('   u⁻¹² digits: %s' % digits(u12q, L, q))
        for y in lifted:
            agree = next((e for e in range(q) if (y - u12q) % L ** (e + 1)), q)
            say('   root %s: agrees with u⁻¹² to %d of %d digits' % (digits(y, L, min(q, 4)), agree, q))

        for d in (2, 3):                                             # the d-term law over ℤ_ℓ, its slope-1 factor at u⁻¹²
            rows = []
            for b in sorted(G):
                if b < 1 or b + d not in G: continue
                for n in range(1, min(len(G[b + i]) for i in range(d + 1))):
                    rows.append([G[b + i][n] % M for i in range(d)] + [G[b + d][n] % M])
            if len(rows) < d + 2: say('   %d-term law: too few equations (%d)' % (d, len(rows))); continue
            found = None
            for Kp in range(K, 1, -1):                              # the largest lane on which the data allow one law
                ker, Np = F.padic_kernel([[x % L ** Kp for x in r] for r in rows], len(rows), d + 1, L, Kp)
                if len(ker) == 1 and ker[0][d] % L and Np >= 2: found = (ker, Np, Kp); break
            if found is None: say('   %d-term law: no unique law on any lane (%d equations)' % (d, len(rows))); continue
            ker, Np, Kp = found
            k0 = ker[0]; Mn = L ** Np; inv = pow(k0[d] % Mn, -1, Mn)
            chi = [k0[i] * inv % Mn for i in range(d)] + [1]       # x^d + Σ c_i x^i, lowest first
            vs = [val(c, L, Np) for c in chi]
            Gf, pr = F.slope_factor(chi, L, 1, Np)
            if Gf is None: say('   %d-term law (%d equations, to ℓ^%d): valuations of its coefficients %s; no slope-1 factor' % (d, len(rows), Np, vs)); continue
            ev = sum(c * pow(u12 % L ** pr, i, L ** pr) for i, c in enumerate(Gf)) % L ** pr
            say('   %d-term law (%d equations, unique on lane ℓ^%d, known to ℓ^%d): coefficient valuations %s; slope-1 factor of degree '
                '%d (to ℓ^%d); at y = u⁻¹² it vanishes to ℓ^%d of %d' % (d, len(rows), Kp, Np, vs, len(Gf) - 1, pr, val(ev, L, pr), pr))
