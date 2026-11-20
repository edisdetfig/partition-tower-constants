# A dig, 12 Nov 2026: CURVES THAT REPEAT THEMSELVES (README "Curves that repeat themselves"; predictions written there first).
# U_ℓ moves a pole by φ, the lift of Frobenius. φ's fixed points are the canonical lifts (CM curves with ℓ split); a CM curve
# whose prime above ℓ has order r in the class group comes back after r steps; a curve without CM drifts toward the
# canonical lift of its reduction.
# Part A: a pole ON the repeating curve, 1/(E₄³ − j_lift·Δ), against the drifting pole of 1/(E₄³ − 2Δ), at the primes where
#   j = 2's lift has class number 1 (5: −32768, 7 and 19: 54000, 11: −3375). Laws P(x)·(x − u⁻¹⁴)·S(x), S the slope-1
#   families at weight −12 (slope1_families_ℓ_w-12_N….txt), P the ordinary ones (forms1_family_w-12_ℓ.txt).
# Part B: 1/(Δ²·H₋₁₅(j)), two CM poles of class number 2. Inert ℓ: P alone. ℓ = x² + xy + 4y² (19, 31): each curve fixed,
#   P(x)·(x − u⁻²⁶). Split but not of that form (17, 23): the curves swap, P(x)·(x² − U⁻²⁶), U the unit root of the curve's
#   Frobenius over 𝔽_ℓ² (counted on the curve, checked against x² + 15y² = 4ℓ²). S (weight −24 slope-1 families) where computed.
# Per law and starting level b ≥ 1: the smallest depth over n ≥ 1. COUNTED.
# Part C: 1/(Δ³·H₋₂₃(j)), three CM poles of class number 3. Where ℓ's prime has order f in the class group (H mod ℓ splits
#   into factors of degree f), the curves cycle with period f: law P(x)·(x^f − U^(w−2))·S(x), U the unit root of the curves'
#   Frobenius over 𝔽_(ℓ^f) (counted on the curve over 𝔽_ℓ[t]/(g), checked against x² + 23y² = 4ℓ^f). For each ℓ the tower
#   with the most levels is used (n < 1.6·10⁷ on Windows, n < 10⁸ in WSL).
#   Part C also takes --disc=−31 (towers level1_tower_poleHD31_N…), the class polynomial read from the tower file.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/repeat_check.py A|B|C [--disc=…] [--keep] [--tag=…]
import os, sys, re, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
PART = next(a for a in sys.argv[1:] if a in ('A', 'B', 'C'))
TAG = next((a[6:] for a in sys.argv[1:] if a.startswith('--tag=')), '')
fh = open(os.path.join(HERE, 'repeat_check_%s%s.txt' % (PART, TAG)), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def digits(x, L, k): return [(x // L ** e) % L for e in range(k)]


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def trace(A, B, L):
    sq = [0] * L
    for y in range(L): sq[y * y % L] += 1
    return L + 1 - (1 + sum(sq[(x ** 3 + A * x + B) % L] for x in range(L)))


def curve_of_j(j, L): return 3 * j * (1728 - j) % L, 2 * j * (1728 - j) ** 2 % L


def unit_root(t, q, L, K):
    """the root of x² − t·x + q (q a power of ℓ) that is a unit, modulo ℓ^K."""
    M = L ** K; u = t % L
    for _ in range(4 * K + 8): u = (u - (u * u - t * u + q) * pow((2 * u - t) % M, -1, M)) % M
    return u


def pmul(A, B, M):
    out = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        for j, b in enumerate(B): out[i + j] = (out[i + j] + a * b) % M
    return out


def family_P(w, L):
    mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(os.path.join(HERE, 'forms1_family_w%d_%d.txt' % (w, L)), encoding='utf-8').read())
    return [int(v) for v in mt.group(2).split(',')], int(mt.group(1))


def family_S(w, L):
    """∏(x − ℓy) over the slope-1 families at weight w (lowest first) and its precision, or ([1], big) if none computed."""
    srcs = sorted(glob.glob(os.path.join(HERE, 'slope1_families_%d_w%d_N*.txt' % (L, w))))
    if not srcs: return None, 0
    ms = re.search(r'the slope-1 factor ∏\(y − λ/ℓ\) modulo ℓ\^(\d+) \(truncated at t = \d+\), lowest first: \[([^\]]*)\]',
                   open(srcs[-1], encoding='utf-8').read())
    G1 = [int(v) for v in ms.group(2).split(',')]; h1 = len(G1) - 1
    return [G1[s] * L ** (h1 - s) for s in range(h1 + 1)], int(ms.group(1)) + 1


def depths(lv, Q, L, K):
    D = len(Q) - 1; top = max(lv); rows = []; checks = 0
    for b in sorted(lv):
        if b < 1 or b + D > top: continue
        nn = min(len(lv[b + i]) for i in range(D + 1))
        ds = [val(sum(Q[i] * lv[b + i][n] for i in range(D + 1)), L, K) for n in range(1, nn)]
        if ds: rows.append((b, min(ds))); checks += len(ds)
    return rows, checks


def show(rows, checks): return ('  '.join('%d: %d' % r for r in rows) or 'no level span') + '  (%d checks)' % checks


def load(pattern, L):
    for f in sorted(glob.glob(os.path.join(HERE, pattern))):
        d = json.load(open(f, encoding='utf-8'))
        if str(L) in d['primes']: return d, os.path.basename(f)
    return None, None


if PART == 'A':
    say('Part A: a pole on the repeating curve (the canonical lift) against the drifting pole at j = 2; weight −12, u⁻¹⁴')
    for L, jl in ((5, -32768), (7, 54000), (11, -3375), (19, 54000)):
        d2, f2 = load('level1_tower_poleJ02_N*.json', L); dl, fl = load('level1_tower_poleJ0%d_N*.json' % jl, L)
        if d2 is None or dl is None: say('ℓ = %d: tower missing' % L); continue
        P, KP = family_P(-12, L); S, KS = family_S(-12, L)
        a2, al = trace(*curve_of_j(2, L), L), trace(*curve_of_j(jl, L), L)
        s = val(2 - jl, L, 60)
        K = min(d2['primes'][str(L)]['K'], dl['primes'][str(L)]['K'], KP, KS); M = L ** K
        P = [c % M for c in P]; S = [c % M for c in S]
        c14 = pow(unit_root(a2, L, L, K), -14, M); lin = [(-c14) % M, 1]
        say(); say('ℓ = %d, lane ℓ^%d: the lift j = %d (aₗ = %d) and j = 2 (aₗ = %d) reduce to the same curve; s = v_ℓ(2 − j_lift) = %d; '
                   'slope-1 families at weight −12: %d' % (L, K, jl, al, a2, s, len(S) - 1))
        for name, Q in (('P(x)·(x − u⁻¹⁴)', pmul(P, lin, M)), ('P(x)·(x − u⁻¹⁴)·S(x)', pmul(pmul(P, lin, M), S, M))):
            for tag, d in (('pole ON the lift (repeats)', dl), ('pole at j = 2 (drifts)', d2)):
                lv = {int(b): [x % M for x in v] for b, v in d['primes'][str(L)]['levels'].items()}
                say('   %-22s %-27s depth by starting level: %s' % (name, tag, show(*depths(lv, Q, L, K))))
elif PART == 'B':
    bestB = {}                                                       # for each ℓ the tower with the most levels
    for f in sorted(glob.glob(os.path.join(HERE, 'level1_tower_poleHD15_N*.json'))):
        d = json.load(open(f, encoding='utf-8'))
        for Ls, rec in d['primes'].items():
            if Ls not in bestB or d['N'] > bestB[Ls][0]: bestB[Ls] = (d['N'], rec, os.path.basename(f), d['jpoly'])
    H = next(iter(bestB.values()))[3]; W = -12 * (len(H) - 1)
    say('Part B: the towers level1_tower_poleHD15_N…, H = %s (lowest first), weight %d, so u^(w−2) = u^%d' % (H, W, W - 2))
    for Ls, (Nn, rec, fn, _) in sorted(bestB.items(), key=lambda kv: int(kv[0])):
        L = int(Ls); P, KP = family_P(W, L); S, KS = family_S(W, L)
        say(); say('[ℓ = %d from %s, n < %d]' % (L, fn, Nn))
        K = min(rec['K'], KP, KS if S else 99); M = L ** K; P = [c % M for c in P]
        S = [c % M for c in S] if S else [1]
        lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}
        sym = leg(-15, L); princ = any(x * x + x * y + 4 * y * y == L for x in range(-L, L + 1) for y in range(-L, L + 1))
        if sym <= 0:
            say('ℓ = %d, lane ℓ^%d: %s; P on lane ℓ %s' % (L, K, 'ramified' if sym == 0 else 'inert (both curves supersingular)', [c % L for c in P]))
            say('   P alone (PREDICTED): %s' % show(*depths(lv, P, L, K)))
            continue
        roots = [r for r in range(L) if (r * r * H[2] + r * H[1] + H[0]) % L == 0]
        if princ:
            tr = [trace(*curve_of_j(r, L), L) for r in roots]
            t_norm = next(x for x in range(1, 2 * L) for y in range(1, L) if x * x + 15 * y * y == 4 * L)
            u = unit_root(tr[0], L, L, K); c = pow(u, W - 2, M)
            say('ℓ = %d, lane ℓ^%d: split, ℓ = x² + xy + 4y²: each curve fixed; roots of H on lane ℓ %s, their traces %s (norm form: ±%d); '
                'P on lane ℓ %s; slope-1 families %d' % (L, K, roots, tr, t_norm, [c_ % L for c_ in P], len(S) - 1))
            say('   u^%d digits %s' % (W - 2, digits(c, L, K)))
            for name, Q in (('P alone', P), ('P·(x − u^%d)' % (W - 2), pmul(P, [(-c) % M, 1], M)),
                            ('P·S', pmul(P, S, M)), ('P·(x − u^%d)·S (PREDICTED)' % (W - 2), pmul(pmul(P, [(-c) % M, 1], M), S, M))):
                say('   %-32s %s' % (name, show(*depths(lv, Q, L, K))))
        else:
            # the roots of H lie in 𝔽_ℓ² (√5 ∉ 𝔽_ℓ): count the curve over 𝔽_ℓ² = 𝔽_ℓ[t]/(t² − nr)
            nr = next(z for z in range(2, L) if leg(z, L) == -1); q = L * L
            mul2 = lambda a, b: ((a[0] * b[0] + nr * a[1] * b[1]) % L, (a[0] * b[1] + a[1] * b[0]) % L)
            add2 = lambda a, b: ((a[0] + b[0]) % L, (a[1] + b[1]) % L)
            def pw(a, e):
                r = (1, 0)
                while e:
                    if e & 1: r = mul2(r, a)
                    a = mul2(a, a); e >>= 1
                return r
            disc = (H[1] * H[1] - 4 * H[0] * H[2]) % L                 # √disc = c·t with c² = disc/nr
            cc = next(z for z in range(L) if z * z % L == disc * pow(nr, -1, L) % L)
            j1 = ((-H[1]) * pow(2, -1, L) % L, cc * pow(2, -1, L) % L)
            jj = j1; m1728 = ((1728 - jj[0]) % L, (-jj[1]) % L)
            A = mul2((3, 0), mul2(jj, m1728)); B = mul2((2, 0), mul2(jj, mul2(m1728, m1728)))
            assert add2(mul2(jj, jj), add2(mul2((H[1] % L, 0), jj), (H[0] % L, 0))) == (0, 0)
            s_ = 0
            for x0 in range(L):
                for x1 in range(L):
                    x = (x0, x1); z = add2(add2(mul2(x, mul2(x, x)), mul2(A, x)), B)
                    if z == (0, 0): continue
                    s_ += 1 if pw(z, (q - 1) // 2) == (1, 0) else -1
            aq = -s_
            t_norm = next(x for x in range(1, 2 * q) for y in range(1, 2 * L) if x * x + 15 * y * y == 4 * q and (x % L or y % L))
            U = unit_root(aq, q, L, K); c2 = pow(U, W - 2, M); c1 = pow(U, (W - 2) // 2, M)
            say('ℓ = %d, lane ℓ^%d: split, not of the form x² + xy + 4y²: the two curves SWAP; the curve over 𝔽_ℓ² has trace %d '
                '(norm form: ±%d); P on lane ℓ %s; slope-1 families %d' % (L, K, aq, t_norm, [c_ % L for c_ in P], len(S) - 1))
            say('   U^%d digits %s;  U^%d digits %s' % (W - 2, digits(c2, L, K), (W - 2) // 2, digits(c1, L, K)))
            sq = [(-c2) % M, 0, 1]
            for name, Q in (('P alone', P), ('P·S', pmul(P, S, M)),
                            ('P·(x − U^%d)·S  (one of the pair)' % ((W - 2) // 2), pmul(pmul(P, [(-c1) % M, 1], M), S, M)),
                            ('P·(x + U^%d)·S  (the other)' % ((W - 2) // 2), pmul(pmul(P, [c1 % M, 1], M), S, M)),
                            ('P·(x² − U^%d)   (the pair, no S)' % (W - 2), pmul(P, sq, M)),
                            ('P·(x² − U^%d)·S (PREDICTED)' % (W - 2), pmul(pmul(P, sq, M), S, M))):
                say('   %-34s %s' % (name, show(*depths(lv, Q, L, K))))
            G = {b: [sum(P[i] * lv[b + i][n] for i in range(len(P))) % M for n in range(min(len(lv[b + i]) for i in range(len(P))))]
                 for b in lv if b + len(P) - 1 <= max(lv)}
            for step in (1, 2):
                reads = []
                for b in sorted(G):
                    if b < 1 or b + step not in G: continue
                    for n in range(1, min(len(G[b]), len(G[b + step]))):
                        v0 = val(G[b][n], L, K)
                        if v0 >= K - 1 or val(G[b + step][n], L, K) != v0: continue
                        mk = L ** (K - v0); reads.append((b, (G[b + step][n] // L ** v0) * pow((G[b][n] // L ** v0) % mk, -1, mk) % mk))
                say('   the remainders\' ratio %d level%s apart, on lane ℓ, by level: %s' % (step, '' if step == 1 else 's',
                    {b: sorted({r % L for bb, r in reads if bb == b}) for b in sorted({bb for bb, _ in reads})}))
            say('   (on lane ℓ: U^%d ≡ %d, ±U^%d ≡ ±%d)' % (W - 2, c2 % L, (W - 2) // 2, c1 % L))

if PART == 'C':
    from flint import nmod_poly
    DISC = next((int(a[7:]) for a in sys.argv[1:] if a.startswith('--disc=')), -23)
    best = {}
    for f in sorted(glob.glob(os.path.join(HERE, 'level1_tower_poleHD%d_N*.json' % (-DISC)))):
        d = json.load(open(f, encoding='utf-8'))
        for Ls, rec in d['primes'].items():
            if Ls not in best or d['N'] > best[Ls][0]: best[Ls] = (d['N'], rec, os.path.basename(f), d['jpoly'])
    say('Part C: the towers level1_tower_poleHD%d_N…, discriminant %d' % (-DISC, DISC))
    for Ls, (Nn, rec, fn, H) in sorted(best.items(), key=lambda kv: int(kv[0])):
        L = int(Ls); W = -12 * (len(H) - 1)
        P, KP = family_P(W, L); S, KS = family_S(W, L)
        K = min(rec['K'], KP, KS if S else 99); M = L ** K; P = [c % M for c in P]; S = [c % M for c in S] if S else [1]
        lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}
        sym = leg(DISC, L)
        say(); say('ℓ = %d (%s, n < %d), lane ℓ^%d, weight %d: P on lane ℓ %s; slope-1 families at this weight %s' % (
            L, fn[:30] + '…', Nn, K, W, [c % L for c in P], len(S) - 1 if S != [1] else 'not computed'))
        if sym <= 0:
            say('   %s: P alone (PREDICTED): %s' % ('ramified' if sym == 0 else 'inert, the curves supersingular', show(*depths(lv, P, L, K))))
            continue
        fac = nmod_poly([c % L for c in H], L).factor()[1]; fdeg = sorted({g_.degree() for g_, e in fac})
        assert len(fdeg) == 1; fo = fdeg[0]; g = next(g_ for g_, e in fac)
        q = L ** fo; t = nmod_poly([0, 1], L)
        A = (t * (nmod_poly([1728], L) - t) * 3) % g; B = (t * (nmod_poly([1728], L) - t) * (nmod_poly([1728], L) - t) * 2) % g
        chi = 0
        import itertools
        for xs in itertools.product(range(L), repeat=fo):
            x = nmod_poly(list(xs), L)
            z = (x * x * x + A * x + B) % g
            if z == 0: continue
            chi += 1 if z.pow_mod((q - 1) // 2, g) == nmod_poly([1], L) else -1
        aq = -chi
        tn = next(x for x in range(1, 2 * q) for y in range(1, 2 * q) if x * x - DISC * y * y == 4 * q and (x % L or y % L))
        U = unit_root(aq, q, L, K); c = pow(U, W - 2, M)
        say('   split; the prime above ℓ has order %d: H mod ℓ factors into degree %d, so the curves %s; the curve over 𝔽_ℓ^%d has trace %d '
            '(norm form ±%d); U^%d digits %s' % (fo, fo, 'are each fixed' if fo == 1 else 'CYCLE with period %d' % fo, fo, aq, tn, W - 2, digits(c, L, K)))
        roots = []                                                   # the f-th roots of c in ℤ_ℓ (single factors x − r)
        for r0 in range(1, L):
            if (pow(r0, fo, L) - c) % L == 0 and (fo * pow(r0, fo - 1, L)) % L:
                r = r0
                for _ in range(4 * K + 8): r = (r - (pow(r, fo, M) - c) * pow(fo * pow(r, fo - 1, M) % M, -1, M)) % M
                roots.append(r)
        cyc = [(-c) % M] + [0] * (fo - 1) + [1]
        laws = [('P alone', P), ('P·S', pmul(P, S, M))]
        for r in roots: laws.append(('P·(x − r)·S, r ≡ %d one %s root' % (r % L, {2: 'square', 3: 'cube'}.get(fo, 'f-th')), pmul(pmul(P, [(-r) % M, 1], M), S, M)))
        laws += [('P·(x^%d − U^%d)' % (fo, W - 2), pmul(P, cyc, M)), ('P·(x^%d − U^%d)·S  PREDICTED' % (fo, W - 2), pmul(pmul(P, cyc, M), S, M))]
        for name, Q in laws: say('   %-40s %s' % (name, show(*depths(lv, Q, L, K))))
        G = {b: [sum(P[i] * lv[b + i][n] for i in range(len(P))) % M for n in range(min(len(lv[b + i]) for i in range(len(P))))]
             for b in lv if b + len(P) - 1 <= max(lv)}
        for step in range(1, fo + 1):
            reads = []
            for b in sorted(G):
                if b < 1 or b + step not in G: continue
                for n in range(1, min(len(G[b]), len(G[b + step]))):
                    v0 = val(G[b][n], L, K)
                    if v0 >= K - 1 or val(G[b + step][n], L, K) != v0: continue
                    mk = L ** (K - v0); reads.append((b, (G[b + step][n] // L ** v0) * pow((G[b][n] // L ** v0) % mk, -1, mk) % mk))
            say('   the ratio of the remainders %d level%s apart, on lane ℓ, by level: %s' % (step, '' if step == 1 else 's',
                {b: sorted({r % L for bb, r in reads if bb == b}) for b in sorted({bb for bb, _ in reads})}))
        say('   (on lane ℓ: U^%d ≡ %d)' % (W - 2, c % L))
