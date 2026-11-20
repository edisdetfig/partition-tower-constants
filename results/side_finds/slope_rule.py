# slope_rule.py, 14 Nov 2026 (README "Why exactly 13, 17, 19 and 23 settle fast", prediction (e)).
# For every prime 13 ≤ ℓ ≤ LMAX and every possible single seed σ = Sym^r ⊗ det^m on the row (r + 2m ≡ ℓ − 5 mod ℓ − 1,
# r even, 2 ≤ r ≤ ℓ − 1), from the representation theory of GL2(F_ℓ) alone:
#   d_z = how often σ occurs in Sym^(k_z − 2) (Glover's layers), D_z = how often in Ind_B^G Sym^(k_z − 2)
#   (serre_weights.py's factors), with k_z = −2 + (ℓ − 1)(z + 1);
#   for the Eisenstein type (r = ℓ − 3, its det^(m−1) carrying the numbers) the factor det^(m−1) is taken away, as in
#   serre_weights.py;
#   its zeros (Bergdall–Pollack's first coefficient): steps z with d_z = 0 and D_z ≥ 2; its slope at the tower's weight
#   −2 (step −1): Σ over zeros of 1 + v_ℓ(z + 1).
# Checked: the slope is min(m, ℓ − 4 − m) for the generic seeds with m ≤ ℓ − 5, and slope 1 occurs exactly for
# σ = Sym^(ℓ−7) ⊗ det and Sym⁴ ⊗ det^(ℓ−5). Then the families actually found (serre_weights_ℓ.txt) against the slopes
# counted at their arrival (ghost_series_ℓ_*.jsonl). DERIVED (representation theory), compared with COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/slope_rule.py [LMAX] [--keep]
import os, sys, re, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; LMAX = nums[0] if nums else 200; KEEP = '--keep' in sys.argv
OUT = os.path.join(HERE, 'slope_rule_%d.txt' % LMAX)
if KEEP and os.path.exists(OUT): sys.exit('stopped before computing anything: %s exists' % OUT)
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def is_prime(n): return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def vl(x, L):
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def jh_ind(s, t, LM):
    s, t = s % LM, t % LM
    if s == t: return [(0, t), (LM, t)]
    r = (s - t) % LM
    return [(r, t), (LM - r, s)]


def count_sym(n, sig, L):
    """How often sig occurs among the factors of Sym^n (Glover's recursion, layer by layer)."""
    LM = L - 1; c = 0; j = 0
    while True:
        if n <= LM: return c + ((n, j % LM) == sig)
        c += sum(1 for g in jh_ind(n, 0, LM) if (g[0], (g[1] + j) % LM) == sig)
        if n < L + 1: return c
        n -= L + 1; j += 1


def count_ind(k, sig, L):
    """How often sig occurs among the factors of Ind_B^G Sym^(k−2) = Σ_{i=0}^{k−2} Ind(a^(k−2−i) d^i)."""
    LM = L - 1; n = k - 2; c = 0
    full, rest = divmod(n + 1, LM)                      # i runs over full periods of LM, then a partial one
    per = [sum(1 for g in jh_ind(n - i, i, LM) if g == sig) for i in range(LM)]
    return full * sum(per) + sum(per[:rest])


def family_slope(L, sig, eis_m=None, zmax=None):
    """(arrival step, zeros, slope at weight −2) for a family whose one seed is sig; eis_m: the det^(m′) taken away."""
    zmax = zmax or L + 5                                 # every seed has arrived by layer ℓ − 2, so by step ℓ
    d, D = [], []
    for z in range(zmax + 1):
        k = -2 + (L - 1) * (z + 1)
        dz, Dz = count_sym(k - 2, sig, L), count_ind(k, sig, L)
        if eis_m is not None:
            dz -= count_sym(k - 2, (0, eis_m), L); Dz -= count_ind(k, (0, eis_m), L)
        d.append(dz); D.append(Dz)
    a = next((z for z in range(zmax + 1) if d[z] > 0), None)
    if a is None: return None, None, None
    zeros = [z for z in range(a) if D[z] >= 2]
    return a, zeros, sum(1 + vl(z + 1, L) for z in zeros)


say('slopes at the tower\'s weight of every single-seed family, from GL2(F_ℓ) alone, at the primes 13 … %d' % LMAX)
bad_formula = []; slope1 = {}; min_other = {}
for L in range(13, LMAX + 1):
    if not is_prime(L): continue
    LM = L - 1; ones = []
    for r in range(2, L, 2):
        for m in [m for m in range(LM) if (r + 2 * m - (L - 5)) % LM == 0]:
            kinds = [('generic' if r <= L - 3 else 'Steinberg', None)]
            if r == L - 3: kinds.append(('Eisenstein type', (m - 1) % LM))
            for kind, em in kinds:
                a, zeros, s = family_slope(L, (r, m), em)
                if s is None: continue
                if kind == 'generic' and m <= L - 5 and s != min(m, L - 4 - m): bad_formula.append((L, r, m, s))
                if s == 1: ones.append((kind, r, m))
                if not (kind == 'generic' and m in (1, L - 5)) and s >= 1:
                    min_other[L] = min(min_other.get(L, 99), s)
    slope1[L] = ones
    expect = sorted([('generic', L - 7, 1), ('generic', 4, L - 5)])
    if sorted(ones) != expect or L <= 31:
        say('   ℓ = %3d: slope 1 for %s%s; the smallest other positive slope: %s' % (
            L, ones, '' if sorted(ones) == expect else '   NOT the two first-step seeds', min_other.get(L)))
say('   slope = min(m, ℓ − 4 − m) for every generic seed with m ≤ ℓ − 5: %s' % ('yes' if not bad_formula else 'NO at %s' % bad_formula[:5]))
say('   slope 1 exactly for Sym^(ℓ−7) ⊗ det and Sym⁴ ⊗ det^(ℓ−5), at every prime: %s' % (
    'yes' if all(sorted(v) == sorted([('generic', L - 7, 1), ('generic', 4, L - 5)]) for L, v in slope1.items()) else 'NO'))
say('   the smallest positive slope of every other seed, over all primes: %d' % min(min_other.values()))

say(); say('the families found (seeds from serre_weights_ℓ.txt) against the slopes counted where they arrive (ghost_series)')
agree = total = 0
for L in (13, 17, 19, 23):
    txt = open(os.path.join(HERE, 'serre_weights_%d.txt' % L), encoding='utf-8').read()
    seeds = {}
    for lab, w, m in re.findall(r'system (\(.*?\))(?: \(Eisenstein type[^)]*\))?: seeds new weight (\d+)·ω\^(\d+) ×1', txt):
        seeds[lab] = (int(w) - 2, int(m), 'Eisenstein type' in txt.split('system ' + lab)[1].split('\n')[0])
    gs = [json.loads(t) for f in glob.glob(os.path.join(HERE, 'ghost_series_%d_x0-*.jsonl' % L)) for t in open(f, encoding='utf-8')]
    gs.sort(key=lambda r: r['x'])
    for lab, (r, m, eis) in sorted(seeds.items()):
        first = next((row for row in gs if row['systems'].get(lab, {}).get('d', 0) > 0), None)
        a, zeros, s = family_slope(L, (r, m), (m - 1) % (L - 1) if eis else None)
        if first is None: say('   ℓ = %d %s: seed Sym^%d ⊗ det^%d; predicted arrival %s, slope %s; not reached in the counted steps' % (L, lab, r, m, a, s)); continue
        cs = first['systems'][lab]['counted']; got = cs[0] if cs else None
        ok = (first['x'] == a and got is not None and int(got) == s); total += 1; agree += ok
        say('   ℓ = %d %s: seed Sym^%d ⊗ det^%d%s; predicted arrival x = %s, slope %s; counted x = %d, slope %s%s' % (
            L, lab, r, m, ' (Eisenstein type)' if eis else '', a, s, first['x'], got, '' if ok else '   DIFFERENT'))
say('   agreeing: %d of %d' % (agree, total))
