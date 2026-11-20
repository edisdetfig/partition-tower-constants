# A dig, 12 Nov 2026: 3-COLORED PARTITIONS, THE FIRST DIGITS AGAINST THE FORMS. The towers H_b(n) = p₃(ℓᵇn + δ_b) come from
# colors3_tower_*.json (exact values on lane ℓ^K, WSL). The forms: the class of colors3_census.py, the newforms of weight
# ℓ − 5 on Γ0(2) with U₂ = (2|ℓ)·2^((ℓ−7)/2) (exactly, forms2_flint), and T_ℓ on it on lane ℓ.
# PREDICTED before the run (derived, as for ordinary partitions, where α ≡ (3|ℓ)·a_ℓ(f)):
#   - one form (R = 1): one constant α on lane ℓ for every n and every level pair, and α ≡ t·a_ℓ(f) with t = ±1;
#   - two forms (R = 2): H_(b+4) ≡ T·H_(b+2) − N·H_b on lane ℓ with T = t·(a_ℓ(f) + a_ℓ(f′)), N = a_ℓ(f)·a_ℓ(f′), for one
#     sign t; no single constant fits;
#   - the signs t(ℓ) are the values at ℓ of ONE quadratic character (read off, then checked on every prime);
#   - no form (R = 0: 11, 17): laws, H_b ≡ 0 modulo a power of ℓ that grows with b (BLMT's Corollary 1.2).
# COUNTED (agreements), with the character read off the data.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors3_compare.py <tower json> [--keep]
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import nmod_mat
import forms2_flint as F2
import forms6_flint as F6
SRC = next(a for a in sys.argv[1:] if a.endswith('.json')); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'colors3_compare.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def R(L):
    k = 2 * (L - 1)
    return k // 12 - (1 if k % 12 == 2 else 0) - (3 * (L * L - 1)) // (24 * L)


data = json.load(open(os.path.join(HERE, SRC) if not os.path.isabs(SRC) else SRC, encoding='utf-8'))
say('towers from %s (n < %d)' % (os.path.basename(SRC), data['N']))
signs = {}
for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
    L = int(Ls); K = rec['K']; lv = {int(b): v for b, v in rec['levels'].items()}; r = R(L); say()
    val = lambda x: next((e for e in range(K) if x % L ** (e + 1)), K)
    if r == 0:
        say('ℓ = %d (R = 0): the smallest ℓ-adic valuation at each level b: %s' % (L, {b: min(val(x) for x in v) for b, v in lv.items()}))
        continue
    k = L - 5; s = leg(2, L); dz = F2.class_basis_exact2(k, s)[1]
    Xc, Rc, d = F2.class_on_lane2(k, L, s)                      # the class on lane ℓ (exactly separated)
    if d != dz: say('ℓ = %d: the class has %d forms exactly but %d on lane ℓ; skipped' % (L, dz, d)); continue
    H = F2.hecke_matrices2(k, L, [('T', L)]); Tl = F6.restrict(H[('T', L)], Xc, Rc, d)
    cp = [int(c) for c in Tl.charpoly().coeffs()]               # lowest first
    roots = [int(x) for x, e in Tl.charpoly().roots()]
    say('ℓ = %d (R = %d): the class (weight %d, level 2, U₂ sign %+d) has %d form%s; T_ℓ on lane ℓ: polynomial %s, roots on lane ℓ %s' % (
        L, r, k, s, d, '' if d == 1 else 's', cp, roots))
    if d != r: say('   the class has %d forms, R = %d: NOT as the census' % (d, r))
    if d == 1:
        a = (-cp[0]) % L; al = None; seen = []
        for b in lv:
            if b + 2 not in lv: continue
            for x, y in zip(lv[b], lv[b + 2]):
                if x % L: seen.append(y * pow(x, -1, L) % L)
        one = len(set(seen)) == 1
        al = seen[0] if seen else None
        t = None if al is None else (al * pow(a, -1, L) % L)
        tt = 1 if t == 1 else (-1 if t == L - 1 else None)
        signs[L] = tt
        say('   the tower: H_(b+2)/H_b on lane ℓ over %d (b, n): %s  →  %s;  α = %s, a_ℓ(f) = %d, α/a_ℓ(f) = %s  →  %s' % (
            len(seen), sorted(set(seen)), 'ONE constant' if one else 'NOT one constant', al, a, t,
            ('t = %+d' % tt) if tt else 'NOT ±1'))
        # the constant's digits from the deepest pair
        b2 = max(b for b in lv if b - 2 in lv); b1 = b2 - 2; MK = L ** K
        n0 = next((i for i, x in enumerate(lv[b1]) if x % L), None)
        if n0 is not None and n0 < len(lv[b2]):
            rr = lv[b2][n0] * pow(lv[b1][n0], -1, MK) % MK
            say('   the deepest pair (%d, %d), n = %d: H_%d/H_%d, digits lowest first: %s' % (b1, b2, n0, b2, b1, [(rr // L ** j) % L for j in range(6)]))
    elif d == 2:
        res = {}
        for t in (1, -1):
            T = t * (-cp[1]) % L; Nn = cp[0] % L; ok = tot = 0
            for b in lv:
                if b + 4 not in lv: continue
                for x, y, z in zip(lv[b], lv[b + 2], lv[b + 4]):
                    tot += 1; ok += (z - T * y + Nn * x) % L == 0
            res[t] = (ok, tot, T, Nn)
        single = 0; tot1 = 0
        for b in lv:
            if b + 2 not in lv: continue
            pairs = [(x % L, y % L) for x, y in zip(lv[b], lv[b + 2]) if x % L]
            tot1 += 1; single += len({y * pow(x, -1, L) % L for x, y in pairs}) <= 1
        good = [t for t in (1, -1) if res[t][1] and res[t][0] == res[t][1]]
        signs[L] = good[0] if len(good) == 1 else None
        say('   the three-level law on lane ℓ: t = +1 (T, N) = (%d, %d): %d of %d;  t = −1 (T, N) = (%d, %d): %d of %d;  one constant for all n: %d of %d level pairs' % (
            res[1][2], res[1][3], res[1][0], res[1][1], res[-1][2], res[-1][3], res[-1][0], res[-1][1], single, tot1))
say()
say('the signs t(ℓ): %s' % signs)
for name, a in (('1', 1), ('(−1|ℓ)', -1), ('(2|ℓ)', 2), ('(−2|ℓ)', -2), ('(3|ℓ)', 3), ('(−3|ℓ)', -3), ('(6|ℓ)', 6), ('(−6|ℓ)', -6)):
    fit = [L for L, t in signs.items() if t is not None and (t == (1 if a == 1 else leg(a, L)))]
    sfit = [L for L, t in signs.items() if t is not None and (t == -(1 if a == 1 else leg(a, L)))]
    say('   t = +%s fits %d of %d;  t = −%s fits %d of %d' % (name, len(fit), sum(1 for t in signs.values() if t), name, len(sfit), sum(1 for t in signs.values() if t)))
