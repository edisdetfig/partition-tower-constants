# A dig, 12 Nov 2026: r-COLORED PARTITIONS, THE FORMS AGAINST THE EXACT TOWERS (any r; colors3_compare.py and
# colors3_digits_check.py did r = 3 on level 2). Reads a tower file colors_tower_r*.json (exact values on lane ℓ^K) and the
# family files (colors6_family_r<r>_ℓ.txt on level 6, colors3_family_ℓ.txt on level 2).
# For each prime: the rank R_ℓ(r) (Belmont–Lee–Musat–Trebat-Leder, QUOTED formula); for R = 0 the laws (the smallest
# valuation at each level); else the class on lane ℓ, T_ℓ on it, and
#   - one form: one constant for every level pair and n? its first digit against a_ℓ(f): the twist t = α/a_ℓ(f) (±1?);
#     the digits shared with t·u(−(r+1)) at each level pair;
#   - two forms: H_(b+4) + t·c₁·H_(b+2) + c₀·H_b on lane ℓ (the class's numbers) and its depth with the family's polynomial;
# and at the end the one quadratic character that gives every twist.
# PREDICTED for 5 colours before the run (derived): everything as for 1 and 3 colours, with the twist t(ℓ) = (3|ℓ), the
# character of η⁵ as of η (3 colours had (−1|ℓ) because η³ is a theta series of 8). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors_check.py <tower json> --level=6 --e2=-1 --e3=-1 [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import forms6_flint as F6
import forms2_flint as F2
SRC = next(a for a in sys.argv[1:] if a.endswith('.json')); KEEP = '--keep' in sys.argv
opt = lambda name, dflt: next((int(a[len(name) + 3:]) for a in sys.argv[1:] if a.startswith('--%s=' % name)), dflt)
LEVEL, E2, E3 = opt('level', 6), opt('e2', -1), opt('e3', -1)
data = json.load(open(os.path.join(HERE, SRC) if not os.path.isabs(SRC) else SRC, encoding='utf-8')); RC = data['r']
fh = open(os.path.join(HERE, 'colors_check_%s.txt' % os.path.basename(SRC)[len('colors_tower_'):-5]), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def R(L, r):
    k = ((r + 1) // 2) * (L - 1)
    return k // 12 - (1 if k % 12 == 2 else 0) - (r * (L * L - 1)) // (24 * L)


say('%d colours; towers %s (n < %d); level %d' % (RC, os.path.basename(SRC), data['N'], LEVEL))
signs = {}
for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
    L = int(Ls); K = rec['K']; MK = L ** K; lv = {int(b): [x % MK for x in v] for b, v in rec['levels'].items()}; rr = R(L, RC)
    val = lambda x: next((e for e in range(K) if x % L ** (e + 1)), K)
    say()
    if rr == 0:
        say('ℓ = %d (R = 0): the smallest ℓ-adic valuation at each level b: %s' % (L, {b: min(val(x) for x in v) for b, v in lv.items()}))
        continue
    k = L - RC - 2
    if LEVEL == 6:
        s2, s3 = E2 * leg(2, L), E3 * leg(3, L)
        H = F6.hecke_matrices(k, L ** 4, [('U', 2), ('U', 3), ('T', L)])
        W_, RW_, r_, Wn, Rn, d, Np = F6.class_padic(H, k, L, 4, s2, s3)
        if Wn is None: say('ℓ = %d: the class is empty on lane ℓ' % L); continue
        Xc, Rc = F6.reduce_mat(Wn, L), Rn; Tl = F6.restrict(F6.reduce_mat(H[('T', L)], L), Xc, Rc, d)
        fam = os.path.join(HERE, 'colors6_family_r%d_%d.txt' % (RC, L))
    else:
        Xc, Rc, d = F2.class_on_lane2(k, L, leg(2, L))
        Tl = F6.restrict(F2.hecke_matrices2(k, L, [('T', L)])[('T', L)], Xc, Rc, d)
        fam = os.path.join(HERE, 'colors3_family_%d.txt' % L)
    cp = [int(c) for c in Tl.charpoly().coeffs()]
    say('ℓ = %d (R = %d): the class on lane ℓ has %d form%s (weight %d); T_ℓ polynomial %s' % (L, rr, d, '' if d == 1 else 's', k, cp))
    poly = None
    if os.path.exists(fam):
        mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\]', open(fam, encoding='utf-8').read())
        if mt: KF = int(mt.group(1)); poly = [int(v) for v in mt.group(2).split(',')]
    if d == 1:
        a = (-cp[0]) % L; ratios = []
        for b in lv:
            if b + 2 in lv:
                ratios += [lv[b + 2][i] * pow(x, -1, L) % L for i, x in enumerate(lv[b]) if x % L and i < len(lv[b + 2])]
        al = ratios[0] if ratios else None; t = al * pow(a, -1, L) % L if al is not None else None
        tt = 1 if t == 1 else (-1 if t == L - 1 else None); signs[L] = tt
        say('   one constant for all %d (b, n): %s;  α = %s, a_ℓ(f) = %d  →  %s' % (
            len(ratios), 'yes' if len(set(ratios)) == 1 else 'NO', al, a, ('t = %+d' % tt) if tt else 'NOT ±1'))
        if poly and tt:
            KK = min(K, KF); M2 = L ** KK; u = tt * (-poly[0]) % M2; rows = []
            for b in sorted(lv):
                if b + 2 not in lv: continue
                n0 = next((i for i, x in enumerate(lv[b]) if x % L and i < len(lv[b + 2])), None)
                if n0 is None: continue
                rq = lv[b + 2][n0] * pow(lv[b][n0], -1, M2) % M2
                rows.append('(%d,%d): %d' % (b, b + 2, next((e for e in range(KK) if (rq - u) % L ** (e + 1)), KK)))
            say('   t·u(−%d) = %s…; digits shared with H_(b+2)/H_b: %s' % (RC + 1, [(u // L ** e) % L for e in range(min(KK, 8))], ', '.join(rows)))
    else:
        res = {}
        for t in (1, -1):
            T = t * (-cp[1]) % L; Nn = cp[0] % L; okc = tot = 0
            for b in lv:
                if b + 4 not in lv: continue
                for x, y, z in zip(lv[b], lv[b + 2], lv[b + 4]): tot += 1; okc += (z - T * y + Nn * x) % L == 0
            res[t] = (okc, tot)
        good = [t for t in (1, -1) if res[t][1] and res[t][0] == res[t][1]]
        signs[L] = good[0] if len(good) == 1 else None
        say('   the three-level law on lane ℓ: t = +1: %d of %d;  t = −1: %d of %d' % (res[1][0], res[1][1], res[-1][0], res[-1][1]))
        if poly and len(poly) == 3:
            KK = min(K, KF); M2 = L ** KK; c0, c1 = poly[0] % M2, poly[1] % M2
            for t in (1, -1):
                depths = []
                for b in sorted(lv):
                    if b + 4 not in lv: continue
                    for n, (x, y, z) in enumerate(zip(lv[b], lv[b + 2], lv[b + 4])):
                        w = (z + t * c1 * y + c0 * x) % M2
                        depths.append('(%d) n=%d: %d' % (b, n, next((e for e in range(KK) if w % L ** (e + 1)), KK)))
                say('   with the family\'s polynomial, t = %+d: depth %s' % (t, ', '.join(depths) or '—'))
say()
say('the twists t(ℓ): %s' % signs)
n_ = sum(1 for t in signs.values() if t is not None)
for name, a in (('1', 1), ('(−1|ℓ)', -1), ('(2|ℓ)', 2), ('(−2|ℓ)', -2), ('(3|ℓ)', 3), ('(−3|ℓ)', -3), ('(6|ℓ)', 6), ('(−6|ℓ)', -6)):
    for sg, lab in ((1, '+'), (-1, '−')):
        fit = sum(1 for L, t in signs.items() if t is not None and t == sg * (1 if a == 1 else leg(a, L)))
        say('   t = %s%s fits %d of %d' % (lab, name, fit, n_))
