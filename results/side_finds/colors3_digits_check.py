# A dig, 12 Nov 2026: 3-COLORED PARTITIONS, THE FORMS' CONSTANTS AGAINST THE EXACT TOWERS. colors3_family_ℓ.txt gives, from
# modular forms of level 2 alone, the polynomial of u(−4) (one root for 13, 19, 23, 29, 41; two for 31, 37, 43, 47).
# colors3_tower_N150000000_*.json gives the exact towers H_b(n) = p₃(ℓᵇn + δ_b) on lane ℓ^K.
# PREDICTED before the run (derived, as for ordinary partitions):
#   - one form: the ratio H_(b+2)/H_b shares ever more digits with t·u(−4) as b grows (about one per level), t = ±1 the
#     twist, the same for every level and n;
#   - two forms: the law H_(b+4) + t·c₁·H_(b+2) + c₀·H_b ≡ 0 (x² + c₁x + c₀ the polynomial of u(−4)) holds on lane ℓ at
#     every level triple available (1, 3, 5), for one t and every n;
#   - the twists t(ℓ) are the values of one quadratic character.
# COUNTED.  usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors3_digits_check.py <tower json> [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = next(a for a in sys.argv[1:] if a.endswith('.json')); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'colors3_digits_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


data = json.load(open(os.path.join(HERE, SRC) if not os.path.isabs(SRC) else SRC, encoding='utf-8'))
say('towers: %s (n < %d)' % (os.path.basename(SRC), data['N']))
signs = {}
for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
    L = int(Ls); f = os.path.join(HERE, 'colors3_family_%d.txt' % L)
    if not os.path.exists(f): continue
    mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\]', open(f, encoding='utf-8').read())
    KF = int(mt.group(1)); poly = [int(v) for v in mt.group(2).split(',')]; d = len(poly) - 1
    KT = rec['K']; K = min(KF, KT); MK = L ** K
    lv = {int(b): [x % MK for x in v] for b, v in rec['levels'].items()}
    say()
    if d == 1:
        u = (-poly[0]) % MK; best = None
        for t in (1, -1):
            a = t * u % MK; rows = []; ok = True
            for b in sorted(lv):
                if b + 2 not in lv: continue
                n0 = next((i for i, x in enumerate(lv[b]) if x % L and i < len(lv[b + 2])), None)
                if n0 is None: continue
                r = lv[b + 2][n0] * pow(lv[b][n0], -1, MK) % MK
                same = next((e for e in range(K) if (r - a) % L ** (e + 1)), K)
                rows.append((b, n0, same))
            if rows and rows[0][2] >= 1 and (best is None or sum(s for _, _, s in rows) > sum(s for _, _, s in best[1])): best = (t, rows)
        if best is None: say('ℓ = %d (one form): no twist gives the first digit; NOT as predicted' % L); continue
        t, rows = best; signs[L] = t
        say('ℓ = %d (one form), twist t = %+d, t·u(−4) = %s…: digits shared with H_(b+2)/H_b at (b, b+2): %s' % (
            L, t, [((t * u % MK) // L ** e) % L for e in range(min(K, 8))], ', '.join('(%d,%d) n=%d: %d' % (b, b + 2, n0, s) for b, n0, s in rows)))
    else:
        c0, c1 = poly[0] % MK, poly[1] % MK; res = {}
        for t in (1, -1):
            out = []
            for b in sorted(lv):
                if b + 4 not in lv: continue
                for n, (x, y, z) in enumerate(zip(lv[b], lv[b + 2], lv[b + 4])):
                    w = (z + t * c1 * y + c0 * x) % MK
                    out.append((b, n, next((e for e in range(K) if w % L ** (e + 1)), K)))
            res[t] = out
        good = [t for t in (1, -1) if res[t] and all(v >= 1 for _, _, v in res[t])]
        signs[L] = good[0] if len(good) == 1 else None
        for t in (1, -1):
            say('ℓ = %d (two forms), t = %+d: the law\'s ℓ-adic depth at (b, b+2, b+4), n: %s' % (
                L, t, ', '.join('(%d) n=%d: %d' % (b, n, v) for b, n, v in res[t]) or 'no level triple below N'))
say()
say('the twists t(ℓ): %s' % signs)
for name, a in (('1', 1), ('(−1|ℓ)', -1), ('(2|ℓ)', 2), ('(−2|ℓ)', -2), ('(3|ℓ)', 3), ('(−3|ℓ)', -3), ('(6|ℓ)', 6), ('(−6|ℓ)', -6)):
    for sg, lab in ((1, '+'), (-1, '−')):
        fit = [L for L, t in signs.items() if t is not None and t == sg * (1 if a == 1 else leg(a, L))]
        n_ = sum(1 for t in signs.values() if t is not None)
        if len(fit) == n_: say('   t = %s%s fits all %d' % (lab, name, n_))
