# A dig, 12 Nov 2026: 3-COLORED PARTITIONS, THE LAWS AT 11 AND 17 AGAINST DEEPER EXACT TOWERS. From the forms alone
# (colors3_slopes_11.txt, _17.txt): 17 follows one family of slope 1, constant (−1|17)·u(−4)/17 = 4, 12, 7, 2, 10, 4, 8, …;
# 11 follows a PAIR of slope-2 families whose constants agree on lane 11 (a repeated root), so its tower obeys the three-term
# law H′_(b+4) + c₁·H′_(b+2) + c₀·H′_b ≡ 0, x² + c₁x + c₀ the pair's polynomial (to 11⁸), and no single constant.
# The towers: colors_tower_r3_N100000000_11_17.json (exact, on 60-bit lanes, WSL); H′ = H divided by the laws' power
# (11ᵇ; 17^⌈b/2⌉).
# PREDICTED before the run (derived): for 17 the ratio H′_(b+2)/H′_b shares ever more digits with the forms' constant (about
# one per two levels, as Watson's 7 does for ordinary partitions); for 11 the three-term law holds ever deeper as b grows.
# COUNTED.  usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors3_laws_check.py [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'colors3_laws_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


data = json.load(open(os.path.join(HERE, 'colors_tower_r3_N100000000_11_17.json'), encoding='utf-8'))
say('towers: colors_tower_r3_N100000000_11_17.json (n < %d)' % data['N'])
LAW = {11: lambda b: b, 17: lambda b: (b + 1) // 2}
for L in (17, 11):
    rec = data['primes'][str(L)]; K = rec['K']; lv = {int(b): v for b, v in rec['levels'].items()}
    norm = {b: [(x // L ** LAW[L](b)) for x in v] for b, v in lv.items() if all(x % L ** LAW[L](b) == 0 for x in v)}
    prec = {b: K - LAW[L](b) for b in norm}
    say(); say('ℓ = %d: levels %s; every value divisible by the law\'s power: %s' % (L, sorted(lv), sorted(norm) == sorted(lv)))
    text = open(os.path.join(HERE, 'colors3_slopes_%d.txt' % L), encoding='utf-8').read()
    if L == 17:
        cst = [int(v) for v in re.search(r'digits lowest first[^:]*: \[([^\]]*)\]', text).group(1).split(',')]
        reach = int(re.search(r'the series reaches (\d+)', text).group(1)); cst = cst[:reach]
        a = sum(dg * L ** j for j, dg in enumerate(cst)); say('   the forms\' constant: %s (the series\' reach)' % cst)
        for b in sorted(norm):
            if b + 2 not in norm: continue
            P = min(prec[b], prec[b + 2], len(cst)); M = L ** P; rows = []
            for n in range(min(len(norm[b]), len(norm[b + 2]))):
                x, y = norm[b][n] % M, norm[b + 2][n] % M
                if x % L == 0: continue
                rr = y * pow(x, -1, M) % M
                rows.append(next((e for e in range(P) if (rr - a) % L ** (e + 1)), P))
            say('   (%d,%d): digits shared with the forms\' constant, n = 0..: %s (of %d available)' % (b, b + 2, rows, P))
    else:
        mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\]', text); KF = int(mt.group(1)); c0, c1 = [int(v) for v in mt.group(2).split(',')][:2]
        say('   the pair\'s polynomial (forms, lane 11^%d): x² + %d·x + %d' % (KF, c1, c0))
        for b in sorted(norm):
            if b + 4 not in norm: continue
            P = min(prec[b], prec[b + 2], prec[b + 4], KF); M = L ** P; rows = []
            for n in range(min(len(norm[b]), len(norm[b + 2]), len(norm[b + 4]))):
                w = (norm[b + 4][n] + c1 * norm[b + 2][n] + c0 * norm[b][n]) % M
                rows.append(next((e for e in range(P) if w % L ** (e + 1)), P))
            say('   (%d,%d,%d): the law holds to 11^…, n = 0..: %s (of %d available)' % (b, b + 2, b + 4, rows, P))
        for b in sorted(norm):                                  # a single constant, for contrast
            if b + 2 not in norm: continue
            P = min(prec[b], prec[b + 2]); M = L ** P
            firsts = sorted({norm[b + 2][n] * pow(norm[b][n], -1, L) % L for n in range(min(len(norm[b]), len(norm[b + 2]))) if norm[b][n] % L})
            say('   (%d,%d): H′_(b+2)/H′_b on lane 11 takes the values %s (one constant would give one)' % (b, b + 2, firsts))
