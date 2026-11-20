# A dig, 12 Nov 2026: 24 COLOURS AT 43, A 607-TYPE EXCEPTION PREDICTED BY THE FORMS. forms1_family_w-12_43.txt: on lane 43
# the constants are 1, 30, 30 — the two cusp forms of weight 30 (a conjugate pair, Hecke field Q(√51349), 51349 a non-square
# mod 43, so inert) share their number at 43, as 607's pair did (README, 607 one digit deeper).
# PREDICTED before the run: if T_43 acts on the pair as one number on lane 43, the tower fills only 2 directions there, so
# the SHORTER law H_(b+2) − 31·H_(b+1) + 30·H_b ≡ 0 (constants 1 and 30) holds on lane 43 for every n and level, while at 41
# (constants 1, 33, 34, all different) no law of degree 2 built from two of its constants holds. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors24_check43.py [--keep]
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'colors24_check43.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


d = json.load(open(os.path.join(HERE, 'colors_tower_r24_N100000000_29_31_37_41_43.json'), encoding='utf-8'))
for L, pairs in ((43, [(1, 30)]), (41, [(1, 33), (1, 34), (33, 34)])):
    lv = {int(b): v for b, v in d['primes'][str(L)]['levels'].items()}
    for a, c in pairs:
        s, p = (a + c) % L, (a * c) % L; res = []
        for b in sorted(lv):
            if b + 2 not in lv: continue
            nn = min(len(lv[b]), len(lv[b + 1]), len(lv[b + 2]))
            res.append((b, [(lv[b + 2][n] - s * lv[b + 1][n] + p * lv[b][n]) % L for n in range(nn)]))
        ok = all(all(x == 0 for x in r) for _, r in res)
        say('ℓ = %d, the law with the constants %d and %d (x² − %d·x + %d): %s  →  %s' % (
            L, a, c, s, p, '; '.join('from level %d: %s' % (b, r) for b, r in res), 'holds' if ok else 'does NOT hold'))
