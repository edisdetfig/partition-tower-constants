# A dig, 12 Nov 2026: 5 COLOURS AT 23, BELMONT–LEE–MUSAT–TREBAT-LEDER'S "SUPEREXCEPTIONAL" PAIR (5, 23). Every value of the
# exact tower p₅(23ᵇn + δ_b) is divisible by 23 exactly once at every level (colors_tower_r5_N120000000_*.txt), yet the
# rank is 1 (QUOTED: their Remark 2 names (5, 23) and (7, 19) superexceptional, rank 1). So the tower should be 23 times a
# tower with one constant, (3|23)·u(−6) = u(−6) from the family (colors6_family_r5_23.txt).
# PREDICTED before the run: after dividing by 23, one constant on lane 23 from some level on, and its digits agree with
# u(−6) as far as the towers settle. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors5_check23.py [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'colors5_check23.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


d = json.load(open(os.path.join(HERE, 'colors_tower_r5_N120000000_13_17_19_23_29_31_11_37_41.json'), encoding='utf-8'))
L = 23; rec = d['primes']['23']; K = rec['K']; lv = {int(b): v for b, v in rec['levels'].items()}; M = L ** (K - 1)
fam = open(os.path.join(HERE, 'colors6_family_r5_23.txt'), encoding='utf-8').read()
u = [int(v) for v in re.search(r'u\(−6\), digits lowest first: \[([^\]]*)\]', fam).group(1).split(',')]
a = sum(dg * L ** j for j, dg in enumerate(u[:K - 1]))
say('the forms: (3|23)·u(−6) = u(−6) = %s…' % u[:6])
say('every tower value divisible by 23: %s' % all(x % L == 0 for v in lv.values() for x in v))
for b in sorted(lv):
    if b + 2 not in lv: continue
    firsts = set(); shared = []
    for n in range(min(len(lv[b]), len(lv[b + 2]))):
        x, y = (lv[b][n] // L) % M, (lv[b + 2][n] // L) % M
        if x % L == 0: continue
        rr = y * pow(x, -1, M) % M; firsts.add(rr % L)
        shared.append(next((e for e in range(K - 1) if (rr - a) % L ** (e + 1)), K - 1))
    say('   (%d,%d): (H/23)_(b+2)/(H/23)_b on lane 23 takes %s; digits shared with u(−6), n = 0..: %s' % (b, b + 2, sorted(firsts), shared))
