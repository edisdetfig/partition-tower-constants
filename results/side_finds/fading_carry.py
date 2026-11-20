# A dig, 13 Nov 2026: THE FADING FORM CARRIED TO WEIGHT −2 (README "The fading form"). fading_forms2.py --sys carried its
# slope-2 root to weight −2 through every weight, the last ones known to only 4, 2 or 0 digits, which capped its answer at 2
# digits. Here only the weights known to 8 digits or more are used (Lagrange interpolation at x = −1 through x = 0 … n − 1,
# n = 2 … 5), read from the kept outputs fading_forms2_13_s2_sys.txt and fading_forms2_17_s2_sys.txt (8 digits printed per
# weight), and compared with the fading constants measured on partition numbers (tower_fading.txt). DERIVED from the kept
# outputs.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/fading_carry.py [--keep]
import os, sys, re
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'fading_carry.txt'), 'x', encoding='utf-8') if KEEP else None
MEASURED = {13: [10, 11, 6, 2, 10, 10], 17: [3, 9, 1, 13, 4, 3]}


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def lagrange_at_minus1(vals, L, P):
    m = L ** P; v = 0
    for i in range(len(vals)):
        num = den = 1
        for j in range(len(vals)):
            if j != i: num *= (-1 - j); den *= (i - j)
        v = (v + vals[i] * num * pow(den % m, -1, m)) % m
    return [(v // L ** k) % L for k in range(P)]


for L in (13, 17):
    txt = open(os.path.join(HERE, 'fading_forms2_%d_s2_sys.txt' % L), encoding='utf-8').read()
    rows = []
    for k, known, dg in re.findall(r'weight (\d+) \(dimension \d+\).*?known to ℓ\^(\d+); \(3\|ℓ\)·u/ℓ\^2, digits lowest first: \[([^\]]*)\]', txt):
        if int(known) >= 8: rows.append((int(k), [int(x) for x in dg.split(',')]))
    say(); say('ℓ = %d: weights known to 8 digits or more: %s' % (L, [k for k, _ in rows]))
    vals = [sum(d * L ** i for i, d in enumerate(dg)) for _, dg in rows]
    for n in range(2, len(rows) + 1):
        dg = lagrange_at_minus1(vals[:n], L, 8)
        agree = next((i for i in range(len(MEASURED[L])) if dg[i] != MEASURED[L][i]), len(MEASURED[L]))
        say('   through %d weights: (3|ℓ)·u(−2)/ℓ² = %s …;  measured on partition numbers: %s  →  %d of %d digits agree' % (
            n, dg[:6], MEASURED[L], agree, len(MEASURED[L])))
