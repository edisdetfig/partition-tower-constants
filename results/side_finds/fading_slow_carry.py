# A dig, 13 Nov 2026: THE SLOW PRIMES' QUIET NOTES CARRIED THROUGH THREE WEIGHTS (README "Counted at 29 and 31"; predicted
# there first). The runs of fading_forms3.py at 29 and 31 were stopped by low memory at their third weight; the first two
# weights are read from the kept outputs fading_forms3_ℓ_s1_split.txt, the third from the resumed run
# fading_forms3_ℓ_s1_split_N8_X2_from2.txt (--from=2). The slope-1 root's (3|ℓ)·u/ℓ, 8 digits printed per weight, is
# carried to weight −2 by Lagrange at x = −1: through 2 weights 2·w(0) − w(1), through 3 weights 3·w(0) − 3·w(1) + w(2), and
# compared with the fading notes measured on partition numbers (tower_fading.txt). DERIVED from the kept outputs.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/fading_slow_carry.py [--keep]
import os, sys, re
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'fading_slow_carry.txt'), 'x', encoding='utf-8') if KEEP else None
MEASURED = {29: [23, 3, 10], 31: [2, 28, 13]}


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def roots(name):
    """(weight, system, digits) for each slope-1 root line, in file order."""
    txt = open(os.path.join(HERE, name), encoding='utf-8').read(); out = []
    for m in re.finditer(r'the system (\{[^}]*\}): \d+ forms?; slopes of T_ℓ: [^\n]*?1 root of slope 1[^\n]*\n\s*the root with u/ℓ\^1 ≡ \d+: '
                         r'known to ℓ\^(\d+); \(3\|ℓ\)·u/ℓ\^1, digits lowest first: \[([^\]]*)\]', txt):
        out.append((m.group(1), int(m.group(2)), [int(x) for x in m.group(3).split(',')]))
    ws = [int(w) for w in re.findall(r'\n   weight (\d+) \(dimension', txt)]
    return ws, out


for L in (29, 31):
    w01, r01 = roots('fading_forms3_%d_s1_split.txt' % L)
    w2, r2 = roots('fading_forms3_%d_s1_split_N8_X2_from2.txt' % L)
    rows = r01[:2] + r2[:1]; ws = w01[:2] + w2[:1]
    say(); say('ℓ = %d: weights %s; the slope-1 root\'s system at each: %s; known to ℓ^%s' % (L, ws, [r[0] for r in rows], [r[1] for r in rows]))
    assert len(rows) == 3 and len({r[0] for r in rows}) == 1, 'not one system with one slope-1 root at each weight'
    m = L ** 8; w = [sum(d * L ** i for i, d in enumerate(r[2])) for r in rows]
    for n, v in ((2, 2 * w[0] - w[1]), (3, 3 * w[0] - 3 * w[1] + w[2])):
        dg = [((v % m) // L ** i) % L for i in range(8)]
        same = next((i for i in range(len(MEASURED[L])) if dg[i] != MEASURED[L][i]), len(MEASURED[L]))
        say('   through %d weights: (3|ℓ)·u(−2)/ℓ = %s …;  measured on partition numbers: %s  →  %d of %d digits agree' % (
            n, dg[:5], MEASURED[L], same, len(MEASURED[L])))
