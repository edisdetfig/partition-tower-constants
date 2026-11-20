# A dig, 12 Nov 2026: ARE THE LOOPS RANDOM? Reads partitions_loops_survey_13_700.txt (the constants of every prime's tower
# on lane ℓ, and their orders) and compares with what random units would give (DERIVED):
#   - how often a constant has each "index" (ℓ − 1)/order: chance φ((ℓ − 1)/i)/(ℓ − 1) for index i;
#   - how often a constant is a small special number (±1, ±2, ±3, ±6 and their inverses): chance 1/(ℓ − 1) each;
#   - for towers with two constants on lane ℓ, how often the two have the same order.
# PREDICTED before the run: everything within chance (the run of primitive roots from 7 to 29, and 6 for 13, 17 and 43,
# luck). COUNTED (the data), DERIVED (the expectations).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_loops_stats.py [--keep]
import os, sys, re, math
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_loops_stats.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def phi(n):
    r, p = n, 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0: n //= p
            r -= r // p
        p += 1
    if n > 1: r -= r // n
    return r


text = open(os.path.join(HERE, 'partitions_loops_survey_13_700.txt'), encoding='utf-8').read()
consts = []                                                    # (ℓ, value, order)
towers = []
for line in text.split('\n'):
    mm = re.match(r'\s+ℓ =\s+(\d+):', line)
    if not mm: continue
    L = int(mm.group(1)); here = [(L, int(v), int(o)) for v, o in re.findall(r'(\d+) \(order (\d+) of \d+', line)]
    consts += here; towers.append((L, here))
say('constants on lane ℓ: %d, from %d primes (13 to 700)' % (len(consts), len(towers)))

say(); say('1  THE INDEX (ℓ − 1)/order (1 = a primitive root, the longest loop)')
for i in (1, 2, 3, 4, 5, 6):
    cnt = sum(1 for L, v, o in consts if (L - 1) // o == i)
    exp = sum(phi((L - 1) // i) / (L - 1) for L, v, o in consts if (L - 1) % i == 0)
    say('   index %d: counted %3d, expected %6.1f (derived)' % (i, cnt, exp))
big = sum(1 for L, v, o in consts if (L - 1) // o > 6)
expb = sum(1 - sum(phi((L - 1) // i) / (L - 1) for i in range(1, 7) if (L - 1) % i == 0) for L, v, o in consts)
say('   index above 6: counted %3d, expected %6.1f' % (big, expb))

say(); say('2  SMALL SPECIAL NUMBERS (each has chance 1/(ℓ − 1))')
for name, f in (('1', lambda L: 1), ('−1', lambda L: L - 1), ('2', lambda L: 2), ('−2', lambda L: L - 2), ('3', lambda L: 3),
                ('−3', lambda L: L - 3), ('6', lambda L: 6), ('−6', lambda L: L - 6), ('1/6', lambda L: pow(6, -1, L)),
                ('−1/6', lambda L: (-pow(6, -1, L)) % L), ('1/2', lambda L: pow(2, -1, L)), ('1/3', lambda L: pow(3, -1, L))):
    hits = [L for L, v, o in consts if v == f(L) % L]
    exp = sum(1 / (L - 1) for L, v, o in consts)
    say('   %-5s counted %2d %s, expected %.2f' % (name, len(hits), hits, exp))

say(); say('3  TWO CONSTANTS ON LANE ℓ: THE SAME ORDER?')
pairs = [(L, c) for L, c in towers if len(c) == 2]
same = [L for L, c in pairs if c[0][2] == c[1][2]]
exps = 0.0
for L, c in pairs:                                             # chance that two independent random units share their order
    exps += sum((phi(dv) / (L - 1)) ** 2 for dv in range(1, L) if (L - 1) % dv == 0)
say('   towers with exactly two constants on lane ℓ: %d; the same order: %d %s; expected %.1f (derived)' % (len(pairs), len(same), same, exps))
