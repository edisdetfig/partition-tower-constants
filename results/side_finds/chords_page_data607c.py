# A dig, 13 Nov 2026: THE 607 PAGE'S SUPPORTING DATA, TO BE CARRIED IN THE PAGE ITSELF. Read from the
# kept outputs, nothing recomputed:
#   - partitions_flint_exceptions_new_37_1297.txt: for each of the 200 primes, the class's size on lane ℓ, ⌊(ℓ + 12)/24⌋,
#     the exact class (499: 21 of 22, one congruent copy of lower level), zero roots and repeated factors of T_ℓ;
#   - partitions_flint_607_deeper.txt: at weights 604, 1210, 1816, the class's size and its ordinary part, T₅'s and C's
#     polynomials on the pair, form A's δ; the second difference; δ(−2); the tower's constants modulo 607²;
#   - partitions_flint_607_span_599_601_607_613_499.txt: the degree of T₆₀₇'s minimal polynomial on the class.
# Writes chords_page_data607c.json. COUNTED (read from the kept outputs).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/chords_page_data607c.py [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
rd = lambda n: open(os.path.join(HERE, n), encoding='utf-8').read()

sweep = []
for line in rd('partitions_flint_exceptions_new_37_1297.txt').splitlines():
    m = re.match(r'\s*ℓ =\s*(\d+): on lane ℓ the class has\s*(\d+) forms \(⌊\(ℓ\+12\)/24⌋ =\s*(\d+)\)(.*?)T_ℓ on the class: zero roots (\d+), repeated factors (\d+)', line)
    if not m: continue
    ex = re.search(r'exactly: (\d+) of the', m.group(4))
    sweep.append({'l': int(m.group(1)), 'lane': int(m.group(2)), 'R': int(m.group(3)), 'exact': int(ex.group(1)) if ex else int(m.group(2)),
                  'zero': int(m.group(5)), 'rep': int(m.group(6))})
assert len(sweep) == 200 and sweep[0]['l'] == 37 and sweep[-1]['l'] == 1297
flagged = [s['l'] for s in sweep if s['rep'] or s['zero']]
assert all(s['exact'] == s['R'] for s in sweep)

t = rd('partitions_flint_607_deeper.txt'); weights = []
for blk in re.finditer(r'weight (\d+) \(signs U₂ ([+-]\d), U₃ ([+-]\d)\): the class has (\d+) forms \(to 607\^\d\), (\d+) of them ordinary; '
                       r'T₆₀₇ ≡ 367 on a piece of dimension (\d+).*?\n.*?T₅ on the pair, on lane 607: \[([^\]]*)\].*?C = \(T₆₀₇ − 367\)/607: \[([^\]]*)\]\n'
                       r'\s*form A \(T₅ ≡ t\): δ = (\d+) \+ (\d+)·t; form B: δ = (\d+) \+ (\d+)·t', t):
    g = blk.groups()
    weights.append({'k': int(g[0]), 's2': int(g[1]), 's3': int(g[2]), 'forms': int(g[3]), 'ordinary': int(g[4]), 'piece': int(g[5]),
                    'T5': [int(x) for x in g[6].split(',')], 'C': [int(x) for x in g[7].split(',')],
                    'dA': [int(g[8]), int(g[9])], 'dB': [int(g[10]), int(g[11])]})
assert [w['k'] for w in weights] == [604, 1210, 1816]
sec = [int(x) for x in re.search(r'second difference .*? = (\d+) \+ (\d+)·t', t).groups()]
dm2 = [int(x) for x in re.search(r'δ\(−2\) = 2·δ\(604\) − δ\(1210\) = (\d+) \+ (\d+)·t', t).groups()]
tc = [int(x) for x in re.search(r'modulo 607²: 240 \+ 607·\((\d+) \+ (\d+)·t\) and 240 \+ 607·\((\d+) \+ (\d+)·t\)', t).groups()]
span = re.search(r'ℓ = 607 \(weight 604\).*?characteristic polynomial of degree (\d+), MINIMAL polynomial of degree (\d+)',
                 rd('partitions_flint_607_span_599_601_607_613_499.txt'))
OUT = {'sweep': [[s['l'], s['R'], s['rep']] for s in sweep], 'flagged': flagged, 'copy499': [s['lane'] for s in sweep if s['l'] == 499][0],
       'weights': weights, 'second_difference': sec, 'delta_m2': dm2, 'tower_const': [tc[:2], tc[2:]],
       'charpoly_degree': int(span.group(1)), 'minpoly_degree': int(span.group(2))}
print('sweep: %d primes, flagged %s; weights %s; second difference %s; δ(−2) %s; minimal polynomial degree %d of %d' % (
    len(sweep), flagged, [w['k'] for w in weights], sec, dm2, OUT['minpoly_degree'], OUT['charpoly_degree']))
if KEEP:
    json.dump(OUT, open(os.path.join(HERE, 'chords_page_data607c.json'), 'x'), separators=(',', ':'))
    print('wrote chords_page_data607c.json')
