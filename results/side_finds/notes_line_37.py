# A dig, 13 Nov 2026: 37'S LINES THROUGH THE LEVELS, READ CORRECTLY (README "Notes in a line"). notes_line_data.py read the
# kept 37 data (chords_page_data.json, "A") as values; they are discrete logarithms to the root 2 (None for 0), so its 37
# line ("the two-step rule holds 1 of 38 times") was a misreading, not a result. Here the logarithms are turned back into
# values mod 37, and the same checks are made: one multiplier per n through the levels; the two-step rule
# H_(b+4) ≡ 12·H_(b+2) − 29·H_b (mod 37); each note's part times its own constant (15 and 34) per two levels. Writes
# notes_line_37.json. COUNTED (from the kept values).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/notes_line_37.py [--keep]
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'notes_line_37.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


A = json.load(open(os.path.join(HERE, 'chords_page_data.json'), encoding='utf-8'))['A']
L, root, lv = A['ell'], A['root'], A['levels']; T, N = 12, 29
val = lambda t: 0 if t is None else pow(root, t, L)
LV = [1, 3, 5, 7, 9]; ns = lv['1']['n']
H = {b: [val(t) for t in lv[str(b)]['chord']] for b in LV}
P = {key: {b: [val(t) for t in lv[str(b)][key]] for b in LV} for key in ('A', 'B')}
cst = {'A': A['alpha']['A'][0], 'B': A['alpha']['B'][0]}
say('ℓ = 37, root %d; the notes %s; %d values of n: %s' % (root, cst, len(ns), ns))
one = rule_ok = rule_n = 0
for i, n in enumerate(ns):
    vv = [H[b][i] for b in LV if i < len(H[b])]
    if len(vv) == 5 and all(vv):
        one += len({vv[j + 1] * pow(vv[j], -1, L) % L for j in range(4)}) == 1
    for j in range(len(vv) - 2):
        rule_n += 1; rule_ok += (T * vv[j + 1] - N * vv[j] - vv[j + 2]) % L == 0
parts = {}
for key in ('A', 'B'):
    ok = tot = 0
    for i in range(len(ns)):
        v = [P[key][b][i] for b in LV if i < len(P[key][b])]
        for j in range(len(v) - 1):
            if v[j] == 0 and v[j + 1] == 0: continue
            tot += 1; ok += (v[j + 1] - cst[key] * v[j]) % L == 0
    parts[key] = [ok, tot]
say('one multiplier fits all five levels for %d of %d values of n' % (one, len(ns)))
say('the two-step rule H(b+4) ≡ %d·H(b+2) − %d·H(b) holds %d of %d times' % (T, N, rule_ok, rule_n))
say("each note's part times its own constant per two levels: note %d %d of %d, note %d %d of %d" % (
    cst['A'], parts['A'][0], parts['A'][1], cst['B'], parts['B'][0], parts['B'][1]))
OUT = {'n': ns, 'H': {str(b): H[b] for b in LV}, 'A': {str(b): P['A'][b] for b in LV}, 'B': {str(b): P['B'][b] for b in LV},
       'T': T, 'N': N, 'single_line': one, 'rule': [rule_n, rule_ok], 'parts': parts, 'const': cst}
if KEEP: json.dump(OUT, open(os.path.join(HERE, 'notes_line_37.json'), 'x'), separators=(',', ':')); say('wrote notes_line_37.json')
