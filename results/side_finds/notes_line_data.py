# A dig, 13 Nov 2026: NOTES IN A LINE (README "Notes in a line: patterns across rows, and one line through the levels";
# predicted there first). Written for the page "Notes in a Line":
#   1. the rows: every ordinary note on lane ℓ (13 to 97 from chords_page_data2/4.json, 607 from chords_page_data607.json)
#      walks 1, α, α², … back to 1; the braids' sizes (braid_polar_data.json, chords_page_data4.json) as extra rows;
#   2. the shift scan: for rows-apart Δr = 1 … 6 and steps-over Δc = −12 … 12, the circular correlation of the rows' places
#      on the clock (cos 2π(x − y), averaged over 60 steps and over the row pairs), against 2,000 shuffles of each row's
#      steps; z-scores, and the largest |z| against the largest |z| of each shuffle (the multiple-testing yardstick);
#   3. the stretches within a prime: row b at step j against row a at step e·j, for every e (perfect where b = a^e);
#   4. the lines through the levels: H_b(n) mod ℓ from exact partition numbers (FLINT) at 13 … 31, levels 1 … 9,
#      n = 0 … 15 (N up to 4.5·10¹⁴), checked against log H₁(n) + j·log α; at 37 the kept values (chords_page_data.json),
#      checked for one multiplier per n and for the two-step rule with (T, N) = (12, 29).
# Writes notes_line_data.json. COUNTED (the rows, the tower values, the checks), DERIVED (the expectations).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/notes_line_data.py [--keep]
import os, sys, json, math, time
HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import numpy as np
import flint
from sympy import primefactors
flint.ctx.threads = 4
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'notes_line_data.txt'), 'x', encoding='utf-8') if KEEP else None
rd = lambda n: json.load(open(os.path.join(HERE, n), encoding='utf-8'))


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def prim_root(L): return next(a for a in range(2, L) if all(pow(a, (L - 1) // q, L) != 1 for q in primefactors(L - 1)))


def walk(a, L):
    w = [1]; x = a % L
    while x != 1: w.append(x); x = x * a % L
    return w + [1]                                                    # from 1 back to 1


t0 = time.time()
# ---- 1. the rows ----
S = dict(rd('chords_page_data2.json')['S']); S.update(rd('chords_page_data4.json')['S'])
notes = []
for L in sorted(S, key=int):
    for a in sorted({nt['first'] for nt in S[L]['notes'] if nt['e'] == 1}): notes.append((int(L), a))
d607 = rd('chords_page_data607.json')
for a in sorted({nt['first'] for nt in d607['notes'] if nt['e'] == 1}): notes.append((607, a))
rows = []
for L, a in notes:
    g = prim_root(L); lg = {}; x = 1
    for i in range(L - 1): lg[x] = i; x = x * g % L
    w = walk(a, L); rows.append({'l': L, 'a': a, 'K': len(w) - 1, 'g': g, 'theta': lg[a], 'walk': w})
bp = rd('braid_polar_data.json'); bp.update(rd('chords_page_data4.json')['polar'])
brows = []
for L in sorted(bp, key=int):
    Na = bp[L]['Na']; w = walk(Na, int(L)); brows.append({'l': int(L), 'a': Na, 'K': len(w) - 1, 'e': bp[L]['e'], 'walk': w})
say('rows: %d ordinary notes (%s), %d braid sizes' % (len(rows), ', '.join('%d:%d' % (r['l'], r['a']) for r in rows), len(brows)))

# ---- 2. the shift scan ----
J, PAD, DRS, DCS = 60, 12, range(1, 7), range(-12, 13)
def zmat(rs, rng=None):
    out = []
    for r in rs:
        base = r['walk'][:-1]
        if rng is not None: base = list(rng.permutation(base))
        seq = [base[(j - PAD) % len(base)] for j in range(J + 2 * PAD)]
        out.append(np.exp(2j * np.pi * np.array(seq, dtype=float) / r['l']))
    return np.array(out)
def scan(Z):
    C = np.zeros((len(DRS), len(DCS)))
    for i, dr in enumerate(DRS):
        A = Z[:-dr, PAD:PAD + J]
        for k, dc in enumerate(DCS):
            B = Z[dr:, PAD + dc:PAD + dc + J]
            C[i, k] = np.mean(np.real(A * np.conj(B)))
    return C
C0 = scan(zmat(rows)); rng = np.random.default_rng(20261113)
null = np.array([scan(zmat(rows, rng)) for _ in range(2000)])
mu, sd = null.mean(axis=0), null.std(axis=0)
Zs = (C0 - mu) / sd; maxz = float(np.abs(Zs).max()); arg = np.unravel_index(np.abs(Zs).argmax(), Zs.shape)
nullmax = np.abs((null - mu) / sd).reshape(2000, -1).max(axis=1); p95 = float(np.percentile(nullmax, 95))
pmax = float((nullmax >= maxz).mean())
say('shift scan: %d offsets (rows apart 1–6, steps over −12…12); the largest |z| is %.2f at rows apart %d, steps over %+d; '
    'the shuffles\' largest |z| is above that %.0f%% of the time (95th percentile %.2f)' % (Zs.size, maxz, DRS[arg[0]], DCS[arg[1]], 100 * pmax, p95))

# ---- 3. the stretches within a prime ----
stretch = []
for i, ra in enumerate(rows):
    for rb in rows:
        if rb is ra or rb['l'] != ra['l']: continue
        wa, wb, Ka = ra['walk'][:-1], rb['walk'][:-1], ra['K']
        for e in range(1, Ka):
            if all(wb[j % rb['K']] == wa[(e * j) % Ka] for j in range(max(Ka, rb['K']))):
                stretch.append({'l': ra['l'], 'a': ra['a'], 'b': rb['a'], 'e': e}); break
say('perfect stretches (row b at step j = row a at step e·j): %s' % ['%d: %d = %d^%d' % (s['l'], s['b'], s['a'], s['e']) for s in stretch])

# ---- 4. the lines through the levels ----
PC = {}
def c(M):
    N = (M + 1) // 24
    if N not in PC: PC[N] = int(flint.fmpz(N).partitions_p())
    return PC[N]
lines = {}
for L in (13, 17, 19, 23, 29, 31):
    R = next(r for r in rows if r['l'] == L); g, th = R['g'], R['theta']
    lg = {}; x = 1
    for i in range(L - 1): lg[x] = i; x = x * g % L
    H = {}
    for b in (1, 3, 5, 7, 9):
        kk = (-pow(L ** b, -1, 24)) % 24; vals = []
        for n in range(16):
            M = L ** b * (24 * n + kk)
            if (M + 1) // 24 > 4.5e14: break
            vals.append(c(M) % L)
        H[b] = vals
    on, off, zero = 0, 0, 0
    for n in range(16):
        h1 = H[1][n]
        for j, b in enumerate((1, 3, 5, 7, 9)):
            if n >= len(H[b]): continue
            if h1 == 0: zero += (H[b][n] == 0); off += (H[b][n] != 0); continue
            if H[b][n] != 0 and lg[H[b][n]] == (lg[h1] + j * th) % (L - 1): on += 1
            else: off += 1
    lines[L] = {'alpha': R['a'], 'g': g, 'theta': th, 'H': {str(b): v for b, v in H.items()}, 'on': on, 'off': off, 'zero': zero}
    say('ℓ = %d: constant %d = %d^%d; tower values at levels 1–9 for n = 0…15 (%s values per level): %d on the line from the start, '
        '%d off it, %d zeros that stay 0  (%.0f s)' % (L, R['a'], g, th, [len(H[b]) for b in (1, 3, 5, 7, 9)], on, off, zero, time.time() - t0))
A = rd('chords_page_data.json')['A']; lv = A['levels']; L = 37; T, Nn = 12, 29
ns = lv['1']['n']; ch = {b: lv[str(b)]['chord'] for b in (1, 3, 5, 7, 9)}
one = 0; rule = [0, 0]
for i, n in enumerate(ns):
    vs = [ch[b][i] if i < len(ch[b]) else None for b in (1, 3, 5, 7, 9)]
    vv = [v for v in vs if v is not None]
    if len(vv) >= 3 and all(v != 0 for v in vv):
        ratios = {vv[j + 1] * pow(vv[j], -1, L) % L for j in range(len(vv) - 1)}
        one += (len(ratios) == 1)
    for j in range(len(vv) - 2):
        rule[0] += 1; rule[1] += ((T * vv[j + 1] - Nn * vv[j] - vv[j + 2]) % L == 0)
lines['37'] = {'n': ns, 'H': {str(b): ch[b] for b in ch}, 'T': T, 'N': Nn, 'single_line': one, 'rule': rule,
               'A': {str(b): lv[str(b)]['A'] for b in ch}, 'B': {str(b): lv[str(b)]['B'] for b in ch}}
say('ℓ = 37: %d values of n; one multiplier fits all levels for %d of them; the two-step rule (12, 29) holds %d of %d times' % (len(ns), one, rule[1], rule[0]))
OUT = {'rows': rows, 'braids': brows, 'scan': {'C': C0.round(4).tolist(), 'z': Zs.round(2).tolist(), 'dr': list(DRS), 'dc': list(DCS),
       'maxz': round(maxz, 2), 'at': [DRS[arg[0]], DCS[arg[1]]], 'p95': round(p95, 2), 'pmax': round(pmax, 3), 'J': J},
       'stretch': stretch, 'lines': {str(k): v for k, v in lines.items()}, 'maxN': max(PC) if PC else 0}
if KEEP: json.dump(OUT, open(os.path.join(HERE, 'notes_line_data.json'), 'x'), separators=(',', ':'))
say('the largest partition number used: p(%d); seconds: %.0f' % (max(PC), time.time() - t0))
