# A dig, 13 Nov 2026: 607'S PAIR, DIRECTIONS FROM A FIXED REFERENCE (README "607 on a page"). chords_page_data607.py counted
# the pair's directions from form A's own place at 5 (its field was F_607[t]/(T₅'s polynomial), whose first generator is
# t = form A's number at 5), which puts the two forms at +1 and −1 at 5, next to each other only by that choice. Here the
# same numbers (read from the kept chords_page_data607.json) are carried into the standard model F_607(√3) (3 is not a
# square mod 607) by t ↦ τ = (−561 + s·√3)/2, s the smaller square root of (561² − 4·66)/3, and their directions are
# counted from E = G^606, G the first generator of F_(607²)^* in the order a + b·√3 (a, then b). Writes
# chords_page_data607b.json. DERIVED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/chords_page_data607b.py [--keep]
import os, sys, json, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'chords_page_data607b.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


src = json.load(open(os.path.join(HERE, 'chords_page_data607.json'), encoding='utf-8'))
L = src['L']; c0, c1 = src['chi5'][0], src['chi5'][1]                    # t² + c1·t + c0 = 0
assert pow(3, (L - 1) // 2, L) == L - 1                                    # 3 is not a square mod 607
mul = lambda x, y: ((x[0] * y[0] + 3 * x[1] * y[1]) % L, (x[0] * y[1] + x[1] * y[0]) % L)   # (a, b) = a + b·√3


def pw(x, n):
    r = (1, 0)
    while n:
        if n & 1: r = mul(r, x)
        x = mul(x, x); n >>= 1
    return r


disc = (c1 * c1 - 4 * c0) % L; q = disc * pow(3, -1, L) % L
s = min(v for v in range(L) if v * v % L == q)
half = pow(2, -1, L); tau = ((-c1) * half % L, s * half % L)
assert mul(tau, tau) == tuple(((-c1 * tau[0] - c0) % L, (-c1 * tau[1]) % L))   # τ² = −c1·τ − c0
emb = lambda ab: ((ab[0] + ab[1] * tau[0]) % L, (ab[1] * tau[1]) % L)          # a + b·t ↦ a + b·τ
NN = L * L - 1
G = next(g for g in itertools.product(range(L), repeat=2) if g != (0, 0) and all(pw(g, NN // p) != (1, 0) for p in (2, 3, 19, 101)))
E = pw(G, L - 1); dtab = {}; cur = (1, 0)
for i in range(L + 1): dtab[cur] = i; cur = mul(cur, E)
size = lambda u: pw(u, L + 1)[0]
direc = lambda u: dtab[pw(u, L - 1)]
say('F_607(√3); G = %d + %d·√3; τ = %d + %d·√3 (a root of t² + %d·t + %d)' % (G[0], G[1], tau[0], tau[1], c1, c0))
pair = {}
for key in ('at5', 'at607', 'deeper'):
    uA = emb(src['pair'][key]['A_el']); uB = emb(src['pair'][key]['B_el'])
    pair[key] = {'A': [size(uA), direc(uA)], 'B': [size(uB), direc(uB)]}
    assert pair[key]['A'][0] == src['pair'][key]['A'][0] and (pair[key]['A'][1] + pair[key]['B'][1]) % (L + 1) == 0
    say('   %s: form A size %d, direction %d of 608; form B size %d, direction %d of 608' % (key, *pair[key]['A'], *pair[key]['B']))
if KEEP: json.dump({'pair': pair, 'G': list(G), 'tau': list(tau)}, open(os.path.join(HERE, 'chords_page_data607b.json'), 'x'), separators=(',', ':'))
say('wrote chords_page_data607b.json' if KEEP else '(not kept)')
