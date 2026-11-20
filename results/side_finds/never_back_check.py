# A dig, 12 Nov 2026: CURVES THAT NEVER COME BACK (README "Curves that never come back"; predictions written there first).
# At ℓ = 5, weight −48, three poles reducing to j ≡ 2 (mod 5), whose returning curve is j = −32768 (CM by √−11):
#   (a) D275: 1/(Δ⁴·H₋₂₇₅(j)), the four quasi-canonical lifts of level 1, which land on the returning curve after one step;
#   (b) NC4:  1/(Δ⁴·(j − 2)(j − 7)(j − 12)(j − 17)), four curves without CM, which approach it forever;
#   (c) LIFT4: 1/(E₄⁹·(E₄³ + 32768Δ)), a simple pole on the returning curve (and one at ρ, supersingular at 5).
# Law P(x)·(x − u⁻⁵⁰)·S(x): u the unit root of x² + 3x + 5 (trace −3 at 5), P the weight −48 ordinary families
# (forms1_family_w-48_5.txt), S the slope-1 family (slope1_families_5_w-48_N20.txt). By level: the smallest depth over n ≥ 1.
# COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/never_back_check.py [--keep]
import os, sys, re, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'never_back_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def pmul(A, B, M):
    out = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        for j, b in enumerate(B): out[i + j] = (out[i + j] + a * b) % M
    return out


L = 5; W = -48
mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(os.path.join(HERE, 'forms1_family_w%d_%d.txt' % (W, L)), encoding='utf-8').read())
P = [int(v) for v in mt.group(2).split(',')]; KP = int(mt.group(1))
ms = re.search(r'the slope-1 factor ∏\(y − λ/ℓ\) modulo ℓ\^(\d+) \(truncated at t = \d+\), lowest first: \[([^\]]*)\]',
               open(os.path.join(HERE, 'slope1_families_%d_w%d_N20.txt' % (L, W)), encoding='utf-8').read())
G1 = [int(v) for v in ms.group(2).split(',')]; h1 = len(G1) - 1; KS = int(ms.group(1)) + 1
a = -3; M0 = L ** 40; u = a % L
for _ in range(200): u = (u - (u * u - a * u + L) * pow((2 * u - a) % M0, -1, M0)) % M0
say('ℓ = 5, weight %d: P on lane ℓ %s; %d slope-1 family; law P(x)·(x − u^%d)·S(x)' % (W, [c % L for c in P], h1, W - 2))
for tagn, what in (('LIFT4', '(c) a pole ON the returning curve j = −32768'),
                   ('D275', '(a) the four curves of discriminant −275 (fall onto it after one step)'),
                   ('NC4', '(b) four curves without CM (approach it forever)')):
    f = glob.glob(os.path.join(HERE, 'level1_tower_poleH%s_N16000000_5.json' % tagn))
    if not f: say('   %s: tower missing' % tagn); continue
    rec = json.load(open(f[0], encoding='utf-8'))['primes']['5']
    K = min(rec['K'], KP, KS); M = L ** K
    lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}; top = max(lv)
    c = pow(u, W - 2, M)
    S = [G1[s] * L ** (h1 - s) % M for s in range(h1 + 1)]
    laws = [('P alone', [x % M for x in P]), ('P·(x − u^%d)' % (W - 2), pmul([x % M for x in P], [(-c) % M, 1], M)),
            ('P·(x − u^%d)·S' % (W - 2), pmul(pmul([x % M for x in P], [(-c) % M, 1], M), S, M))]
    say(); say('%s (lane 5^%d, levels 0..%d):' % (what, K, top))
    for name, Q in laws:
        D = len(Q) - 1; rows = []; checks = 0
        for b in sorted(lv):
            if b < 1 or b + D > top: continue
            nn = min(len(lv[b + i]) for i in range(D + 1))
            ds = [val(sum(Q[i] * lv[b + i][n] for i in range(D + 1)), L, K) for n in range(1, nn)]
            if ds: rows.append('%d: %d' % (b, min(ds))); checks += len(ds)
        say('   %-22s depth by level  %s   (%d checks)' % (name, '  '.join(rows), checks))
