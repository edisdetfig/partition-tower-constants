# A dig, 12 Nov 2026: THE DRIFT IS THE NEXT FROBENIUS ROOT (README "The drift is the next Frobenius root"; predicted there
# first). At ℓ = 5, poles of 1/(E₄³ − j₀Δ) (weight −12) with j₀ ≡ 2 (mod 5), s = v₅(j₀ + 32768) digits from the twin j = −32768.
# Laws P(x)·∏_(i<I)(x − 5ⁱ·u^(−14−2i)), I = 1, 2, 3, u the unit root of x² + 3x + 5. By level: the smallest depth over n ≥ 1.
# COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/drift_roots_check.py [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'drift_roots_check.txt'), 'x', encoding='utf-8') if KEEP else None


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


L = 5; JL = -32768
mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(os.path.join(HERE, 'forms1_family_w-12_5.txt'), encoding='utf-8').read())
P0 = [int(v) for v in mt.group(2).split(',')]; KP = int(mt.group(1))
a = -3; M0 = L ** 40; u = a % L
for _ in range(200): u = (u - (u * u - a * u + L) * pow((2 * u - a) % M0, -1, M0)) % M0
say('ℓ = 5, weight −12: laws P(x)·∏_(i<I) (x − 5ⁱ·u^(−14−2i)); depth by starting level (lane 5^min(K))')
for j0, f in ((JL, 'level1_tower_poleJ0-32768_N16000000_5.json'), (2, 'level1_tower_poleJ02_N16000000_5_7_11_13_17_19_23_29_31_37_41_43.json'),
              (-32743, 'level1_tower_poleJ0-32743_N16000000_5.json'), (-32643, 'level1_tower_poleJ0-32643_N16000000_5.json'),
              (-32143, 'level1_tower_poleJ0-32143_N16000000_5.json'), (-17143, 'level1_tower_poleJ0-17143_N16000000_5.json')):
    rec = json.load(open(os.path.join(HERE, f), encoding='utf-8'))['primes']['5']; K = min(rec['K'], KP); M = L ** K
    lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}; top = max(lv)
    s = val(j0 - JL, L, 60); rows = []
    for I in (1, 2, 3):
        Q = [c % M for c in P0]
        for i in range(I): Q = pmul(Q, [(-(L ** i) * pow(u, -14 - 2 * i, M)) % M, 1], M)
        D = len(Q) - 1; ds = []
        for b in sorted(lv):
            if b < 1 or b + D > top: continue
            nn = min(len(lv[b + t]) for t in range(D + 1))
            dd = [val(sum(Q[t] * lv[b + t][n] for t in range(D + 1)), L, K) for n in range(1, nn)]
            if dd: ds.append(min(dd))
        rows.append('I = %d: %s' % (I, ds))
    say('   j₀ = %7d (%s):  %s' % (j0, 'the twin' if j0 == JL else 's = %d' % s, ';   '.join(rows)))
