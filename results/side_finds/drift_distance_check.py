# A dig, 12 Nov 2026: THE DISTANCE FROM A CURVE TO ITS RETURNING TWIN, READ FROM A TOWER (README "What the tower view
# reveals"). At ℓ = 5, the poles of 1/(E₄³ − j₀Δ) with j₀ ≡ 2 (mod 5) all drift toward the canonical lift j = −32768 (CM by
# √−11). A curve s = v₅(j₀ + 32768) digits from it is s + b digits from it at level b (Serre–Tate; derived).
# PREDICTED before the runs: the law P(x)·(x − u⁻¹⁴) holds to depth min(2b + 3, b + s, lane) from level b.
#   - 2b + 3 is the pole ON the lift (measured, repeat_check_A.txt);
#   - b + s is the drift, with the constant fixed by j = 2 (s = 1, measured b + 1).
# COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/drift_distance_check.py [--keep]
import os, sys, re, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'drift_distance_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


L = 5; JL = -32768
mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(os.path.join(HERE, 'forms1_family_w-12_5.txt'), encoding='utf-8').read())
P = [int(v) for v in mt.group(2).split(',')]; KP = int(mt.group(1))
a = -3; M0 = L ** 40; u = a % L                                     # the curve j ≡ 2 mod 5 has trace −3 (counted earlier)
for _ in range(200): u = (u - (u * u - a * u + L) * pow((2 * u - a) % M0, -1, M0)) % M0
say('ℓ = 5: poles of 1/(E₄³ − j₀Δ) with j₀ ≡ 2 (mod 5), drifting toward the lift j = %d; law P(x)·(x − u⁻¹⁴)' % JL)
runs = [(JL, 'level1_tower_poleJ0-32768_N16000000_5.json')]
for j0 in (2, -32743, -32643, -32143, -17143):
    f = 'level1_tower_poleJ0%d_N16000000_%s.json' % (j0, '5_7_11_13_17_19_23_29_31_37_41_43' if j0 == 2 else '5')
    runs.append((j0, f))
for j0, f in runs:
    d = json.load(open(os.path.join(HERE, f), encoding='utf-8')); rec = d['primes']['5']
    K = min(rec['K'], KP); M = L ** K
    c14 = pow(u, -14, M); Q = [0] * (len(P) + 1)
    for i, c in enumerate(P): Q[i + 1] = (Q[i + 1] + c) % M; Q[i] = (Q[i] - c14 * c) % M
    lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}; top = max(lv); D = len(Q) - 1
    meas = []
    for b in sorted(lv):
        if b < 1 or b + D > top: continue
        nn = min(len(lv[b + i]) for i in range(D + 1))
        ds = [val(sum(Q[i] * lv[b + i][n] for i in range(D + 1)), L, K) for n in range(1, nn)]
        if ds: meas.append((b, min(ds), len(ds)))
    s = val(j0 - JL, L, 60)
    pred = [(b, min(2 * b + 3, b + s, K)) if j0 != JL else (b, min(2 * b + 3, K)) for b, _, _ in meas]
    agree = sum(1 for (b, m, _), (_, p) in zip(meas, pred) if m == p)
    say('   j₀ = %7d (%s): measured %s;  predicted %s;  %d of %d levels as predicted (%d checks)' % (
        j0, 'ON the lift' if j0 == JL else 's = %d' % s, [m for _, m, _ in meas], [p for _, p in pred], agree, len(meas),
        sum(c for _, _, c in meas)))
