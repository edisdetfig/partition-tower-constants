# A dig, 12 Nov 2026: 1/E₄², THE WHOLE LAW FIXED IN ADVANCE. At ℓ ≡ 1 mod 3 the tower of 1/E₄² (weight −8, a double pole at
# ρ) should obey
#     P(x)·(x − ℓ·u⁻¹²)·S(x)
#   - P: the ordinary families at weight −8 (forms1_family_w-8_ℓ.txt; only the Eisenstein 1 at 7, 13, 19);
#   - ℓ·u⁻¹²: the double pole's constant (derived, README "Double poles"; u the unit root of y² = x³ + 1 at ℓ);
#   - S(x) = ∏(x − ℓ·y): the slope-1 families at weight −8, from the forms (slope1_families_ℓ_w-8_N….txt).
# No constant is read from the tower. PREDICTED before the run (README "What 13 needs"): the law holds to a depth growing with
# the level; dropping the pole's factor or S makes it fail; at 7 it decides between the forms' family constant (unit part
# 4, 0, 1, 6, …) and the earlier data fit's (4, 0, 1, 3). Per law and starting level b ≥ 1: the smallest depth over n ≥ 1.
# COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/double_pole_full_law.py <invE4sq tower json> [--keep]
import os, sys, re, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRC = next(a for a in sys.argv[1:] if a.endswith('.json'))
fh = open(os.path.join(HERE, 'double_pole_full_law_%s.txt' % os.path.basename(SRC)[13:-5]), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def unit_root(L, K):
    M = L ** K
    a = next(a for a in range(1, 2 * L) if 4 * L - a * a > 0 and (4 * L - a * a) % 3 == 0
             and round(((4 * L - a * a) // 3) ** 0.5) ** 2 == (4 * L - a * a) // 3)
    u = a % L
    for _ in range(4 * K + 8): u = (u - (u * u - a * u + L) * pow((2 * u - a) % M, -1, M)) % M
    return u


def pmul(A, B, M):
    out = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        for j, b in enumerate(B): out[i + j] = (out[i + j] + a * b) % M
    return out


data = json.load(open(os.path.join(HERE, SRC), encoding='utf-8')); assert data['form'] == 'invE4sq'
say('1/E₄² towers: %s (n < %d); every law\'s constants fixed by the forms and the pole, none read from the tower' % (SRC, data['N']))
for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
    L = int(Ls)
    if L % 3 != 1: continue
    mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(os.path.join(HERE, 'forms1_family_w-8_%d.txt' % L), encoding='utf-8').read())
    srcs = sorted(glob.glob(os.path.join(HERE, 'slope1_families_%d_w-8_N*.txt' % L)))
    if not srcs: say('ℓ = %d: no slope-1 families computed' % L); continue
    ms = re.search(r'the slope-1 factor ∏\(y − λ/ℓ\) modulo ℓ\^(\d+) \(truncated at t = \d+\), lowest first: \[([^\]]*)\]',
                   open(srcs[-1], encoding='utf-8').read())
    M1 = int(ms.group(1)); G1 = [int(v) for v in ms.group(2).split(',')]; h1 = len(G1) - 1
    K = min(rec['K'], int(mt.group(1)), M1 + 1); M = L ** K
    P = [int(v) % M for v in mt.group(1 + 1).split(',')]
    S = [G1[s] * L ** (h1 - s) % M for s in range(h1 + 1)]          # ∏(x − ℓy) = Σ_s G1[s]·ℓ^(h1 − s)·x^s
    pole = [(-L * pow(unit_root(L, K), -12, M)) % M, 1]
    lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}; top = max(lv)
    say(); say('ℓ = %d, lane ℓ^%d: P on lane ℓ %s; %d slope-1 families (from %s, to ℓ^%d)' % (
        L, K, [c % L for c in P], h1, os.path.basename(srcs[-1]), M1))
    laws = [('P alone', P), ('P·(x − ℓu⁻¹²): the pole only', pmul(P, pole, M)), ('P·S: the families only', pmul(P, S, M)),
            ('P·(x − ℓu⁻¹²)·S: PREDICTED', pmul(pmul(P, pole, M), S, M))]
    if L == 7:
        yfit = 4 + 0 * 7 + 1 * 49 + 3 * 343                         # the earlier two-term fit's family (README "Double poles")
        laws.append(('P·(x − ℓu⁻¹²)·(x − 7·y_fit), y_fit = 4, 0, 1, 3 (the earlier fit)', pmul(pmul(P, pole, M), [(-7 * yfit) % M, 1], M)))
    for name, Q in laws:
        D = len(Q) - 1; rows = []; checks = 0
        for b in sorted(lv):
            if b < 1 or b + D > top: continue
            nn = min(len(lv[b + i]) for i in range(D + 1))
            ds = [val(sum(Q[i] * lv[b + i][n] for i in range(D + 1)), L, K) for n in range(1, nn)]
            if ds: rows.append('%d: %d' % (b, min(ds))); checks += len(ds)
        say('   %s (degree %d): smallest depth over n ≥ 1 from level  %s   (%d checks)' % (name, D, ';  '.join(rows) or 'no level span', checks))
