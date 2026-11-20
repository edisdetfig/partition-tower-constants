# A dig, 12 Nov 2026: THE MOCK THETA FUNCTION f(q), THE TWO-LEVEL CONSTANT TO EVERY DIGIT AT 5 AND 7 (README "A mock theta
# function"). PREDICTED before the first run: the constant is the shadow's T(ℓ²)-eigenvalue χ₁₂(ℓ)·(ℓ + (−1|ℓ)) (−6 at 5 and 7);
# WRONG at 5 (the constant is −1). Corrected, before 7 was read: the unit root (−3|ℓ), and the shadow's whole polynomial.
# Read here from the tower alone:
#   - the ratio H_(b+2)(n)/H_b(n) over every n with H_b(n) a unit, level by level, and the digits all n agree on;
#   - then the law x² − λ, for each candidate λ, by level: the smallest depth over n ≥ 1.
# COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/mock_theta_digits.py <mock_theta_tower json> [--keep]
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRC = next(a for a in sys.argv[1:] if a.endswith('.json'))
fh = open(os.path.join(HERE, 'mock_theta_digits.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def digits(x, L, k): return [(x // L ** e) % L for e in range(k)]


data = json.load(open(os.path.join(HERE, SRC), encoding='utf-8'))
for L in (5, 7):
    rec = data['primes'][str(L)]; K = rec['K']; M = L ** K
    lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}; top = max(lv)
    c12, m1, m3 = leg(3, L), leg(-1, L), leg(-3, L)
    lam = c12 * (L + m1)
    cands = {'y − χ₁₂(ℓ)(ℓ + (−1|ℓ)) (the first prediction)': [(-lam), 1], 'y − (−3|ℓ)(1 + ℓ)': [-(m3 * (1 + L)), 1],
             'y − (−3|ℓ)(1 − ℓ)': [-(m3 * (1 - L)), 1], 'y − (−3|ℓ), the unit root (corrected)': [-m3, 1],
             'y² − χ₁₂(ℓ)(ℓ + (−1|ℓ))y + (−1|ℓ)ℓ, the whole polynomial of the shadow': [m1 * L, -lam, 1]}
    say(); say('ℓ = %d, lane ℓ^%d, levels 0..%d: χ₁₂(ℓ) = %+d, (−1|ℓ) = %+d, (−3|ℓ) = %+d' % (L, K, top, c12, m1, m3))
    for b in sorted(lv):
        if b < 1 or b + 2 > top: continue
        rs = []
        for n in range(1, min(len(lv[b]), len(lv[b + 2]))):
            if lv[b][n] % L == 0: continue
            rs.append(lv[b + 2][n] * pow(lv[b][n], -1, M) % M)
        if not rs: continue
        agree = next((e for e in range(K) if len({r % L ** (e + 1) for r in rs}) > 1), K)
        say('   levels %d → %d: the ratio over %d values of n agrees to %d digits: %s' % (b, b + 2, len(rs), agree, digits(rs[0], L, agree)))
    for name, Q in cands.items():
        Q = [c % M for c in Q]; D = len(Q) - 1; rows = []
        for b in sorted(lv):
            if b < 1 or b + 2 * D > top: continue
            nn = min(len(lv[b + 2 * i]) for i in range(D + 1))
            if nn <= 1: continue
            rows.append('%d: %d' % (b, min(val(sum(Q[i] * lv[b + 2 * i][n] for i in range(D + 1)), L, K) for n in range(1, nn))))
        say('   %-62s %s: depth by level  %s' % (name, [c if c <= M // 2 else c - M for c in Q], '  '.join(rows)))
