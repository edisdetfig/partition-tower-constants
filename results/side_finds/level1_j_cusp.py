# A dig, 12 Nov 2026: j's TOWER FOR n ≥ 1 FOLLOWS THE CUSP FORMS ALONE. At weight 0 the Eisenstein family's member is the
# constant 1 (E_k* = 1 + (2/ζ*(1−k))·Σσ*_(k−1)(n)qⁿ, and ζ*(1−k) has its pole at k = 0), so the Eisenstein family carries
# only j's constant term 744. PREDICTED before the run (derived): for n ≥ 1, Σ_i [xⁱ]C·H_(b+i)(n) ≡ 0, deeper level by
# level, with C(x) = P(x)/(x − 1) the cusp forms' polynomial at weight 0 (forms1_family_w0_ℓ.txt); at 5, 7 and 11, where C
# = 1, this says c(ℓᵇn) → 0: Lehner's congruences c(5ᵃn) ≡ 0 (mod 5^(a+1)), c(7ᵃn) ≡ 0 (mod 7ᵃ), c(11ᵃn) ≡ 0 (mod 11ᵃ).
# Also printed: the ratio c(ℓ^(b+1)n)/c(ℓᵇn) on lane ℓ where C has degree 1 (one cusp form: 13, 17, 19, 23). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/level1_j_cusp.py <j tower json> [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRC = next(a for a in sys.argv[1:] if a.endswith('.json'))
fh = open(os.path.join(HERE, 'level1_j_cusp_%s.txt' % os.path.basename(SRC)[13:-5]), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


data = json.load(open(os.path.join(HERE, SRC), encoding='utf-8')); assert data['form'] == 'j'
say('j towers: %s (n < %d); n = 1…11 at each level' % (SRC, data['N']))
for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
    L = int(Ls); f = os.path.join(HERE, 'forms1_family_w0_%d.txt' % L)
    mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(f, encoding='utf-8').read())
    KF = int(mt.group(1)); P = [int(v) for v in mt.group(2).split(',')]; D = len(P) - 1
    K = min(KF, rec['K']); M = L ** K
    C = [0] * D; acc = 0
    for i in range(D, 0, -1): acc = (P[i] + acc) % M; C[i - 1] = acc        # P / (x − 1), highest first
    lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}; top = max(lv)
    rows = []
    for b in sorted(lv):
        if b + D - 1 > top: continue
        nn = min(len(lv[b + i]) for i in range(D))
        ds = [val(sum(C[i] * lv[b + i][n] for i in range(D)), L, K) for n in range(1, nn)]
        if ds: rows.append('%d: %d' % (b, min(ds)))
    say('ℓ = %d, lane ℓ^%d: the cusp forms\' polynomial C on lane ℓ %s (%d constant%s); the smallest depth over n ≥ 1, '
        'from level b:  %s' % (L, K, [c % L for c in C], D - 1, '' if D == 2 else 's', ';  '.join(rows) or 'no level span'))
    if D == 2:
        a = (-C[0]) % M
        rs = set()
        for b in sorted(lv):
            if b >= 1 and b + 1 in lv:
                for n in range(1, min(len(lv[b]), len(lv[b + 1]))):
                    if lv[b][n] % L: rs.add(lv[b + 1][n] * pow(lv[b][n], -1, L) % L)
        say('   c(ℓ^(b+1)n)/c(ℓᵇn) on lane ℓ for b ≥ 1, n ≥ 1: %s;  the cusp form\'s constant, digits lowest first: %s' % (
            sorted(rs), [(a // L ** e) % L for e in range(K)]))
