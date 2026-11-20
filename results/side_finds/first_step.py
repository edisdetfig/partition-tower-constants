# first_step.py, 14 Nov 2026 (README "Why exactly 13, 17, 19 and 23 settle fast").
# Exact dimensions c(w; ε₂, ε₃) of the newforms of level 6 and weight w with U_2 = ε₂·2^((w−2)/2) and
# U_3 = ε₃·3^((w−2)/2), for every even w = 4 … WMAX and all four sign pairs. Each is the nullity of the stacked matrix
# [U_2 − c₂; U_3 − c₃] on M_w(Γ0(6)) (forms6_flint's basis); no old form and no Eisenstein series has these numbers
# exactly (Deligne's bound). The nullity is taken modulo two large primes, which must agree (it equals the dimension over
# Q unless the prime divides a minor). Then README's predictions:
#   (a) c(w + 24; ε) = c(w; ε) + 1;
#   (b) the class's first-step gain c(2ℓ − 4; +, +) − c(ℓ − 3; (2|ℓ), (3|ℓ)) equals c(ℓ − 5; (2|ℓ), (3|ℓ)) + c(6; +, +);
#   (c) both equal ⌊ℓ/24⌋;
#   (d) c(6; +, +) = 0;
# and the residue classes mod 24 that (a) carries to every prime.
# COUNTED (the dimensions), DERIVED (the comparisons).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/first_step.py [WMAX] [--keep]
import os, sys, time, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat, fmpz
import forms6_flint as F

flint.ctx.threads = 2
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; WMAX = nums[0] if nums else 610; KEEP = '--keep' in sys.argv
base = os.path.join(HERE, 'first_step_w%d' % WMAX); OUT, JS = base + '.txt', base + '.jsonl'
if KEEP and (os.path.exists(OUT) or os.path.exists(JS)): sys.exit('stopped before computing anything: %s exists' % OUT)
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None
jh = open(JS, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def prev_prime(n):
    while not fmpz(n).is_prime(): n -= 1
    return n


PRIMES = [2 ** 61 - 1, prev_prime(2 ** 60)]
SIGNS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
sg = lambda s: '+' if s > 0 else '−'


def dims_at(w, P):
    n = w + 1
    H = F.hecke_matrices(w, P, [('U', 2), ('U', 3)])
    out = {}
    for s2, s3 in SIGNS:
        c2 = s2 * pow(2, (w - 2) // 2, P) % P; c3 = s3 * pow(3, (w - 2) // 2, P) % P
        A2 = H[('U', 2)] - F.identity(n, P) * c2; A3 = H[('U', 3)] - F.identity(n, P) * c3
        S = nmod_mat(2 * n, n, [int(v) for v in A2.entries()] + [int(v) for v in A3.entries()], P)
        R, rk = S.rref(); out[(s2, s3)] = n - rk
    return out


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def is_prime(n): return n > 1 and fmpz(n).is_prime()


say('the newforms of level 6 by their signs at 2 and 3: c(w; ε₂, ε₃) for even w = 4 … %d, as nullities modulo the '
    'primes %s (they must agree)' % (WMAX, PRIMES))
c = {}; t0 = time.time(); disagree = []
for w in range(4, WMAX + 1, 2):
    d1 = dims_at(w, PRIMES[0]); d2 = dims_at(w, PRIMES[1])
    if d1 != d2: disagree.append(w)
    c[w] = d1
    if jh: jh.write(json.dumps({'w': w, 'c': {'%d,%d' % k: v for k, v in d1.items()}, 'c_prime2': {'%d,%d' % k: v for k, v in d2.items()}}) + '\n'); jh.flush()
    if w <= 52 or w % 100 == 0:
        say('   w = %3d: %s  (%.0f s)' % (w, ', '.join('(%s,%s) %d' % (sg(a), sg(b), c[w][(a, b)]) for a, b in SIGNS), time.time() - t0))
say('   the two primes disagree at: %s' % (disagree or 'no weight'))
# a newforms-of-level-6 check: the four classes add up to the new part's dimension
def dimS(N, k):
    import fractions, math
    ps = [p for p in (2, 3) if N % p == 0]
    mu = fractions.Fraction(N)
    for p in ps: mu *= fractions.Fraction(p + 1, p)
    nu2 = 0 if N % 4 == 0 else math.prod(1 + (0 if p == 2 else (1 if p % 4 == 1 else -1)) for p in ps)
    nu3 = 0 if N % 9 == 0 else math.prod(1 + (0 if p == 3 else (1 if p % 3 == 1 else -1)) for p in ps)
    cc = sum(sum(1 for j in range(1, g + 1) if math.gcd(j, g) == 1) for d in range(1, N + 1) if N % d == 0 for g in [math.gcd(d, N // d)])
    gg = 1 + mu / 12 - fractions.Fraction(nu2, 4) - fractions.Fraction(nu3, 3) - fractions.Fraction(cc, 2)
    return int((k - 1) * (gg - 1) + (k // 2 - 1) * cc + nu2 * (k // 4) + nu3 * (k // 3))
bad_sum = [w for w in c if sum(c[w].values()) != dimS(6, w) - 2 * dimS(3, w) - 2 * dimS(2, w) + 4 * dimS(1, w)]
say('   the four classes add up to dim S_w(6) − 2 dim S_w(3) − 2 dim S_w(2) + 4 dim S_w(1) at every weight: %s' % (
    'yes' if not bad_sum else 'NO at %s' % bad_sum[:10]))

say(); say('(a) c(w + 24; ε) − c(w; ε), over every w with w + 24 ≤ %d: %s' % (
    WMAX, sorted({c[w + 24][e] - c[w][e] for w in c if w + 24 in c for e in SIGNS})))
say('(d) c(6; +, +) = %d   (all signs at weight 6: %s)' % (c[6][(1, 1)], ', '.join('(%s,%s) %d' % (sg(a), sg(b), c[6][(a, b)]) for a, b in SIGNS)))
say(); say('(b), (c) the first-step gain at every prime ℓ ≥ 11 (direct: c(2ℓ−4; +,+) − c(ℓ−3; (2|ℓ),(3|ℓ)) where 2ℓ − 4 ≤ %d; '
           'reduced: c(ℓ−5; (2|ℓ),(3|ℓ)) + c(6;+,+)), against ⌊ℓ/24⌋' % WMAX)
rows = []; ok_b = ok_c = n_b = n_c = 0
for l in range(11, WMAX + 6):
    if not is_prime(l): continue
    e = (leg(2, l), leg(3, l))
    red = c[l - 5][e] + c[6][(1, 1)]
    direct = c[2 * l - 4][(1, 1)] - c[l - 3][e] if 2 * l - 4 <= WMAX else None
    loud = c[l - 3][e] if l - 3 in c else None
    if direct is not None: n_b += 1; ok_b += (direct == red)
    n_c += 1; ok_c += (red == l // 24)
    rows.append((l, e, loud, direct, red, l // 24))
for l, e, loud, direct, red, fl in rows:
    if l <= 113 or l in (607,) or red != fl or (direct is not None and direct != red):
        say('   ℓ = %3d: signs (%s,%s); forms at ℓ − 3: %s; gain direct %s, reduced %d; ⌊ℓ/24⌋ = %d%s' % (
            l, sg(e[0]), sg(e[1]), loud, '—' if direct is None else direct, red, fl,
            '' if red == fl and (direct is None or direct == red) else '   DIFFERENT'))
say('   (b) direct = reduced at %d of %d primes;  (c) reduced = ⌊ℓ/24⌋ at %d of %d primes (11 … %d)' % (
    ok_b, n_b, ok_c, n_c, rows[-1][0]))
say(); say('the residue classes of ℓ mod 24 (where (2|ℓ) and (3|ℓ) are fixed): the first prime of each, which with (a) '
           'carries gain = ⌊ℓ/24⌋ to every prime of its class')
for r in (1, 5, 7, 11, 13, 17, 19, 23):
    l0 = next((l for l, *_ in rows if l % 24 == r), None)
    if l0 is None: say('   ℓ ≡ %2d: none in range' % r); continue
    rec = next(t for t in rows if t[0] == l0)
    say('   ℓ ≡ %2d (mod 24): first prime %d, gain %d, ⌊ℓ/24⌋ = %d: %s' % (r, l0, rec[4], rec[5], 'equal' if rec[4] == rec[5] else 'DIFFERENT'))
say(); say('seconds: %.0f' % (time.time() - t0))
