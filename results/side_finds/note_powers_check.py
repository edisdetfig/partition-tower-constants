# A dig, 13 Nov 2026: IS ONE NOTE A POWER OF ANOTHER, IN A PATTERN ACROSS PRIMES? (README section of that name; predicted
# there first). For every prime 37 ≤ ℓ ≤ 500: the class at weight ℓ − 3 modulo ℓ (as partitions_flint_607.py: U₂, U₃ at the
# class's values), T_ℓ on it, its characteristic polynomial factored; the ordinary notes on lane ℓ are (3|ℓ)·(its roots in
# F_ℓ). For each ordered pair of distinct ordinary notes (a, b): whether b is a power of a (b ∈ ⟨a⟩), and if so the exponent
# e with b = a^e, 0 ≤ e < ord(a), as the fraction e/ord(a). Compared with random pairs in the cyclic group of order ℓ − 1:
# the expected share with b ∈ ⟨a⟩ is (1/(ℓ−1)²)·Σ_a |⟨a⟩| = Σ_{d | ℓ−1} φ(d)·d/(ℓ−1)², and the fractions should spread
# evenly. COUNTED (the notes), DERIVED (the expectations).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/note_powers_check.py [--keep] [--max=ℓ]
import os, sys, time, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
import forms6_flint as F
from sympy import primerange, totient, divisors, primefactors
flint.ctx.threads = 2
KEEP = '--keep' in sys.argv
LMAX = next((int(a[6:]) for a in sys.argv[1:] if a.startswith('--max=')), 500)
fh = open(os.path.join(HERE, 'note_powers_check_%d.txt' % LMAX), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def class_constants(L):
    k = L - 3; n = k + 1
    H = F.hecke_matrices(k, L, [('U', 2), ('U', 3), ('T', L)])
    c2 = leg(2, L) * pow(2, (k - 2) // 2, L) % L; c3 = leg(3, L) * pow(3, (k - 2) // 2, L) % L
    K1, r1 = F.generalized_kernel(H[('U', 2)] - F.identity(n, L) * c2); X1, R1 = F.normalized_basis(K1, r1, L)
    M3 = F.restrict(H[('U', 3)], X1, R1, r1); K2, r0 = F.generalized_kernel(M3 - F.identity(r1, L) * c3)
    Xc0 = X1 * flint.nmod_mat(r1, r0, [int(K2[i, j]) for i in range(r1) for j in range(r0)], L)
    Xc, Rc = F.normalized_basis(Xc0, r0, L)
    fac = F.restrict(H[('T', L)], Xc, Rc, r0).charpoly().factor()[1]
    rat = sorted(leg(3, L) * (-int(f.coeffs()[0])) % L for f, e in fac if f.degree() == 1 for _ in range(e))
    return r0, rat, sorted(f.degree() for f, e in fac for _ in range(e))


t0 = time.time(); fracs = []; pairs = 0; inside = 0; expect = 0.0; rows = 0
say('primes 37 to %d: the ordinary notes on lane ℓ, and which are powers of which' % LMAX)
for L in primerange(37, LMAX + 1):
    R, rat, degs = class_constants(L)
    rat = sorted(set(a for a in rat if a))                          # distinct ordinary notes (a zero root would be non-ordinary)
    if len(rat) < 2: continue
    g = next(a for a in range(2, L) if all(pow(a, (L - 1) // q, L) != 1 for q in primefactors(L - 1)))
    lg = {}; x = 1
    for i in range(L - 1): lg[x] = i; x = x * g % L
    n = L - 1; ex = sum(int(totient(d)) * d for d in divisors(n)) / n ** 2
    out = []
    for a in rat:
        for b in rat:
            if a == b: continue
            ta, tb = lg[a], lg[b]; oa = n // math.gcd(ta, n); pairs += 1; expect += ex
            if tb % math.gcd(ta, n) == 0:                           # b = a^e: e·ta ≡ tb (mod n)
                d = math.gcd(ta, n); e = (tb // d) * pow(ta // d, -1, n // d) % (n // d); inside += 1
                fracs.append(e / oa); out.append('%d = %d^%d (of %d)' % (b, a, e, oa))
            else: out.append('%d not a power of %d' % (b, a))
    rows += 1
    say('   ℓ = %3d: %d forms (pieces %s); ordinary notes on lane ℓ: %s;  %s' % (L, R, degs, rat, '; '.join(out)))
say()
say('pairs (ordered): %d over %d primes; b a power of a: %d, against %.1f expected for random pairs (derived)' % (pairs, rows, inside, expect))
if fracs:
    bins = [0] * 10
    for f in fracs: bins[min(9, int(f * 10))] += 1
    chi = sum((b - len(fracs) / 10) ** 2 / (len(fracs) / 10) for b in bins)
    say('the exponents as fractions e/ord(a), in tenths: %s;  chi-squared %.1f on 9 degrees of freedom (about 9 ± 4 if even)' % (bins, chi))
say('seconds: %.0f' % (time.time() - t0))
