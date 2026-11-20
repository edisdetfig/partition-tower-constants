# A dig, 12 Nov 2026: 5-COLORED PARTITIONS, THE CENSUS. After 3 colours (colors3_census.txt: level 2, weight ℓ − 5, exactly
# R_ℓ(3) forms for 50 of 50 primes), the next test of the recipe's generality. Belmont, Lee, Musat and Trebat-Leder's bound
# (QUOTED, Theorem 1.1, ℓ ≥ r + 5 = 10): R_ℓ(5) = ⌊k/12⌋ − [k ≡ 2 mod 12] − ⌊5(ℓ² − 1)/(24ℓ)⌋, k = 3(ℓ − 1).
# PREDICTED before the run (derived): η^5 has no theta-series shortcut (5 is prime to 3), so the forms are back on LEVEL 6,
# at weight ℓ − 7 (≡ −6, the Shimura weight of η^(−5)), in one Atkin–Lehner class given by a sign rule of the form
# U₂ = ε₂·(2|ℓ)·2^((k−2)/2), U₃ = ε₃·(3|ℓ)·3^((k−2)/2) with fixed ε₂, ε₃ = ±1, whose dimension is R_ℓ(5) for every ℓ.
# Shown: all four rules, and level 2 for comparison. Method: forms6_flint.class_basis_exact, forms2_flint.class_basis_exact2
# (exact, over the integers). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors5_census.py [HI] [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import forms6_flint as F6
import forms2_flint as F2
from sympy import primerange
HI = next((int(a) for a in sys.argv[1:] if a.isdigit()), 157); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'colors5_census.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def R(L, r=5):
    k = (r + 1) // 2 * (L - 1)
    return k // 12 - (1 if k % 12 == 2 else 0) - (r * (L * L - 1)) // (24 * L)


t0 = time.time(); fits = {}
say('ℓ, weight ℓ − 7, R_ℓ(5) (BLMT); level 6: the class under each rule (ε₂, ε₃) = (+,+), (+,−), (−,+), (−,−); level 2 (U₂ sign (2|ℓ), −(2|ℓ))')
for L in primerange(11, HI + 1):
    k = L - 7; r = R(L); s2, s3 = leg(2, L), leg(3, L); row = []
    for e2, e3 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        dd = F6.class_basis_exact(k, e2 * s2, e3 * s3)[1]; row.append(dd); fits.setdefault(('6', e2, e3), []).append(dd == r)
    two = [F2.class_basis_exact2(k, e * s2)[1] for e in (1, -1)]
    for e, dd in zip((1, -1), two): fits.setdefault(('2', e, 0), []).append(dd == r)
    say('   ℓ = %3d, weight %3d: R = %d;  level 6: %s;  level 2: %s' % (L, k, r, ' '.join(map(str, row)), ' '.join(map(str, two))))
say()
for key, v in sorted(fits.items()):
    say('   level %s, rule %s: %d of %d primes match R_ℓ(5)' % (key[0], '(%+d, %+d)' % (key[1], key[2]) if key[0] == '6' else '(%+d)' % key[1], sum(v), len(v)))
say('seconds: %.0f' % (time.time() - t0))
