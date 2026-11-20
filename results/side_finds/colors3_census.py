# A dig, 12 Nov 2026: 3-COLORED PARTITIONS, THE CENSUS. Direction 1 (the user's choice): is the recipe a general law?
# For p₃(n) (coefficients of ∏(1 − qⁿ)^(−3)) Belmont, Lee, Musat and Trebat-Leder proved that the towers
# p₃((ℓᵇn + 3)/24) settle like the partition tower, with rank at most (QUOTED, their Theorem 1.1, for ℓ ≥ r + 5 = 8)
#     R_ℓ(3) = ⌊k/12⌋ − ⌊3(ℓ² − 1)/(24ℓ)⌋,  k = 2(ℓ − 1)   (k ≢ 2 mod 12 for every odd ℓ),
# which is 0 for 11 and 17: laws there, like 5, 7, 11 for ordinary partitions.
# PREDICTED before the run (derived): the recipe's forms live at weight ℓ − 5 (≡ −4, the Shimura weight of η^(−3)), and since
# η³ is a theta series of 8, not of 24, on LEVEL 2, not 6. The count that should equal R_ℓ(3): the newforms of weight ℓ − 5
# on Γ0(2) with U₂ = (2|ℓ)·2^((ℓ−7)/2), the sign rule of the ordinary partitions' class (there U₂ = (2|ℓ)·2^((k−2)/2) too).
# Shown for comparison: the other sign at level 2, and the four classes at level 6 (ℓ ≤ 151).
# Method: forms2_flint.class_basis_exact2 and forms6_flint.class_basis_exact (exact, over the integers). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors3_census.py [HI] [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import forms2_flint as F2
import forms6_flint as F6
from sympy import primerange
HI = next((int(a) for a in sys.argv[1:] if a.isdigit()), 251); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'colors3_census.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def R(L):
    k = 2 * (L - 1)
    return k // 12 - (1 if k % 12 == 2 else 0) - (3 * (L * L - 1)) // (24 * L)


t0 = time.time(); hits = 0; tot = 0; six = {}
say('ℓ, weight ℓ − 5, R_ℓ(3) (BLMT), level 2: newforms with U₂ = +2^((k−2)/2) and with −2^((k−2)/2); the predicted sign (2|ℓ); level 6 (ℓ ≤ 151): the four classes (U₂, U₃ signs ++, +−, −+, −−)')
for L in primerange(11, HI + 1):
    k = L - 5; r = R(L); d = {s: F2.class_basis_exact2(k, s)[1] for s in (1, -1)}; s = leg(2, L)
    ok = d[s] == r; hits += ok; tot += 1
    line = '   ℓ = %3d, weight %3d: R = %d;  level 2: %d (+), %d (−);  predicted sign %+d: %d  →  %s' % (
        L, k, r, d[1], d[-1], s, d[s], 'as predicted' if ok else 'NOT as predicted')
    if L <= 151:
        c6 = {(a, b): F6.class_basis_exact(k, a, b)[1] for a in (1, -1) for b in (1, -1)}
        line += ';  level 6: %s' % ' '.join('%d' % c6[ab] for ab in ((1, 1), (1, -1), (-1, 1), (-1, -1)))
        for ab in c6: six.setdefault(ab, []).append(c6[ab] == r)
    say(line)
say(); say('level 2 with the predicted sign: %d of %d primes match R_ℓ(3)' % (hits, tot))
for ab, v in sorted(six.items()): say('level 6, signs %+d %+d: %d of %d match' % (ab[0], ab[1], sum(v), len(v)))
say('seconds: %.0f' % (time.time() - t0))
