# A dig, 12 Nov 2026: THE LAW AT Q = 2999 FOR 37, CHECKED THROUGH ITS CAUSE. partitions_split.txt found Q = 2999 to be the one
# prime below 10,000 where both of 37's twin numbers vanish on lane 37 (t ≡ n ≡ 0). DERIVED from the step rule, that gives
#     p(N) ≡ 0 (mod 37) whenever 24N − 1 = 37·2999³·m with 2999 ∤ m
# (an Ono-type congruence in the literature's naming: Yang, IMRN 2011, Remark 3.2 and Theorem 7). Its smallest member,
# p(540,584,055,355), is beyond sympy's partition() (which refuses n > 149,832,547,102), and PARI's numbpart was found wrong
# once at a tenth of that size (partitions_twins_43_sympy.py). So the cause is checked instead. Both numbers are 0 exactly when
# the step at 2999 multiplies by ε·δ alone, for every M with 37 | M and 2999 ∤ M:
#     c(2999²·M) ≡ ε·δ·c(M) (mod 37),  c(M) = p(N) with M = 24N − 1,  ε = (M | 2999),  δ = −(−3|2999)·2999⁻²,
# while if either number were not 0 this would fail for most M. PREDICTED before the run: it holds for every M tested (the
# first 8 values of M with 37 | M, 2999 ∤ M, both classes; exact values from sympy). A trial with the direct member stopped at
# sympy's limit. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_split_2999.py [--keep]
import os, sys, time
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
from sympy import jacobi_symbol
from sympy.functions.combinatorial.numbers import partition
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_split_2999.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); L, Q = 37, 2999
leg = lambda a, q: int(jacobi_symbol(a % q, q))
dl = (-leg(-3, Q) * pow(Q, -2, L)) % L
say('ℓ = %d, Q = %d: δ = −(−3|Q)·Q⁻² ≡ %d on lane %d' % (L, Q, dl, L))
d1 = pow(24, -1, L); ok = tot = 0; j = 0; seen = {1: 0, -1: 0}
while tot < 8:
    Ns = L * j + d1; M = 24 * Ns - 1; j += 1
    if M % Q == 0: continue
    e = leg(M, Q); N = Q * Q * Ns - (Q * Q - 1) // 24
    t = time.time(); c0, c1 = int(partition(Ns)) % L, int(partition(N)) % L; want = e * dl * c0 % L
    tot += 1; ok += c1 == want; seen[e] += 1
    say('   M = 24·%d − 1 (class %+d): p(%d) ≡ %d;  ε·δ·p(%d) = %+d·%d·%d ≡ %d  →  %s  (%.0f s)' % (
        Ns, e, N, c1, Ns, e, dl, c0, want, 'equal' if c1 == want else 'NOT equal', time.time() - t))
say('%d of %d equal (classes +1: %d, −1: %d)  →  %s' % (ok, tot, seen[1], seen[-1],
    'both of 37\'s numbers are 0 at 2999, so the law with 2999³ holds (derived from the step rule)' if ok == tot else 'not both 0'))
say(); say('seconds: %.0f' % (time.time() - t0))
