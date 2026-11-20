# A dig, 11 Nov 2026: DO THE EXTRA FACTORS PILE UP? Where the step rule's multiplier at a prime Q is 0 on lane ℓ, the law
# gains one factor of ℓ (Lovejoy–Ono for 5, Jameson for 7 and 11; partitions_tower.txt, partitions_jameson.txt). Both papers
# stop at one extra factor. Here: (1) whether the step rule also holds on lanes ℓ² and ℓ³ (one multiplier per class), and
# (2) whether two or three zero steps in a row give two or three extra factors, checked with exact values of p(N)
# (Hardy–Ramanujan–Rademacher formula, sympy), with controls in the opposite classes and a contrast where the rule fails.
# A step: 24N − 1 = Q²·(24N* − 1); its class: the sign of 24N* − 1 at Q. Multipliers, members and powers COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_stack.py [--keep]
import os, sys, time, collections
import numpy as np
import sympy as sp
from sympy.functions.combinatorial.numbers import partition
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_stack.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); NMAX = 2_400_000; L = 5 ** 7 * 7 ** 5 * 11 ** 5
G, S, k = [], [], 1
while True:
    a, b = k * (3 * k - 1) // 2, k * (3 * k + 1) // 2
    if a > NMAX: break
    s = 1 if k % 2 else -1; G.append(a); S.append(s)
    if b <= NMAX: G.append(b); S.append(s)
    k += 1
o = np.argsort(G); G = np.array(G, dtype=np.int64)[o]; S = np.array(S, dtype=np.int64)[o]
p = np.zeros(NMAX + 1, dtype=np.int64); p[0] = 1; c = 0
for n in range(1, NMAX + 1):
    while c < len(G) and G[c] <= n: c += 1
    p[n] = int(np.dot(S[:c], p[n - G[:c]])) % L
say('p(n) for every n up to %d on the lane 5⁷·7⁵·11⁵ (%.0f s); p(100) right: %s' % (NMAX, time.time() - t0, int(p[100]) == 190569292 % L))
VL = lambda l, b: b if l in (5, 11) else b // 2 + 1
leg = lambda a, q: int(sp.jacobi_symbol(a % q, q))


def vl(x, l):
    v = 0
    while x % l == 0: x //= l; v += 1
    return v


def square_rule(l, b, M):
    """the classes where c is a multiple of ℓ already (5 at even levels, 7 at odd levels: partitions_tower.txt, section 2)"""
    return (l == 5 and b % 2 == 0 and leg(M // l ** b, 5) == 1) or (l == 7 and b % 2 == 1 and leg(M // l ** b, 7) == 1)


def multipliers(l, b, Q, depth):
    """per class ε = (24N* − 1 | Q): the ratios c(N)/c(N*) on lane ℓ^depth, c = p/ℓ^law, N* at exact level b"""
    g = collections.defaultdict(list); A = l ** b; r0 = pow(24, -1, A); mod = l ** depth
    for Ns in range(r0, (NMAX + (Q * Q - 1) // 24) // (Q * Q) + 1, A):
        M = 24 * Ns - 1
        if vl(M, l) != b or M % Q == 0 or square_rule(l, b, M): continue
        N = Q * Q * Ns - (Q * Q - 1) // 24; a = int(p[Ns]) // l ** VL(l, b); cc = int(p[N]) // l ** VL(l, b)
        if a % l == 0: continue
        g[leg(M, Q)].append(cc * pow(a, -1, mod) % mod)
    return g


say(); say('1  THE STEP RULE ON DEEPER LANES: one multiplier per class (ε = ±1)? Primes Q from 5 to 37 with at least 6 pairs in each class;')
say('   classes where the square rule already makes c a multiple of ℓ are left out')
ZERO = {}
for l, b in ((5, 1), (5, 2), (5, 3), (7, 1), (7, 2), (11, 1), (11, 2)):
    parts = []
    for depth in (1, 2, 3):
        ok = tot = 0
        for Q in sp.primerange(5, 38):
            if Q == l: continue
            g = multipliers(l, b, Q, depth)
            if len(g) < 2 or min(len(v) for v in g.values()) < 6: continue
            tot += 1; ok += all(len(set(v)) == 1 for v in g.values())
            if depth == 1:
                for e_, v in g.items():
                    if set(v) == {0}: ZERO[(l, b, Q)] = e_
        parts.append('lane %d^%d: %d of %d primes' % (l, depth, ok, tot))
    say('   ℓ = %2d, level %d: %s' % (l, b, ';  '.join(parts)))
say('   the zero classes found on lane ℓ (Q: class): ' + ';  '.join('ℓ = %d, level %d: %s' % (l, b, ', '.join(
    '%d: %+d' % (Q, e_) for (l2, b2, Q), e_ in sorted(ZERO.items()) if (l2, b2) == (l, b))) for l, b in ((5, 1), (5, 2), (7, 1), (7, 2), (11, 1), (11, 2))))
CURVE = {7: (1, 1, 1, -4, 5), 11: (1, 0, 1, -14, 20)}                    # Jameson's E7 and E11 (partitions_jameson.txt)


def trace(curve, Q):
    a1, a2, a3, a4, a6 = curve; xs = np.arange(Q, dtype=np.int64); x2 = xs * xs % Q; x3 = x2 * xs % Q
    D = ((a1 * xs + a3) ** 2 + 4 * (x3 + a2 * x2 + a4 * xs + a6)) % Q; sq = np.zeros(Q, dtype=bool); sq[x2] = True
    return -int(np.where(D == 0, 0, np.where(sq[D], 1, -1)).sum())


def zero_class(l, b, Q):
    """the class ε = (24N* − 1 | Q) whose multiplier the theorems (Lovejoy–Ono for 5, Jameson for 7 and 11) make 0 on lane ℓ"""
    for s_ in (1, -1):                                                      # s_ = (m | Q), m = (24N* − 1)/ℓᵇ
        if l == 5: f = 1 + Q - (Q * Q if b % 2 else 1) * leg(-1, Q) * s_
        else:
            cQ = trace(CURVE[l], Q); f = cQ - (s_ if b % 2 else Q ** (3 if l == 7 else 5) * leg(-1, Q) * s_)
        if f % l == 0: return s_ * leg(l, Q) ** b
    return None


agree_ = [(key, e_, zero_class(*key)) for key, e_ in sorted(ZERO.items())]
say('   the same classes from the theorems\' formulas: %d of %d agree' % (sum(1 for _, a_, b_ in agree_ if a_ == b_), len(agree_)))
for l, b, Q in ((11, 1, 67), (11, 2, 67), (5, 1, 43)): ZERO[(l, b, Q)] = zero_class(l, b, Q)
say('   from the formulas, beyond Q = 37: ℓ = 11, Q = 67: class %+d at level 1, %+d at level 2;  ℓ = 5, Q = 43: class %+d at level 1' % (
    ZERO[(11, 1, 67)], ZERO[(11, 2, 67)], ZERO[(5, 1, 43)]))


def run(l, b, Qs, count):
    """N with 24N − 1 = (∏Q)²·M*, M* at exact level b, in the zero class of every Q (and the control: every class flipped)"""
    A = l ** b; r0 = pow(24, -1, A); P2 = 1
    for Q in Qs: P2 *= Q * Q
    law, ctl, Ns = [], [], r0
    while len(law) < count or len(ctl) < count:
        M = 24 * Ns - 1
        if vl(M, l) == b and all(M % Q for Q in Qs) and not square_rule(l, b, M):
            sg = [leg(M, Q) * ZERO[(l, b, Q)] for Q in Qs]
            if all(s_ == 1 for s_ in sg) and len(law) < count: law.append(P2 * Ns - (P2 - 1) // 24)
            elif all(s_ == -1 for s_ in sg) and len(ctl) < count: ctl.append(P2 * Ns - (P2 - 1) // 24)
        Ns += A
    t = time.time(); vals = {N: int(partition(N)) for N in law + ctl}
    pw = [vl(vals[N], l) for N in law]; pc = [vl(vals[N], l) for N in ctl]; big = max(law + ctl)
    say('   ℓ = %2d, level %d (law: %d^%d), Q = %s: powers of %d in p(N): %s;  every class flipped: %s  (largest N = %d, p(N) of %d digits; %.0f s)' % (
        l, b, l, VL(l, b), ' and '.join(map(str, Qs)), l, pw, pc, big, len(str(vals[big])), time.time() - t))
    return pw


say(); say('2  TWO ZERO STEPS IN A ROW (12 members each, the smallest; exact values of p(N))')
for l, b, Qs in ((5, 1, (13, 23)), (5, 2, (13, 23)), (7, 2, (13, 23)), (11, 2, (5, 67)), (7, 1, (13, 23)), (11, 1, (5, 67))):
    if all((l, b, Q) in ZERO for Q in Qs): run(l, b, Qs, 12)
    else: say('   ℓ = %2d, level %d, Q = %s: no zero class for one of them' % (l, b, Qs))
say(); say('3  THREE ZERO STEPS FOR 5 (4 members each)')
if all((5, 1, Q) in ZERO for Q in (13, 23, 43)) and '--quick' not in sys.argv: run(5, 1, (13, 23, 43), 4)
say(); say('seconds: %.0f' % (time.time() - t0))
