# A dig, 11 Nov 2026: WHAT THE STEP RULE IS BEYOND ONE DIGIT. partitions_stack.txt found one multiplier per class on lanes ℓ²
# and ℓ³ (for 5 at levels 1 to 3; for 7 and 11 on lane ℓ² from level 2). Here the two halves of each pair of multipliers,
# μ₊ and μ₋ (classes +1 and −1 of 24N* − 1 at Q): the half-difference is the middle term of a Hecke operator, the half-sum
# its eigenvalue λ. Question: on lanes ℓ, ℓ², ℓ³, is the middle term −(−3|Q)·Q⁻² (Atkin's operator of weight −1/2), and does
# λ keep the simple shape χ(Q)·Qᵃ·(1 + Q^c) that it has on lane 5 for 5 (partitions_tower.txt, section 2)?
# Multipliers COUNTED (primes 7 to 79 with at least 6 pairs in each class); shapes searched over every exponent.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_multipliers.py [--keep]
import os, sys, time, collections
import numpy as np
import sympy as sp
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_multipliers.txt'), 'x', encoding='utf-8') if KEEP else None


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
say('p(n) for every n up to %d on the lane 5⁷·7⁵·11⁵ (%.0f s)' % (NMAX, time.time() - t0))
VL = lambda l, b: b if l in (5, 11) else b // 2 + 1
leg = lambda a, q: int(sp.jacobi_symbol(a % q, q))


def vl(x, l):
    v = 0
    while x % l == 0: x //= l; v += 1
    return v


def mus(l, b, depth):
    out = {}; mod = l ** depth
    for Q in sp.primerange(7, 80):
        if Q == l: continue
        g = collections.defaultdict(list); A = l ** b; r0 = pow(24, -1, A)
        for Ns in range(r0, (NMAX + (Q * Q - 1) // 24) // (Q * Q) + 1, A):
            M = 24 * Ns - 1
            if vl(M, l) != b or M % Q == 0: continue
            if (l == 5 and b % 2 == 0 or l == 7 and b % 2 == 1) and leg(M // l ** b, l) == 1: continue   # square-rule classes
            N = Q * Q * Ns - (Q * Q - 1) // 24; x = int(p[Ns]) // l ** VL(l, b); y = int(p[N]) // l ** VL(l, b)
            if x % l: g[leg(M, Q)].append(y * pow(x, -1, mod) % mod)
        if len(g) == 2 and min(len(v) for v in g.values()) >= 6 and all(len(set(v)) == 1 for v in g.values()):
            out[Q] = (g[1][0], g[-1][0])
    return out


CHI = (('(3|Q)', 3), ('(−3|Q)', -3), ('(−1|Q)', -1), ('(15|Q)', 15), ('(5|Q)', 5), ('1', 1))
say(); say('   ℓ, level, lane: primes with one multiplier per class;  the half-difference (μ₊ − μ₋)/2 = ±χ(Q)·Q^d ?;  the half-sum λ = χ(Q)·Qᵃ·(1 + Q^c) ?')
for l, b, depths in ((5, 1, (1, 2, 3)), (5, 2, (1, 2, 3)), (5, 3, (1, 2)), (7, 2, (1, 2)), (11, 2, (1, 2))):
    for depth in depths:
        mod = l ** depth; MU = mus(l, b, depth); inv2 = pow(2, -1, mod); per = (l - 1) * l ** (depth - 1)
        if len(MU) < 5: continue
        dd = ['%s%s·Q^%d' % ('−' if sg < 0 else '', nm, d) for nm, ch in CHI for sg in (1, -1) for d in range(per)
              if all((MU[Q][0] - MU[Q][1]) * inv2 % mod == sg * leg(ch, Q) * pow(Q, d, mod) % mod for Q in MU)]
        ss = ['%s·Q^%d·(1 + Q^%d)' % (nm, a_, c_) for nm, ch in CHI for a_ in range(per) for c_ in range(per)
              if all((MU[Q][0] + MU[Q][1]) * inv2 % mod == leg(ch, Q) * pow(Q, a_, mod) * (1 + pow(Q, c_, mod)) % mod for Q in MU)]
        say('   ℓ = %2d, level %d, lane %d^%d: %d primes;  half-difference: %s;  half-sum: %s' % (
            l, b, l, depth, len(MU), ', '.join(dd[:2]) or 'no such shape', ', '.join(ss[:2]) or 'no such shape'))
say('   (Q^(per − 2) is Q⁻² on that lane: Q^2 on lane 5, Q^18 on lane 25, Q^98 on lane 125, Q^4 on lane 7, Q^40 on lane 49, …)')
CURVE = {7: (1, 1, 1, -4, 5), 11: (1, 0, 1, -14, 20)}                    # Jameson's E7 and E11 (partitions_jameson.txt)


def trace(curve, Q):
    a1, a2, a3, a4, a6 = curve; xs = np.arange(Q, dtype=np.int64); x2 = xs * xs % Q; x3 = x2 * xs % Q
    D = ((a1 * xs + a3) ** 2 + 4 * (x3 + a2 * x2 + a4 * xs + a6)) % Q; sq = np.zeros(Q, dtype=bool); sq[x2] = True
    return -int(np.where(D == 0, 0, np.where(sq[D], 1, -1)).sum())


say('   Jameson\'s curve on the deeper lane: is λ = ±χ(Q)·Q^a·(1 + Q − #E(F_Q))?')
for l in (7, 11):
    for depth in (1, 2):
        mod = l ** depth; MU = mus(l, 2, depth); inv2 = pow(2, -1, mod); per = (l - 1) * l ** (depth - 1)
        fits = ['%s%s·Q^%d' % ('−' if sg < 0 else '', nm, a_) for nm, ch in CHI for sg in (1, -1) for a_ in range(per)
                if all((MU[Q][0] + MU[Q][1]) * inv2 % mod == sg * leg(ch, Q) * pow(Q, a_, mod) * trace(CURVE[l], Q) % mod for Q in MU)]
        say('   ℓ = %2d, level 2, lane %d^%d: %d primes;  λ = %s' % (l, l, depth, len(MU), ', '.join(fits[:2]) + '·(1 + Q − #E)' if fits else 'no such shape'))
MU = mus(5, 1, 3)
say('   the eigenvalue λ for 5 on lane 125 (level 1), Q = 7 to 79: %s' % ', '.join('%d: %d' % (Q, (MU[Q][0] + MU[Q][1]) * pow(2, -1, 125) % 125) for Q in sorted(MU)))
say(); say('seconds: %.0f' % (time.time() - t0))
