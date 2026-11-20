# A dig, 11 Nov 2026: CAN WE EXTRAPOLATE THE CURVES FOR 17 TO 31? For 7, 11 and 13 the towers' rule on lane ℓ is the newform
# f of level 6 found in partitions_newforms.txt (weight k: 16, 8, 10), and an elliptic curve E has the same numbers:
# a_Q(f) ≡ Q^b·a_Q(E) (mod ℓ) with b = (k − 2)/2 (42a: b ≡ 1; 726i: b = 3; 1014b: b = 4), or that times (Q|ℓ). The power is not a
# free choice. If a curve E over Q had the tower's numbers up to any twist, i.e. its mod-ℓ representation were that of f
# times a character ψ, then comparing determinants (ω for E, ω^(k−1) for f, ω the cyclotomic character mod ℓ) gives
# ψ² = ω^(2−k), so ψ = ω^((2−k)/2)·χ with χ quadratic. Twisting E by χ removes χ (DERIVED). So a curve exists exactly
# when some E over Q has, at every prime Q not dividing 6ℓ,
#     trace of Frobenius ≡ u_Q := Q^(−(k−2)/2)·a_Q(f)  (mod ℓ).
# At such Q the mod-ℓ representation is unramified, so E has good reduction there (trace = 1 + Q − #E(F_Q), with
# |trace| ≤ 2√Q by Hasse's theorem) or multiplicative reduction (trace ≡ ±(Q + 1)); additive reduction would make it
# ramified (ℓ ≥ 5). So if at one prime Q the residue u_Q is neither a whole number of size at most 2√Q nor ±(Q + 1) on
# lane ℓ, no elliptic curve over Q of ANY conductor can carry the tower's rule. That is checked here at every prime up to 397.
# The newforms are recomputed as in partitions_newforms.py; which one is each tower's rule is read from its .txt.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_nocurve.py [--keep]
import os, re, sys, math, time
import numpy as np
import sympy as sp
from sympy import QQ
from sympy.polys.matrices import DomainMatrix
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_nocurve.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); LEN = 400
E2 = [1] + [-24 * int(sp.divisor_sigma(n)) for n in range(1, LEN)]
def scale(f, d): g = [0] * LEN; g[::d] = f[:(LEN - 1) // d + 1]; return g
GENS = [[a - d * b for a, b in zip(E2, scale(E2, d))] for d in (2, 3, 6)]
def mul(f, g): return list(np.convolve(np.array(f, dtype=object), np.array(g, dtype=object))[:LEN])


def nullspace(rows):
    ns = DomainMatrix([[QQ(int(x)) for x in r] for r in rows], (len(rows), len(rows[0])), QQ).nullspace().to_Matrix(); out = []
    for i in range(ns.rows):
        v = [sp.Rational(x) for x in ns.row(i)]; den = sp.ilcm(*[x.q for x in v]) if len(v) > 1 else v[0].q
        out.append([int(x * den) for x in v])
    return out


def newform(k, a2, a3):
    """the newform of weight k and level 6 with a₂ = a2, a₃ = a3 (U₂, U₃ eigenvalues), from products of E₂(z) − d·E₂(dz)"""
    cur = {(0, 0, 0): [1] + [0] * (LEN - 1)}
    for _ in range(k // 2):
        nxt = {}
        for key, f in cur.items():
            for i in range(3):
                nk = list(key); nk[i] += 1; nk = tuple(nk)
                if nk not in nxt: nxt[nk] = mul(f, GENS[i])
        cur = nxt
    P = 2 ** 61 - 1; chosen, red = [], []
    for f in cur.values():
        v = [x % P for x in f[:3 * (k + 1) + 20]]
        for piv, r in red:
            if v[piv]: c = v[piv] * pow(r[piv], -1, P) % P; v = [(a - c * b) % P for a, b in zip(v, r)]
        nz = next((i for i, x in enumerate(v) if x), None)
        if nz is not None: red.append((nz, v)); chosen.append(f)
        if len(chosen) == k + 1: break
    ns = nullspace([[b[2 * n] - a2 * b[n] for b in chosen] for n in range(min(LEN // 2, 4 * k + 20))])
    B2 = [[sum(v[i] * chosen[i][n] for i in range(len(chosen))) for n in range(LEN)] for v in ns]
    ns3 = nullspace([[b[3 * n] - a3 * b[n] for b in B2] for n in range(min(LEN // 3, 4 * k + 20))])
    assert len(ns3) == 1, 'not one-dimensional'
    f = [sum(v * B2[i][n] for i, v in enumerate(ns3[0])) for n in range(LEN)]
    assert f[1] != 0 and all(x % f[1] == 0 for x in f)
    return [x // f[1] for x in f]


txt = open(os.path.join(HERE, 'partitions_newforms.txt'), encoding='utf-8').read()
WHICH = {int(l): (int(k), int(a2), int(a3)) for l, k, a2, a3 in re.findall(
    r'ℓ = +(\d+), level 1 \(\d+ primes\): weight (\d+), a₂ = ([+-]\d+), a₃ = ([+-]\d+)', txt)}
say('the towers\' rules, read from partitions_newforms.txt (level 1, first match): %s' % ', '.join(
    'ℓ = %d: weight %d (a₂ = %+d, a₃ = %+d)' % (l, *WHICH[l]) for l in sorted(WHICH)))


def trace(curve, Q):
    a1, a2, a3, a4, a6 = curve; xs = np.arange(Q, dtype=np.int64); x2 = xs * xs % Q; x3 = x2 * xs % Q
    D = ((a1 * xs + a3) ** 2 + 4 * (x3 + a2 * x2 + a4 * xs + a6)) % Q; sq = np.zeros(Q, dtype=bool); sq[x2] = True
    return -int(np.where(D == 0, 0, np.where(sq[D], 1, -1)).sum())


CURVES = {7: ('42a1', (1, 1, 1, -4, 5)), 11: ('726i1', (1, 0, 0, -1636, -28588)), 13: ('1014b1', (1, 0, 1, -10482, 722308))}
say(); say('THE TEST at every prime 5 ≤ Q ≤ 397, Q not dividing 6ℓ: is u_Q = Q^(−(k−2)/2)·a_Q(f) on lane ℓ a whole number of size at')
say('most 2√Q, or ±(Q + 1)? "fingerprint": the only possible traces of a curve at the first primes (several if there is a choice)')
for l in sorted(WHICH):
    k, a2, a3 = WHICH[l]; f = newform(k, a2, a3); b0 = (k - 2) // 2; obst = []; fp = []
    for Q in sp.primerange(5, LEN):
        if Q == l: continue
        u = f[Q] * pow(Q, -b0, l) % l; hb = math.isqrt(4 * Q)
        good = [a for a in range(-hb, hb + 1) if (a - u) % l == 0]; mult = [s * (Q + 1) for s in (1, -1) if (s * (Q + 1) - u) % l == 0]
        if len(fp) < 8: fp.append('%d: %s' % (Q, '/'.join(map(str, good + ['±(Q+1)'] * (1 if mult else 0))) or 'none'))
        if not good and not mult: obst.append(Q)
    extra = ''
    if l in CURVES:
        lab, cur = CURVES[l]
        agree = all((trace(cur, Q) - f[Q] * pow(Q, -b0, l)) % l == 0 for Q in sp.primerange(5, LEN) if Q != l and 6 * l * 0 == 0)
        extra = ';  %s has exactly these traces at every prime up to 397: %s' % (lab, agree)
    say('   ℓ = %2d (weight %d, power Q^-%d): %s%s' % (l, k, b0, 'no obstruction at any prime up to 397' if not obst else
        'IMPOSSIBLE at Q = %s (%d primes of %d)' % (', '.join(map(str, obst[:6])), len(obst), sum(1 for Q in sp.primerange(5, LEN) if Q != l)), extra))
    say('      fingerprint: %s' % ';  '.join(fp))
say(); say('seconds: %.0f' % (time.time() - t0))
