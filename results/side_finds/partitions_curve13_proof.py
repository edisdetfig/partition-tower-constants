# A dig, 11 Nov 2026: THE CURVE BEHIND 13, PROVED THE WAY JAMESON PROVED HERS. Claim (the 13-analogue of her Theorem 1.1):
# for every prime Q ≥ 5, Q ≠ 13, and every n ≥ 0,
#   p(13Q²n + β) ≡ (3|Q)·Q⁴·[a_Q − ((−24n − 11)|Q)]·p(13n + 6) − Q⁹·p(13n/Q² + β′)   (mod 13),
#   β = (143Q² + 1)/24,  a_Q = 1 + Q − #E(F_Q),  E = 1014b1: y² + xy + y = x³ − 10482x + 722308,
# the last term only when Q² divides 24n + 11. Her exponents were Q¹ (7) and Q³ (11); this one is Q⁴ = Q^((ℓ−5)/2).
# The chain, each link checked here:
#   (0) p(13n + 6) ≡ 11·[qⁿ]∏(1 − qᵏ)¹¹ (mod 13): classical; checked again on every n with 13n + 6 ≤ 2,400,000.
#   (1) g = η(24z)¹¹ = Σ a(m)qᵐ is a Hecke eigenform (weight 11/2, character χ₁₂); its eigenvalue at Q is
#       λ(Q) = a(11Q²) + (3|Q)·(−11|Q)·Q⁴ (a(11) = 1). Its Shimura lift F has weight 10 and level 288 (Niwa), and
#       F ≡ f ⊗ χ₁₂ (mod 13) for the newform f of weight 10 and level 6 (Ahlgren–Allen–Tang's form) follows, by Sturm's
#       theorem in S₁₀(Γ0(864)), from λ(Q) ≡ (3|Q)·a_Q(f) at every prime 5 ≤ Q ≤ 1440 (both sides vanish at n not prime to 6,
#       and both are eigenforms, so primes suffice, as in Jameson's proof).
#   (2) f ≡ Θ⁴(E's newform) (mod 13) at every n prime to 78: with the 2-, 3- and 13-parts removed, both are forms of weight 58
#       on Γ0(6084) modulo 13 (f·E₁₂⁴, E₁₂ ≡ 1; Θ raises the weight by 14 through the Serre derivative), so Sturm's bound
#       58·13104/12 = 63,336 decides. Both are eigenforms with the same recursion at prime powers (p⁹ = p⁴·p⁴·p), so it is
#       enough that a_p(f) ≡ p⁴·a_p(E) (mod 13) at every prime 5 ≤ p ≤ 63,336, p ≠ 13.
#   (3) Then (1) and (2) give λ(Q) ≡ (3|Q)·Q⁴·a_Q(E), and the Hecke relation of g gives the claim; checked once more
#       directly on partition numbers.
# Standard theorems used and not re-proved: Shimura/Niwa (the lift and its level), Sturm (bounds), the eigenform property of
# g (Ahlgren–Allen–Tang: these spaces are one-dimensional for ℓ ≤ 31), the Θ operator on forms mod ℓ (Serre), modularity
# of E (a_p from point counts). Checks COUNTED (exact integers, FFT products mod 13 with exact rounding).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_curve13_proof.py [--keep]
import os, sys, time
import numpy as np
import sympy as sp
from sympy import QQ
from sympy.polys.matrices import DomainMatrix
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_curve13_proof.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); M13 = 13; E13B = (1, 0, 1, -10482, 722308); E13F = (1, 0, 0, -62, 324)
leg = lambda a, q: int(sp.jacobi_symbol(a % q, q))


def conv(a, b, L):
    """product of two series mod 13, first L terms; FFT in floating point, exact after rounding (sums < 13²·L)"""
    n = 1 << (2 * L - 1).bit_length()
    c = np.fft.irfft(np.fft.rfft(a.astype(np.float64), n) * np.fft.rfft(b.astype(np.float64), n), n)[:L]
    r = np.rint(c); assert np.abs(c - r).max() < 0.2
    return r.astype(np.int64) % M13


def trace(curve, Q):
    a1, a2, a3, a4, a6 = curve; xs = np.arange(Q, dtype=np.int64); x2 = xs * xs % Q; x3 = x2 * xs % Q
    D = ((a1 * xs + a3) ** 2 + 4 * (x3 + (a2 % Q) * x2 + (a4 % Q) * xs + (a6 % Q))) % Q; sq = np.zeros(Q, dtype=bool); sq[x2] = True
    return -int(np.where(D == 0, 0, np.where(sq[D], 1, -1)).sum())


# ---------- the series ∏(1 − qᵏ)¹¹ mod 13 ----------
QA = 1440; LG = 11 * (QA * QA - 1) // 24 + 2
E = np.zeros(LG, dtype=np.int64); E[0] = 1; k = 1
while k * (3 * k - 1) // 2 < LG:
    s = 1 if k % 2 == 0 else -1
    for e_ in (k * (3 * k - 1) // 2, k * (3 * k + 1) // 2):
        if e_ < LG: E[e_] += s
    k += 1
E %= M13; E2 = conv(E, E, LG); E4 = conv(E2, E2, LG); E8 = conv(E4, E4, LG); E11 = conv(conv(E8, E2, LG), E, LG)
say('∏(1 − qᵏ)¹¹ mod 13 to q^%d (%.0f s)' % (LG - 1, time.time() - t0))

say(); say('(0) PARTITIONS TO g: p(13n + 6) ≡ 11·[qⁿ]∏(1 − qᵏ)¹¹ (mod 13), every n with 13n + 6 ≤ 2,400,000')
NMAX = 2_400_000
G, S, kk = [], [], 1
while True:
    x1, x2 = kk * (3 * kk - 1) // 2, kk * (3 * kk + 1) // 2
    if x1 > NMAX: break
    s = 1 if kk % 2 else -1; G.append(x1); S.append(s)
    if x2 <= NMAX: G.append(x2); S.append(s)
    kk += 1
o = np.argsort(G); G = np.array(G, dtype=np.int64)[o]; S = np.array(S, dtype=np.int64)[o]
p = np.zeros(NMAX + 1, dtype=np.int64); p[0] = 1; c = 0
for n in range(1, NMAX + 1):
    while c < len(G) and G[c] <= n: c += 1
    p[n] = int(np.dot(S[:c], p[n - G[:c]])) % M13
ns = np.arange(0, (NMAX - 6) // 13 + 1)
say('   %d values, wrong: %d' % (len(ns), int(((p[13 * ns + 6] - 11 * E11[ns]) % M13 != 0).sum())))

# ---------- the newform f of weight 10, level 6, mod 13 ----------
LB = 63_337; LX = 90
def sig_table(L):
    sg = np.zeros(L, dtype=np.int64)
    for d in range(1, L): sg[d::d] += d
    return sg
sg_exact = sig_table(LX); sg13 = sig_table(LB) % M13


def e2_series(sg, L, exact):
    if exact: return np.array([1] + [-24 * int(x) for x in sg[1:L]], dtype=object)     # Python integers: no overflow
    e = np.zeros(L, dtype=np.int64); e[0] = 1; e[1:] = -24 * sg[1:L]
    return e


def gens(L, exact):
    sg = sg_exact if exact else sg13; base = e2_series(sg, L, exact); out = []
    for d in (2, 3, 6):
        sc = np.zeros(L, dtype=base.dtype); sc[::d] = base[:(L - 1) // d + 1]; g = base - d * sc
        out.append(g if exact else g % M13)
    return out


def monomials(gs, L, mult):
    cur = {(0, 0, 0): np.eye(1, L, dtype=gs[0].dtype)[0]}
    for _ in range(5):
        nxt = {}
        for key, f in cur.items():
            for i in range(3):
                nk = list(key); nk[i] += 1; nk = tuple(nk)
                if nk not in nxt: nxt[nk] = mult(f, gs[i])
        cur = nxt
    return cur


ex = monomials(gens(LX, True), LX, lambda a, b: np.convolve(a, b)[:LX])
allkeys = sorted(ex); P61 = 2 ** 61 - 1; keys, red = [], []          # 11 independent products, picked by elimination mod P
for kk_ in allkeys:
    v = [int(x) % P61 for x in ex[kk_][:60]]
    for piv, rr in red:
        if v[piv]: cc = v[piv] * pow(rr[piv], -1, P61) % P61; v = [(a - cc * b) % P61 for a, b in zip(v, rr)]
    nz = next((i for i, x in enumerate(v) if x), None)
    if nz is not None: red.append((nz, v)); keys.append(kk_)
assert len(keys) == 11, len(keys)
Mx = [[int(ex[kk_][n]) for kk_ in keys] for n in range(LX)]
# the newform: U2 eigenvalue −16, U3 eigenvalue +81 (a₂ = −16, a₃ = 81; Ahlgren–Allen–Tang's q − 16q² + 81q³ + …)
rows = [[Mx[2 * n][j] + 16 * Mx[n][j] for j in range(len(keys))] for n in range(LX // 2)] + \
       [[Mx[3 * n][j] - 81 * Mx[n][j] for j in range(len(keys))] for n in range(LX // 3)]
nsp = DomainMatrix([[QQ(x) for x in r] for r in rows], (len(rows), len(keys)), QQ).nullspace().to_Matrix()
assert nsp.rows == 1, nsp.rows                                          # one newform with these U₂, U₃ eigenvalues
fvec = [sum(sp.Rational(nsp[0, j]) * Mx[n][j] for j in range(len(keys))) for n in range(LX)]
coef = [sp.Rational(nsp[0, j]) / fvec[1] for j in range(len(keys))]; fx = [x / fvec[1] for x in fvec]
assert all(x.q == 1 for x in fx)
say(); say('THE NEWFORM f OF WEIGHT 10 AND LEVEL 6: %s …  (Ahlgren–Allen–Tang print q − 16q² + 81q³ + 256q⁴ + 2694q⁵)' % ' '.join('%+d·q^%d' % (int(fx[n]), n) for n in range(1, 8)))
dens = [x.q for x in coef]; bad13 = [d_ for d_ in dens if d_ % 13 == 0]
say('   written in the 21 products of E₂(z) − d·E₂(dz) (d = 2, 3, 6); denominators prime to 13: %s' % (not bad13))
md = monomials(gens(LB, False), LB, lambda a, b: conv(a, b, LB))
f13 = np.zeros(LB, dtype=np.int64)
for j, kk_ in enumerate(keys):
    if coef[j] != 0: f13 = (f13 + (int(coef[j].p) * pow(int(coef[j].q), -1, M13)) % M13 * md[kk_]) % M13
assert all((f13[n] - int(fx[n])) % M13 == 0 for n in range(LX)), 'mod-13 series disagrees with the exact one'
say('   f mod 13 to q^%d (%.0f s); its first %d coefficients agree with the exact ones' % (LB - 1, time.time() - t0, LX))

say(); say('(1) THE TOWER\'S EIGENFORM AGAINST f: λ(Q) = a(11Q²) + (3|Q)(−11|Q)·Q⁴ ≡ (3|Q)·a_Q(f) (mod 13), every prime 5 ≤ Q ≤ 1440')
bad1 = []; n1 = 0
for Q in sp.primerange(5, QA + 1):
    if Q == 13: continue
    lam = (int(E11[11 * (Q * Q - 1) // 24]) + leg(3, Q) * leg(-11, Q) * pow(Q, 4, M13)) % M13
    n1 += 1
    if lam != leg(3, Q) * int(f13[Q]) % M13: bad1.append(Q)
say('   %d primes, wrong: %d %s  (Sturm bound for S₁₀(Γ0(864)): 10·1728/12 = 1440)' % (n1, len(bad1), bad1[:5] if bad1 else ''))
say('   also, λ(Q) from partition numbers (the ratio p(13Q²n + β)/p(13n + 6) in each class, as in partitions_curves.txt)')
say('   is the same λ: checked at the primes 5 to 113 below, in (3)')

say(); say('(2) f AGAINST THE CURVE: a_p(f) ≡ p⁴·a_p(E) (mod 13), E = 1014b1, every prime 5 ≤ p ≤ 63,336, p ≠ 13')
t = time.time(); bad2 = []; n2 = 0
for Q in sp.primerange(5, LB):
    if Q == 13: continue
    n2 += 1
    if (int(f13[Q]) - pow(Q, 4, M13) * trace(E13B, Q)) % M13: bad2.append(Q)
say('   %d primes, wrong: %d %s  (%.0f s; Sturm bound for weight 58 on Γ0(6084): 58·13104/12 = 63,336)' % (n2, len(bad2), bad2[:5] if bad2 else '', time.time() - t))
tw = all(trace(E13B, Q) == leg(Q, 13) * trace(E13F, Q) for Q in sp.primerange(5, 2000) if Q != 13)
say('   1014b1 is the twist of 1014f1 by 13 (a_p(1014b) = (p|13)·a_p(1014f) for every prime p < 2000): %s; so with 1014f1 the' % tw)
say('   same statement reads Q¹⁰·a_Q(1014f), since (Q|13) ≡ Q⁶ (mod 13)')

say(); say('(3) THE STATEMENT ITSELF on partition numbers: p(13Q²n + β) ≡ (3|Q)·Q⁴·[a_Q − ((−24n − 11)|Q)]·p(13n + 6) − Q⁹·p(N**)')
checks = wrong = nq = 0
for Q in sp.primerange(5, 400):
    if Q == 13: continue
    beta = (143 * Q * Q + 1) // 24; aq = trace(E13B, Q); nq_any = False
    for n in range(0, (NMAX - beta) // (13 * Q * Q) + 1):
        m = 24 * n + 11; rhs = leg(3, Q) * pow(Q, 4, M13) * (aq - leg(-m, Q)) * int(p[13 * n + 6])
        if m % (Q * Q) == 0: rhs -= pow(Q, 9, M13) * int(p[(13 * (m // (Q * Q)) + 1) // 24])
        checks += 1; nq_any = True; wrong += (int(p[13 * Q * Q * n + beta]) - rhs) % M13 != 0
    nq += nq_any
say('   %d primes Q, %d values of n, wrong: %d' % (nq, checks, wrong))

say(); say('WHAT IT GIVES (her corollaries 1.2 and 1.3, for 13). With i = a_Q mod 13:')
say('   i = ±1: p(13Q²(QN + r) + β) ≡ 0 (mod 13) for every r with ((−24r − 11)|Q) = i;')
say('   i = 0:  p(13Q²(Q²N + Qs + r) + β) ≡ 0 (mod 13) for r with Q | 24r + 11 and s with Q² ∤ 24(r + Qs) + 11.')
tab = {i: [Q for Q in sp.primerange(5, 200) if Q != 13 and trace(E13B, Q) % 13 == i % 13] for i in (0, 1, -1)}
say('   primes below 200 with i = 0: %s;  i = +1: %s;  i = −1: %s' % tuple(', '.join(map(str, tab[i])) for i in (0, 1, -1)))
Q = 59; beta = (143 * Q * Q + 1) // 24; r = (-11 * pow(24, -1, Q)) % Q
say('   Q = 59 (i = 0): β = %d, r = %d (24r + 11 = %d), s = 0: p(13·59²·(59²N + %d) + %d) = p(13·59⁴·N + %d), Atkin and O\'Brien\'s and Ono\'s example' % (
    beta, r, 24 * r + 11, r, beta, 13 * Q * Q * r + beta))
say(); say('seconds: %.0f' % (time.time() - t0))
