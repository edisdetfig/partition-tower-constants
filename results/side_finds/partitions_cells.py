# A dig, 11 Nov 2026: WHERE THE EXTRA MULTIPLES COME FROM. partitions_bias.txt: beyond every law, p(N) is a multiple of 5,
# 7 and 13 a little more often than 1 in 5, 7, 13 (not of 11, 17, 19, 23, 29, 31), and the extra multiples crowd into the
# classes where 24N − 1 holds a square Q². Here: write 24N − 1 = Q²·m with Q not dividing m, and N* for the smaller number
# with 24N* − 1 = m.
#   1  Is p(N) ≡ μ·p(N*) (mod ℓ), with one μ for each class? Classes: the sign of −m on lane ℓ (Radu proved in 2013 that a law
#      p(An + B) ≡ 0 (mod ℓ) can only stand where that sign is −1 or where ℓ divides m: here "Radu's class" is the −1 half),
#      and the sign of m on lane Q.
#   2  A formula for μ, fitted on the primes Q below 140, then checked on every prime Q up to 7,589 (the largest whose square
#      fits under 24·2,400,000).
#   3  Where μ is 0 there is a law. Atkin proved these in 1968 (Proc. London Math. Soc. 18, 563–576, Theorems 1 and 2, as
#      quoted by Ahlgren, Allen and Tang 2022, eq. (1.7)): for ℓ = 5, 7, 13 and every prime Q ≡ −2 (mod ℓ). Found here from
#      the counts, then checked.
#   4  How much of the extra multiples these laws account for.
#   5  What is left over.
# p(n) for every n up to 2,400,000 by Euler's pentagonal rule, whole numbers only, kept on the lane 5·7·11·13·17·19·23.
# Shares, members and failures are COUNTED; the μ formula is FITTED (on what: section 2); predictions and z are DERIVED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_cells.py [--keep]
import os, sys, time, math
import numpy as np
import sympy as sp
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_cells.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); N = 2_400_000; LS = [5, 7, 11, 13, 17, 19, 23]; L = math.prod(LS)
G, S, k = [], [], 1
while True:
    a, b = k * (3 * k - 1) // 2, k * (3 * k + 1) // 2
    if a > N: break
    s = 1 if k % 2 else -1; G.append(a); S.append(s)
    if b <= N: G.append(b); S.append(s)
    k += 1
o = np.argsort(G); G = np.array(G, dtype=np.int64)[o]; S = np.array(S, dtype=np.int64)[o]
p = np.zeros(N + 1, dtype=np.int64); p[0] = 1; c = 0
for n in range(1, N + 1):
    while c < len(G) and G[c] <= n: c += 1
    p[n] = int(np.dot(S[:c], p[n - G[:c]])) % L
say('p(n) for every n up to %d, kept on the lane 5·7·11·13·17·19·23 (%.0f s); p(100) comes out right: %s' % (N, time.time() - t0, int(p[100]) == 190569292 % L))

nn = np.arange(N + 1, dtype=np.int64); MM = 24 * nn - 1; QMAX = math.isqrt(24 * N - 1); PRIMES = list(sp.primerange(5, QMAX + 1))
_tabs = {}


def leg(x, q):
    """the sign of x on lane q (Legendre symbol), for a number or an array"""
    if q not in _tabs: _tabs[q] = np.array([0] + [1 if pow(r, (q - 1) // 2, q) == 1 else -1 for r in range(1, q)], dtype=np.int64)
    return _tabs[q][np.asarray(x) % q]


def square_members(Q):
    """every N ≤ N with Q² dividing 24N − 1 exactly once over (Q not dividing m); returns N, m, N*"""
    r = pow(24, -1, Q * Q); Ns = np.arange(r if r else Q * Q, N + 1, Q * Q, dtype=np.int64)
    m = (24 * Ns - 1) // (Q * Q); keep = m % Q != 0
    return Ns[keep], m[keep], (m[keep] + 1) // 24


MEM = {Q: square_members(Q) for Q in PRIMES}

say(); say('1  IS p(N) ≡ μ·p(N*) (mod ℓ) WITH ONE μ PER CLASS?  (every prime Q from 5 to 139 whose classes hold at least 12 members each)')
for l in LS:
    inv = np.array([0] + [pow(r, -1, l) for r in range(1, l)], dtype=np.int64); pl = p % l; got = {'Radu': [], 'other': [], 'ℓ | m': []}; tested = {key: 0 for key in got}
    for Q in [q for q in PRIMES if q < 140 and q != l]:
        Ns, m, Nst = MEM[Q]; e2 = leg(-m, l); e1 = leg(m, Q); a = pl[Ns]; b = pl[Nst]
        for key, E2 in (('Radu', -1), ('other', 1), ('ℓ | m', 0)):
            ok = True; enough = True
            for E1 in (1, -1):
                cell = (e2 == E2) & (e1 == E1)
                if cell.sum() < 12: enough = False; break
                aa, bb = a[cell], b[cell]; nz = bb != 0
                ok &= len(set(((aa[nz] * inv[bb[nz]]) % l).tolist())) == 1 and bool(np.all(aa[~nz] == 0))
            if enough: tested[key] += 1; got[key] += [Q] if ok else []
    def show(key):
        if key == 'ℓ | m' and l in (5, 7, 11): return 'a law already (p(N) is always a multiple)'
        g = got[key]; t = tested[key]
        return '%d of %d%s' % (len(g), t, (' (Q = %s)' % ', '.join(map(str, g))) if 0 < len(g) < t else '')
    say('   ℓ = %2d:  Radu\'s class: %s;  the other class: %s;  ℓ | m: %s' % (l, show('Radu'), show('other'), show('ℓ | m')))

say(); say('2  A FORMULA FOR μ IN RADU\'S CLASS, ℓ = 5, 7, 13. Fitted on Q < 140: μ ≡ Q⁻²·[ (3|Q)·(1 + Q) − (−3|Q)·(m|Q) ]  (mod ℓ),')
say('   where (a|Q) is the sign of a on lane Q. Checked on every member, for every prime Q up to %d:' % QMAX)


def mu_formula(l, Q, m):
    s3 = int(leg(3, Q)); sm3 = int(leg(-3, Q)); q2 = pow(Q * Q % l, -1, l)
    return (q2 * (s3 * (1 + Q) - sm3 * leg(m, Q))) % l


for l in (5, 7, 13):
    pl = p % l; tot = {True: [0, 0, 0], False: [0, 0, 0]}
    for Q in PRIMES:
        if Q == l: continue
        Ns, m, Nst = MEM[Q]; cell = leg(-m, l) == -1
        if not cell.any(): continue
        mu = mu_formula(l, Q, m[cell]); bad = int(((pl[Ns[cell]] - mu * pl[Nst[cell]]) % l != 0).sum())
        t = tot[Q < 140]; t[0] += 1; t[1] += int(cell.sum()); t[2] += bad
    say('   ℓ = %2d:  Q < 140 (where it was fitted): %d primes, %d members, %d wrong;  140 < Q ≤ %d (not seen in the fit): %d primes, %d members, %d wrong' % (
        l, tot[True][0], tot[True][1], tot[True][2], QMAX, tot[False][0], tot[False][1], tot[False][2]))

say(); say('3  WHERE μ IS 0: A LAW. By the formula, μ ≡ 0 exactly when Q ≡ −2 (mod ℓ) and (−m|Q) = −1. Counted:')
for l in (5, 7, 13):
    pl = p % l; mem = zer = 0; qs = []; other_all_zero = []
    for Q in PRIMES:
        if Q == l: continue
        Ns, m, Nst = MEM[Q]; radu = leg(-m, l) == -1
        if Q % l == l - 2:
            cell = radu & (leg(-m, Q) == -1); mem += int(cell.sum()); zer += int((pl[Ns[cell]] == 0).sum()); qs.append(Q)
        for E1 in (1, -1):
            cell = radu & (leg(-m, Q) == E1)
            if cell.sum() >= 12 and not (Q % l == l - 2 and E1 == -1) and bool(np.all(pl[Ns[cell]] == 0)): other_all_zero.append((Q, E1))
    say('   ℓ = %2d:  the %d primes Q ≡ −2 up to %d (%s, …): %d members, %d of them multiples of %d;  any other class all multiples: %s' % (
        l, len(qs), QMAX, ', '.join(map(str, qs[:6])), mem, zer, l, other_all_zero if other_all_zero else 'none'))
say('   Where Q ≡ −1 (mod ℓ), the formula\'s first part, (3|Q)·(1 + Q)·Q⁻², is 0. In the Hecke rule that part alone links p at Q³·m\' to')
say('   p at Q·m\', so p(N) should be a multiple whenever Q³ divides 24N − 1 exactly (N in Radu\'s class; laws of the kind in Ahlgren,')
say('   Allen and Tang\'s Theorem 1.1). DERIVED from the formula, then counted:')
for l in (5, 7, 13):
    pl = p % l; mem = zer = 0; qs = []
    for Q in PRIMES:
        if Q % l != l - 1 or Q ** 3 > 24 * N: continue
        r = pow(24, -1, Q ** 3); Ns = np.arange(r, N + 1, Q ** 3, dtype=np.int64); m3 = (24 * Ns - 1) // Q ** 3
        cell = (m3 % Q != 0) & (leg(-(24 * Ns - 1), l) == -1); mem += int(cell.sum()); zer += int((pl[Ns[cell]] == 0).sum()); qs.append(Q)
    say('   ℓ = %2d:  Q = %s: %d members, %d of them multiples' % (l, ', '.join(map(str, qs)), mem, zer))


def progressions(l, Q):
    A = l * Q ** 3; out = []
    for B in range(A):
        Mv = 24 * B - 1
        if Mv <= 0 or Mv % (Q * Q) or (Mv // (Q * Q)) % Q == 0: continue
        if int(leg(-(Mv // (Q * Q)), Q)) == -1 and int(leg(-Mv, l)) == -1: out.append(B)
    return A, out


say('   the laws as progressions (every member up to %d checked):' % N)
for l, Q in ((7, 5), (5, 13), (13, 11), (5, 23), (7, 19)):
    A, Bs = progressions(l, Q); allN = np.concatenate([np.arange(B, N + 1, A) for B in Bs]); bad = int((p[allN] % l != 0).sum())
    extra = '' if l != 13 else ';  Atkin\'s own example 237 is among them: %s' % (237 in Bs)
    say('   p(%d·n + B) ≡ 0 (mod %d)  [%d = %d·%d³] for %d values of B%s;  members %d, not multiples %d%s' % (
        A, l, A, l, Q, len(Bs), (' (B = %s)' % ', '.join(map(str, Bs))) if len(Bs) <= 12 else ' (B = %s, …)' % ', '.join(map(str, Bs[:6])), len(allN), bad, extra))

say(); say('4  HOW MUCH OF THE EXTRA MULTIPLES THESE LAWS ACCOUNT FOR. "Struck": N in Radu\'s class with a Q ≡ −2 law (Q² exactly, (−m|Q) = −1)')
say('   or a Q ≡ −1 law (Q³ exactly). Prediction: the old laws, plus 1 for struck N, plus chance (1/ℓ) everywhere else; the share of')
say('   struck N from each prime\'s own odds, (Q − 1)/2Q³ and (Q − 1)/Q⁴, multiplied over all primes up to 2,000,000.')
n1 = nn[1:]; M1 = MM[1:]; struck_by = {}
for l in (5, 7, 13):
    st = np.zeros(N + 1, dtype=bool)
    for Q in PRIMES:
        if Q == l: continue
        if Q % l == l - 2:
            Ns, m, _ = MEM[Q]; st[Ns[leg(-m, Q) == -1]] = True
        if Q % l == l - 1 and Q ** 3 <= 24 * N:
            r = pow(24, -1, Q ** 3); Ns = np.arange(r, N + 1, Q ** 3, dtype=np.int64); st[Ns[((24 * Ns - 1) // Q ** 3) % Q != 0]] = True
    struck_by[l] = st
    cls = leg(-MM, l); cls[0] = 99; radu = cls == -1; st &= radu; z = (p % l == 0)
    keep = 1.0
    for Q in sp.primerange(5, 2_000_000):
        if Q % l == l - 2: keep *= 1 - (Q - 1) / (2 * Q ** 3)
        elif Q % l == l - 1: keep *= 1 - (Q - 1) / Q ** 4
    pk = 1 - keep; rshare = z[cls == 0].mean()
    pred = (1 / l) * rshare + (l - 1) / (2 * l) * (1 / l) + (l - 1) / (2 * l) * (pk + (1 - pk) / l)
    sh = lambda msk: '%.5f of %d' % (z[msk].mean(), int(msk.sum()))
    say('   ℓ = %2d:  ℓ | 24N − 1: %s;  the other class: %s;  Radu\'s class, struck: %s;  not struck: %s  (1/ℓ = %.5f)' % (
        l, sh(cls == 0), sh(cls == 1), sh(st), sh(radu & ~st), 1 / l))
    say('           struck share of Radu\'s class: counted %.5f, from the odds %.5f;  all N: counted %.5f, predicted %.5f' % (
        st.sum() / radu.sum(), pk, z[1:].mean(), pred))
    if l in (5, 7):
        excess = z[1:].mean() - (1 / l + (l - 1) / l / l); exp_part = (l - 1) / (2 * l) * pk * (1 - 1 / l)
        say('           extra multiples over the old laws plus chance: %.5f, of which the laws of 3 account for %.5f (%.0f%%)' % (excess, exp_part, 100 * exp_part / excess))

say(); say('5  WHAT IS LEFT OVER')
sqf = np.ones(N + 1, dtype=bool)
for Q in PRIMES:
    r = pow(24, -1, Q * Q); sqf[np.arange(r if r else Q * Q, N + 1, Q * Q)] = False
sqf[0] = False
for l in (5, 7):
    cls = leg(-MM, l); cls[0] = 99; z = (p % l == 0); base = (cls == -1) & ~struck_by[l]; core = base & sqf
    zc = (z[core].mean() - 1 / l) / math.sqrt((1 / l) * (1 - 1 / l) / core.sum())
    shares = []
    for key, mod in ((nn, 8), (nn, 9), (MM, l * l)) + tuple((MM, q) for q in sp.primerange(5, 32) if q != l):
        r = key % mod
        for a in range(mod):
            sel = core & (r == a)
            if sel.sum() > 5000: shares.append(z[sel].mean())
    bands = [(1, 100_000), (100_000, 1_000_000), (1_000_000, N + 1)]
    vals = np.bincount(p[core] % l, minlength=l) / core.sum()
    say('   ℓ = %d, Radu\'s class, not struck, 24N − 1 with no square factor at all: %.5f of %d multiples (z = %+.1f against 1/%d);' % (l, z[core].mean(), int(core.sum()), zc, l))
    say('      with a square factor (not struck): %.5f of %d;  the core share by size: %s' % (
        z[base & ~sqf].mean(), int((base & ~sqf).sum()), ', '.join('%.4f' % z[core & (nn >= lo) & (nn < hi)].mean() for lo, hi in bands)))
    say('      the core share in each of %d classes (N mod 8, N mod 9, 24N − 1 mod %d, 24N − 1 mod each prime 5 to 31): from %.4f to %.4f' % (len(shares), l * l, min(shares), max(shares)))
    say('      the core values on lane %d: %s' % (l, ', '.join('%d: %.4f' % (v, vals[v]) for v in range(l))))
say('   ℓ = 13, the class 13 | 24N − 1: p(13n + 6) ≡ 11·[qⁿ] ∏(1 − qᵏ)¹¹ (mod 13), checked for every n ≤ %d:' % ((N - 6) // 13))
K = (N - 6) // 13; f = np.zeros(K + 1, dtype=np.int64); f[0] = 1
eta = {}
for kk in range(-1000, 1001):
    e = kk * (3 * kk - 1) // 2
    if 0 <= e <= K: eta[e] = eta.get(e, 0) + (-1) ** (kk % 2)
eta3 = {kk * (kk + 1) // 2: (-1) ** kk * (2 * kk + 1) for kk in range(0, 1000) if kk * (kk + 1) // 2 <= K}
for terms in (eta, eta, eta3, eta3, eta3):
    g = np.zeros(K + 1, dtype=np.int64)
    for e, cf in terms.items(): g[e:] = (g[e:] + cf * f[:K + 1 - e]) % 13
    f = g
same = bool(((p[6::13][:K + 1] - 11 * f) % 13 == 0).all())
say('      %s;  multiples of 13 among them: %.5f of %d (1/13 = %.5f)' % (same, (f == 0).mean(), K + 1, 1 / 13))
for l in (17, 19, 23):
    sel = (MM % l == 0); sel[0] = False
    say('   ℓ = %d, the class %d | 24N − 1: multiples %.5f of %d (1/%d = %.5f)' % (l, l, (p[sel] % l == 0).mean(), int(sel.sum()), l, 1 / l))
say(); say('seconds: %.0f' % (time.time() - t0))
