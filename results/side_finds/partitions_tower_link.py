# A dig, 11 Nov 2026: THE TOWER LINK FOR 5, IN ONE LINE. partitions_tower.txt (section 4) found that two levels below
# Atkin's class p is 3 times p on lane 5, and that the other class gains one more factor of 5. Here the whole link as one
# formula,
#     p(25N − 1) ≡ 25·N·(N − 3)·p(N)   (mod 125)   for every N ≥ 1,
# checked on every N and at every depth (laws divided out); why it holds (two identities from Ramanujan's unpublished
# manuscript, QUOTED from Berndt and Ono's edition, plus a short derivation written out in section 3); why no other prime
# up to 31 has a link at level 0; and what it predicts far beyond the recurrence, checked with exact values of p(N)
# (Hardy–Ramanujan–Rademacher formula, sympy). Checks COUNTED; the derivation DERIVED; identities QUOTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_tower_link.py [--keep]
import os, sys, time, collections
import numpy as np
from sympy.functions.combinatorial.numbers import partition
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_tower_link.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); NMAX = 2_400_000
G, S, k = [], [], 1
while True:
    a, b = k * (3 * k - 1) // 2, k * (3 * k + 1) // 2
    if a > NMAX: break
    s = 1 if k % 2 else -1; G.append(a); S.append(s)
    if b <= NMAX: G.append(b); S.append(s)
    k += 1
o = np.argsort(G); G = np.array(G, dtype=np.int64)[o]; S = np.array(S, dtype=np.int64)[o]
LANES = {'5': 5 ** 21, 'rest': 7 ** 4 * 11 ** 4 * 13 * 17 * 19 * 23 * 29 * 31}
P = {}
for name, L in LANES.items():
    p = np.zeros(NMAX + 1, dtype=np.int64); p[0] = 1; c = 0
    for n in range(1, NMAX + 1):
        while c < len(G) and G[c] <= n: c += 1
        p[n] = int(np.dot(S[:c], p[n - G[:c]])) % L
    P[name] = p
p5, pr = P['5'], P['rest']
say('p(n) for every n up to %d on the lanes 5²¹ and 7⁴·11⁴·13·17·19·23·29·31 (%.0f s); p(100) right: %s' % (
    NMAX, time.time() - t0, int(p5[100]) == 190569292 % LANES['5'] and int(pr[100]) == 190569292 % LANES['rest']))


def vl(x, l):
    v = 0
    while x % l == 0: x //= l; v += 1
    return v


NN = (NMAX + 1) // 25; Ns = np.arange(1, NN + 1, dtype=np.int64)
say(); say('1  THE FORMULA: p(25N − 1) ≡ 25·N·(N − 3)·p(N) (mod 125), for every N from 1 to %d' % NN)
lhs = p5[25 * Ns - 1] % 125; rhs = (25 * (Ns * (Ns - 3) % 5) * (p5[Ns] % 5)) % 125
say('   wrong: %d of %d' % (int((lhs != rhs).sum()), NN))
say('   what it says, class by class (N on lane 5; N·(N − 3) ≡ (24N − 1)² − 1 on lane 5, so the classes are those of 24N − 1):')
for r, what in ((0, 'p(125n − 1): always a multiple of 125'), (3, 'p(125n − 51): always a multiple of 125'),
                (4, 'p(125n − 26): always a multiple of 125 (p(N) is a multiple of 5 here: Ramanujan\'s first law)'),
                (1, 'p(25N − 1)/25 ≡ 3·p(N): Atkin\'s class'), (2, 'p(25N − 1)/25 ≡ 3·p(N): Atkin\'s class')):
    sel = (Ns % 5 == r); ff = (r * (r - 3)) % 5
    if ff == 0 or r == 4: ok = int((lhs[sel] == 0).sum())
    else: ok = int(((p5[25 * Ns[sel] - 1] // 25) % 5 == (3 * p5[Ns[sel]]) % 5).sum())
    say('     N ≡ %d: N(N − 3) ≡ %d;  %s  (%d of %d)' % (r, ff, what, ok, int(sel.sum())))
say('   the three multiples of 125 are in Ramanujan\'s manuscript (Berndt–Ono, eq. 21.3: p(125n − 1), p(125n − 26), p(125n − 51) ≡ 0)')

say(); say('2  AT EVERY DEPTH: p(25N − 1)/5^(v+2) ≡ N·(N − 3)·p(N)/5^v (mod 5), where 5^v is the exact power of 5 in 24N − 1')
cnt, bad = collections.Counter(), collections.Counter()
for N in Ns.tolist():
    v = vl(24 * N - 1, 5); a = int(p5[N]); b = int(p5[25 * N - 1])
    cnt[v] += 1; bad[v] += (b // 5 ** (v + 2) - N * (N - 3) * (a // 5 ** v)) % 5 != 0
say('   ' + ';  '.join('v = %d: %d checks, %d wrong' % (v, cnt[v], bad[v]) for v in sorted(cnt)))
say('   for v ≥ 1, N ≡ 4 and N(N − 3) ≡ 4 = −1: the first digit of the constant of partitions_tower.txt (4, 2, 3, 1, …)')

say(); say('3  WHY (level 0). Two identities and one derivation:')
say('   (a) QUOTED, Ramanujan (Berndt–Ono, eq. 21.4): Σ p(25n + 24)qⁿ = 75·(q;q)²³ + 125·(…), i.e. p(25n + 24) ≡ 75·[qⁿ](q;q)²³ (mod 125)')
J3 = [(k_ * (k_ + 1) // 2, (-1) ** k_ * (2 * k_ + 1)) for k_ in range(0, 2000) if k_ * (k_ + 1) // 2 <= NN]
E1 = [(k_ * (3 * k_ - 1) // 2, (-1) ** k_) for k_ in range(-400, 401) if 0 <= k_ * (3 * k_ - 1) // 2 <= NN]


def times(f, sparse, mod):
    g = np.zeros_like(f)
    for e_, c_ in sparse: g[e_:] = (g[e_:] + c_ * f[:len(f) - e_]) % mod
    return g


e23 = np.zeros(NN + 1, dtype=np.int64); e23[0] = 1
for _ in range(7): e23 = times(e23, J3, 125)                      # (q;q)²¹ = ((q;q)³)⁷, Jacobi's series for (q;q)³
for _ in range(2): e23 = times(e23, E1, 125)                      # times (q;q)² by Euler's pentagonal series
n_ = np.arange(0, NN, dtype=np.int64)
say('       checked for n = 0 to %d: wrong %d' % (NN - 1, int(((p5[25 * n_ + 24] - 75 * e23[n_]) % 125 != 0).sum())))
say('   (b) DERIVED: N·(N − 3)·p(N) ≡ 3·[q^(N−1)](q;q)²³ (mod 5). With D = q·d/dq, 1/E = Σ p(n)qⁿ, P, Q, R = Ramanujan\'s E₂, E₄, E₆:')
say('       D(1/E) = (1 − P)/24 · 1/E and DP = (P² − Q)/12 (classical); his eq. 1.2, Q ≡ 1 and R ≡ P (mod 5); his eq. 1.1, Q³ − R² =')
say('       1728·q·(q;q)²⁴. On lane 5 these give (D² − 3D)(1/E) ≡ (1 − P²)/E ≡ (Q³ − R²)/E ≡ 3·q·(q;q)²⁴/E = 3·q·(q;q)²³, and the')
say('       coefficient of q^N on the left is (N² − 3N)·p(N).')
say('       checked for N = 1 to %d: wrong %d' % (NN, int((((Ns * (Ns - 3)) % 5 * (p5[Ns] % 5) - 3 * e23[Ns - 1]) % 5 != 0).sum())))
say('   (a) and (b) together give the formula of section 1 for every N: p(25N − 1) ≡ 75·[q^(N−1)](q;q)²³ ≡ 25·N·(N − 3)·p(N) (mod 125).')
say('   In the language of 24N − 1 = m: multiplying the coefficient at q^m of 1/η(24z) by m² − 1 gives 3·η(24z)²³ on lane 5, and')
say('   η(24z)²³ is (by (a)) the generating function two levels down. Deeper levels (section 2, v ≥ 1) are counted, not derived here.')

say(); say('4  ONLY ONE DIGIT AT LEVEL 0: p(25N − 1)/25 against p(N) on lane 25 (5 not dividing 24N − 1, p(N) not a multiple of 5)')
for modN in (25, 125):
    grp = collections.defaultdict(set)
    for N in Ns.tolist():
        if (24 * N - 1) % 5 == 0 or int(p5[N]) % 5 == 0: continue
        grp[N % modN].add((int(p5[25 * N - 1]) // 25 % 25) * pow(int(p5[N]) % 25, -1, 25) % 25)
    say('   grouped by N mod %d: %d groups; ratios per group: from %d to %d of the 5 possible (one ratio would mean a law mod 625)' % (
        modN, len(grp), min(len(s_) for s_ in grp.values()), max(len(s_) for s_ in grp.values())))
say('   (from level 1 down, the constant holds to more digits: 3, 6, 8, 11, … in partitions_tower.txt)')

say(); say('5  OTHER PRIMES: is p(N\')/ℓ^law ≡ f·p(N)/ℓ^law on lane ℓ, with f fixed by the residue of N, where 24N\' − 1 = ℓ²(24N − 1)?')
LAW = lambda l, b: b if l in (5, 11) else (b // 2 + 1 if (l == 7 and b >= 1) else 0)
for l in (7, 11, 13, 17, 19, 23, 29, 31):
    sh = (l * l - 1) // 24; g0 = collections.defaultdict(set); r13 = collections.Counter(); n0 = n1 = 0
    for N in range(1, (NMAX + sh) // (l * l) + 1):
        v = vl(24 * N - 1, l)
        if v > 1: continue
        a = int(pr[N]); b = int(pr[l * l * N - sh]); ea, eb = LAW(l, v), LAW(l, v + 2)
        a = (a // l ** ea) % l; b = (b // l ** eb) % l
        if a == 0: continue
        if v == 0: g0[N % l].add(b * pow(a, -1, l) % l); n0 += 1
        else: r13[b * pow(a, -1, l) % l] += 1; n1 += 1
    say('   ℓ = %2d: level 0 → 2: %d pairs, residues of N with one ratio: %d of %d (ratios seen per residue: %d to %d);  level 1 → 3: %d pairs, ratios %s' % (
        l, n0, sum(1 for s_ in g0.values() if len(s_) == 1), len(g0), min(len(s_) for s_ in g0.values()), max(len(s_) for s_ in g0.values()),
        n1, ', '.join('%d (%d times)' % kv for kv in r13.most_common(3))))
say('   For 7, two levels down is again (q;q)²³: QUOTED, Jameson Cor. 2.2, p(49n + 47) ≡ 5·49·[qⁿ](q;q)²³ (mod 7³). Checked:')
e23_7 = np.zeros((NMAX - 47) // 49 + 1, dtype=np.int64); e23_7[0] = 1
J3b = [(e_, c_) for e_, c_ in [(k_ * (k_ + 1) // 2, (-1) ** k_ * (2 * k_ + 1)) for k_ in range(0, 2000)] if e_ < len(e23_7)]
E1b = [(e_, c_) for e_, c_ in [(k_ * (3 * k_ - 1) // 2, (-1) ** k_) for k_ in range(-400, 401)] if 0 <= e_ < len(e23_7)]
for _ in range(7): e23_7 = times(e23_7, J3b, 343)
for _ in range(2): e23_7 = times(e23_7, E1b, 343)
nn7 = np.arange(len(e23_7), dtype=np.int64)
say('   n = 0 to %d: wrong %d.  But on lane 7, (q;q)²³ at q^(N−1) against p(N): ratios seen per residue of N: %s' % (
    len(nn7) - 1, int(((pr[49 * nn7 + 47] - 5 * 49 * e23_7) % 343 != 0).sum()),
    [len({(int(e23_7[N - 1]) % 7) * pow(int(pr[N]) % 7, -1, 7) % 7 for N in range(1, len(e23_7)) if N % 7 == r and int(pr[N]) % 7})
     for r in range(7) if (24 * r - 1) % 7]))
say('   so no factor fixed by N carries p(N) to two levels down on lane 7: the bridge of section 3 (b) exists on lane 5 only.')

say(); say('6  PREDICTIONS FAR BEYOND THE RECURRENCE, checked with exact values of p(N)')
ok = 0; nd = 0
for N in range(10 ** 8, 10 ** 8 + 10):
    big = int(partition(25 * N - 1)); nd = max(nd, len(str(big)))
    ok += (big - 25 * N * (N - 3) * int(partition(N))) % 125 == 0
say('   (a) the formula at N = 10⁸ to 10⁸ + 9 (p(25N − 1) has %d digits): %d of 10 right' % (nd, ok))
say('   (b) k steps down: p(N_k)/5^(v+2k) ≡ (−1)^(k−1)·N·(N − 3)·p(N)/5^v (mod 5), where 24N_k − 1 = 25^k·(24N − 1); also checked:')
say('       5^(v+2k) divides p(N_k) (Watson\'s law) and, when N(N − 3) ≡ 0, so does 5^(v+2k+1)')
rows = []; okk = tot = 0
for N in range(1, 13):
    for k_ in range(1, 8 if N <= 2 else 7):
        Nk = 25 ** k_ * N - (25 ** k_ - 1) // 24; v = vl(24 * N - 1, 5); a = int(partition(N)); b = int(partition(Nk))
        f = (N * (N - 3)) % 5; good = b % 5 ** (v + 2 * k_) == 0 and (b // 5 ** (v + 2 * k_) - (-1) ** (k_ - 1) * f * (a // 5 ** v)) % 5 == 0
        if f == 0: good = good and b % 5 ** (v + 2 * k_ + 1) == 0
        okk += good; tot += 1; rows.append((Nk, good, len(str(b))))
big = max(rows)
say('       N = 1 to 12, k = 1 to 6 (and k = 7 for N = 1, 2): %d of %d right; the largest N_k = %d (p of it has %d digits)' % (
    okk, tot, big[0], big[2]))
say(); say('seconds: %.0f' % (time.time() - t0))
