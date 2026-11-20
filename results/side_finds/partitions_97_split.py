# A dig, 13 Nov 2026: 97'S TOWER IS TWO TOWERS: THE SIMPLE PART SPLIT OFF (README "97's tower is two towers"; predicted
# there first). The class of 97 (weight 94, level 6) has, on lane 97, one form g of its own and three conjugate forms.
# Statement B: each step along a prime Q acts on c(M) = p((M + 1)/24) (97 | M) form by form, with the forms' numbers
# λ_f(Q) = (3|Q)·a_Q(f): c_f(Q²M) = (λ_f + εδ)·c_f(M) (Q ∤ M, ε = (M|Q), δ = −(−3|Q)·Q⁻²) and
# c_f(Q^(2j+2)M) = λ_f·c_f(Q^(2j)M) − Q⁻³·c_f(Q^(2j−2)M). So h(T) at Q′ = 5, with h the cubic factor of the polynomial of
# λ(5) on the class, cancels the three forms:
#     G(M) = Σ_j κ_j·c(5^(2j)·M), j = 0 … 3,  κ from expanding h(T) (it depends on (M|5)),  G(M) = h(λ_g(5))·c_g(M).
# THE FORMS (forms6_flint.py, FLINT only): the class basis's q-expansions modulo 97 (S_j = F^k·P_j(w), as in
# partitions_flint_lucky_dual.py), T_Q on the class from the coefficients S_j(Q·n_i); g from the kernel of A − a, A a
# fixed combination of T_5, T_7, T_11, T_13 and a the root of its linear factor.
# THE CHECKS on exact partition numbers (FLINT): along 97, G(97²M) ≡ 28·G(M); along Q = 7 … 23, one-term laws; the simple
# part's lucky primes below 5000 and their zeros. COUNTED (the checks), DERIVED (κ, λ_g, the expected numbers).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_97_split.py [--keep]
import os, sys, time, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly, nmod_mat
import forms6_flint as F
flint.ctx.threads = 4
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_97_split.txt'), 'x', encoding='utf-8') if KEEP else None
L = 97; k = L - 3; n0 = k + 1; QP = 5; QSCAN = 5000


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


t0 = time.time()
# ---- 1. the class and its q-expansions modulo 97 ----
Kz, d = F.class_basis_exact(k, leg(2, L), leg(3, L))
X0 = nmod_mat(n0, d, [int(Kz[i, j]) % L for i in range(n0) for j in range(d)], L)
assert X0.rank() == d
Xc, Rc = F.normalized_basis(X0, d, L)
Nl = max(QSCAN, L) * max(Rc) + 1; mul = lambda a, b: a.mul_low(b, Nl)
E1, E2, E3, E6 = (F.eta_series(dd, Nl, L) for dd in (1, 2, 3, 6))
E1_5 = E1.pow_trunc(5, Nl)
w = mul(mul(E2, E6.pow_trunc(5, Nl)), mul(E1_5, E3).inverse_series_trunc(Nl)).left_shift(1).truncate(Nl)
Fk = mul(mul(mul(E1_5, E1), E6), mul(E2.pow_trunc(3, Nl), E3.pow_trunc(2, Nl)).inverse_series_trunc(Nl)).pow_trunc(k, Nl)
at = lambda p, i: int(p[i]) if i < p.length() else 0
Bs = [[0] * n0 for _ in range(n0)]; h_ = Fk.truncate(n0); ws = w.truncate(n0)
for mm in range(n0):
    for nn in range(n0): Bs[nn][mm] = at(h_, nn)
    h_ = h_.mul_low(ws, n0)
W = nmod_mat(n0, n0, [x for row in Bs for x in row], L).inv() * Xc
S = []
for j in range(d):
    acc = nmod_poly([int(W[k, j])], L)
    for m in range(k - 1, -1, -1): acc = mul(acc, w) + nmod_poly([int(W[m, j])], L)
    Sj = mul(Fk, acc)
    assert all(at(Sj, nn) == int(Xc[nn, j]) for nn in range(n0))
    S.append([at(Sj, i) for i in range(Nl)])


def TQ(Q):
    pk = pow(Q, k - 1, L)
    return nmod_mat(d, d, [(S[j][Q * Rc[i]] + (pk * S[j][Rc[i] // Q] if Rc[i] % Q == 0 else 0)) % L
                           for i in range(d) for j in range(d)], L)


A = None
for c, q in zip((1, 3, 7, 12), (5, 7, 11, 13)): A = TQ(q) * c if A is None else A + TQ(q) * c
fac = A.charpoly().factor()[1]; pieces = sorted(f.degree() for f, e in fac for _ in range(e))
lin = [f for f, e in fac if f.degree() == 1]
say('ℓ = 97 (weight %d): the class has %d forms; pivot coefficients %s; on lane 97 its pieces have degrees %s  (PREDICTED [1, 3])' % (k, d, Rc, pieces))
a = (-int(lin[0].coeffs()[0]) * pow(int(lin[0].coeffs()[1]), -1, L)) % L
Kv, nul = (A - F.identity(d, L) * a).nullspace(); assert nul == 1
v = [int(Kv[i, 0]) for i in range(d)]; i0 = next(i for i in range(d) if v[i])


def lam_g(Q):
    """(3|Q)·a_Q(g): the eigenvalue of T_Q on the piece of degree 1, times (3|Q); checked to be an eigenvector."""
    T = TQ(Q); u = [sum(int(T[i, j]) * v[j] for j in range(d)) % L for i in range(d)]
    e = u[i0] * pow(v[i0], -1, L) % L
    assert all((u[i] - e * v[i]) % L == 0 for i in range(d)), 'not an eigenvector at Q = %d' % Q
    return leg(3, Q) * e % L if Q != 3 else e


alpha = lam_g(L)
say('   (3|97)·a_97(g) ≡ %d  (PREDICTED 28, the first digit of the 97-adic constant α)  →  %s' % (alpha, 'as predicted' if alpha == 28 else 'NOT as predicted'))
# ---- 2. the cubic factor at Q′ = 5, and the combination κ ----
s5 = leg(3, QP); T5 = TQ(QP) * s5; chi = T5.charpoly(); lg5 = lam_g(QP)
hpol, rem = divmod(chi, nmod_poly([(-lg5) % L, 1], L)); assert rem == 0
hc = [int(x) for x in hpol.coeffs()]
say('   at Q′ = 5: λ_g(5) = %d; the cubic factor h (lowest first) = %s; h(λ_g(5)) = %d (must not be 0)' % (lg5, hc, int(hpol(lg5))))
assert int(hpol(lg5)) != 0


def apply_T(x, kap, q3):
    return [(x[1] + kap * x[0]) % L] + [(x[j + 1] + q3 * x[j - 1]) % L for j in range(1, len(x) - 1)]


def law_at_0(poly, x, kap, q3):
    tot, y = 0, list(x)
    for i, p in enumerate(poly):
        tot = (tot + p * y[0]) % L
        if i < len(poly) - 1: y = apply_T(y, kap, q3)
    return tot


KAPPA = {}
for e5 in (1, -1):
    kap = e5 * leg(-3, QP) * pow(QP, -2, L) % L
    KAPPA[e5] = [law_at_0(hc, [1 if jj == j else 0 for jj in range(4)], kap, pow(QP, -3, L)) for j in range(4)]
say('   κ for (M|5) = +1: %s;  for (M|5) = −1: %s   (G(M) = Σ κ_j·c(5^(2j)·M))' % (KAPPA[1], KAPPA[-1]))
say('   (forms: %.0f s)' % (time.time() - t0))

PC = {}


def c(M):
    """c(M) = p((M + 1)/24) modulo 97."""
    N = (M + 1) // 24
    if N not in PC: PC[N] = int(flint.fmpz(N).partitions_p()) % L
    return PC[N]


def G(M):
    e5 = leg(M, QP); assert e5 != 0 and M % L == 0
    return sum(KAPPA[e5][j] * c(QP ** (2 * j) * M) for j in range(4)) % L


Ms = [L * (24 * n + 23) for n in range(40) if (L * (24 * n + 23)) % QP]
# ---- 3. along 97 ----
t = time.time(); rows = []; ok = nz = 0; ratios = set()
for M in Ms[:20]:
    g0, g2 = G(M), G(L * L * M); ok += g2 == alpha * g0 % L; nz += g0 != 0
    rows.append('%d→%d' % (g0, g2))
    if c(M): ratios.add(c(L * L * M) * pow(c(M), -1, L) % L)
say(); say('Along 97 (levels 1 and 3): G(97²M) ≡ 28·G(M) for %d of %d values of M = 97·(24n + 23) (G(M) not 0 for %d of them)  %s' % (
    ok, len(rows), nz, 'as PREDICTED' if ok == len(rows) else 'NOT as predicted'))
say('   G(M) → G(97²M): %s' % ', '.join(rows))
say('   the whole tower: c(97²M)/c(M) takes %d different values over these M (no single number links them)  (%.0f s)' % (len(ratios), time.time() - t))
# ---- 4. along other primes ----
say(); say('Along other primes: G(Q²M) ≡ (λ_g(Q) + ε·δ)·G(M), ε = (M|Q), δ = −(−3|Q)·Q⁻²')
for Q in (7, 11, 13, 17, 19, 23):
    t = time.time(); lg = lam_g(Q); okq = 0; nzq = 0; mult = {}
    use = [M for M in Ms if M % Q][:10]
    for M in use:
        eps = leg(M, Q); dl = -leg(-3, Q) * pow(Q, -2, L) % L
        g0, g1 = G(M), G(Q * Q * M); okq += g1 == (lg + eps * dl) * g0 % L; nzq += g0 != 0
        if c(M): mult.setdefault(eps, set()).add(c(Q * Q * M) * pow(c(M), -1, L) % L)
    say('   Q = %2d: λ_g = %2d;  %d of %d as predicted (G(M) not 0 for %d);  the whole tower: c(Q²M)/c(M) takes %s different values in the classes ε = +1, −1  (%.0f s)' % (
        Q, lg, okq, len(use), nzq, [len(mult.get(1, ())), len(mult.get(-1, ()))], time.time() - t))
# ---- 5. the simple part's lucky primes ----
say(); say('The simple part\'s lucky primes, every prime Q from 7 to %d (5 and 97 left out)' % QSCAN)
TYPES = (('class ε=+1', 1), ('class ε=−1', -1), ('Ono', 0)); found = {ty: [] for ty, _ in TYPES}; expd = {ty: 0.0 for ty, _ in TYPES}
from sympy import primerange
for Q in primerange(7, QSCAN + 1):
    if Q == L: continue
    lg = lam_g(Q); D = pow(Q, -3, L)
    for ty, eps in TYPES:
        cc = eps * leg(-3, Q) * pow(Q, -2, L) % L if eps else 0
        x = (cc * cc - 4 * D) % L; expd[ty] += (L + (0 if x == 0 else leg(x, L))) / (L * L - 1)
        if lg == cc: found[ty].append(Q)
for ty, _ in TYPES:
    say('   %-11s expected %5.1f (derived)   counted %3d   %s' % (ty, expd[ty], len(found[ty]), found[ty]))
for ty, eps in TYPES:
    for Q in found[ty][:2]:
        t = time.time(); dl = -leg(-3, Q) * pow(Q, -2, L) % L
        if eps:
            mem = [M for M in Ms if M % Q and leg(M, Q) == eps][:2]; ctl = next(M for M in Ms if M % Q and leg(M, Q) == -eps)
            vals = [G(Q * Q * M) for M in mem]; vc = G(Q * Q * ctl); pc = (-2 * eps * dl * G(ctl)) % L
            good = vals == [0, 0] and vc == pc
            say('   %-11s Q = %4d: G(Q²M) for M = %s: %s;  control M′ = %d: G(Q²M′) = %d, predicted −2εδ·G(M′) = %d  %s  (%.0f s)' % (
                ty, Q, mem, vals, ctl, vc, pc, 'as PREDICTED' if good else 'NOT as predicted', time.time() - t))
        else:
            ms = []; m = 1
            while len(ms) < 2:
                if (L * Q ** 3 * m) % 24 == 23 and m % Q and m % QP: ms.append(m)
                m += 1
            vals = [G(L * Q ** 3 * m) for m in ms]
            ctl = next(M for M in Ms if M % Q); vc = G(Q * Q * ctl); pc = (leg(ctl, Q) * dl * G(ctl)) % L
            good = vals == [0, 0] and vc == pc
            say('   %-11s Q = %4d: G(97·Q³·m) for m = %s: %s;  control M′ = %d: G(Q²M′) = %d, predicted ε·δ·G(M′) = %d  %s  (%.0f s)' % (
                ty, Q, ms, vals, ctl, vc, pc, 'as PREDICTED' if good else 'NOT as predicted', time.time() - t))
say(); say('the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
