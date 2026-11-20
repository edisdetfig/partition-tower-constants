# A dig, 13 Nov 2026: THE BRAIDED FORMS TAKEN APART, WITH NUMBERS FROM THEIR OWN FIELD (README "The braided forms taken
# apart"; predicted there first). A conjugate group of forms (a pair at 47, 59, 61; a triple at 73, 97) has its constants
# in the unramified extension of degree e of the ℓ-adic numbers. Steps with ℓ-adic numbers cancel a whole group or none of
# it; steps with numbers from that field cancel all forms but one.
#   O = (Z/ℓ^K)[x]/(A(x)), A the group's factor of the constants' polynomial (partitions_flint_family_steps_ℓ.txt, Hensel),
#   so x is the constant α₁ of one member f₁, and the other roots of A (Newton from x^ℓ) are its conjugates.
#   λ₁(5): the root of λ(5)'s polynomial paired with x. On lane ℓ, T₅ = R(T_ℓ) on the base class, so
#   λ₁ ≡ (3|5)·R((3|ℓ)·x); Newton's method lifts it in O.
#   h(y) = (λ(5)'s polynomial)/(y − λ₁);  G₁(M) = Σ κ_j·c(5^(2j)·M), κ from h(T) along the chain (κ in O, by (M|5)).
# Then G₁ = h(λ₁)·c₁ is f₁'s tower alone, and G₁(ℓ²M)/G₁(M) should agree with x in at least (b + 1)/2 digits at levels
# b → b + 2, and with x's conjugates in none. Groups of one form (61's, 97's) run through the same code as a control.
# FLINT's exact partition numbers, at idle priority. COUNTED (the agreements), DERIVED (κ, λ₁).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_braid_split.py [--keep]
import os, sys, re, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly, nmod_mat
import forms6_flint as F
flint.ctx.threads = 4
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_braid_split.txt'), 'x', encoding='utf-8') if KEEP else None
NMAX = 10 ** 12 if "--quick" in sys.argv else 12 * 10 ** 15      # --quick: a short test of the code


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


PC = {}


def pn(N):
    if N not in PC: PC[N] = int(flint.fmpz(N).partitions_p())
    return PC[N]


class Ring:
    """O = (Z/m)[x]/(A), A monic of degree e (coefficients lowest first); elements are lists of e integers."""
    def __init__(self, A, m, L):
        self.A = [a % m for a in A]; self.e = len(A) - 1; self.m = m; self.L = L
    def c(self, k): return [k % self.m] + [0] * (self.e - 1)
    def x(self): return ([0, 1] + [0] * (self.e - 2)) if self.e > 1 else [(-self.A[0]) % self.m]
    def add(self, a, b): return [(u + v) % self.m for u, v in zip(a, b)]
    def sub(self, a, b): return [(u - v) % self.m for u, v in zip(a, b)]
    def scal(self, k, a): return [(k * u) % self.m for u in a]
    def mul(self, a, b):
        e, m, A = self.e, self.m, self.A; pr = [0] * (2 * e - 1)
        for i, u in enumerate(a):
            if u:
                for j, v in enumerate(b): pr[i + j] += u * v
        for t in range(2 * e - 2, e - 1, -1):
            cc = pr[t] % m
            if cc:
                for i in range(e): pr[t - e + i] -= cc * A[i]
            pr[t] = 0
        return [pr[i] % m for i in range(e)]
    def pw(self, a, n):
        r = self.c(1)
        while n:
            if n & 1: r = self.mul(r, a)
            a = self.mul(a, a); n >>= 1
        return r
    def unit(self, a): return any(u % self.L for u in a)
    def inv(self, a):
        e, m, L = self.e, self.m, self.L; cols = []; t = a
        for j in range(e):
            cols.append(t); t = self.mul(t, self.x()) if e > 1 else t
        Mx = [[cols[j][i] for j in range(e)] + [1 if i == 0 else 0] for i in range(e)]
        for cc in range(e):
            piv = next(i for i in range(cc, e) if Mx[i][cc] % L)
            Mx[cc], Mx[piv] = Mx[piv], Mx[cc]; iv = pow(Mx[cc][cc], -1, m); Mx[cc] = [u * iv % m for u in Mx[cc]]
            for i in range(e):
                if i != cc and Mx[i][cc] % m:
                    f = Mx[i][cc]; Mx[i] = [(u - f * v) % m for u, v in zip(Mx[i], Mx[cc])]
        return [Mx[i][e] for i in range(e)]
    def ev(self, poly, y):
        """poly: integer coefficients (lowest first) or ring elements; Horner at y."""
        r = self.c(0)
        for co in reversed(poly): r = self.add(self.mul(r, y), co if isinstance(co, list) else self.c(co))
        return r
    def val(self, a):
        v = min((self.vint(u) for u in a), default=0)
        return v
    def vint(self, u):
        u %= self.m
        if u == 0:
            k = 0; mm = self.m
            while mm > 1: mm //= self.L; k += 1
            return k
        k = 0
        while u % self.L == 0: u //= self.L; k += 1
        return k


def deriv(poly): return [i * poly[i] for i in range(1, len(poly))]


def solve_R(TL, T5, d, L):
    """r with Σ_{i<d} r_i·TL^i = T5 on lane ℓ (Gaussian elimination over the field with ℓ elements)."""
    P = [F.identity(d, L)]
    for i in range(1, d): P.append(P[-1] * TL)
    rows = [[int(P[i][a, b]) for i in range(d)] + [int(T5[a, b])] for a in range(d) for b in range(d)]
    rank = 0; piv = []
    for cc in range(d):
        k = next((i for i in range(rank, len(rows)) if rows[i][cc] % L), None)
        if k is None: return None
        rows[rank], rows[k] = rows[k], rows[rank]; iv = pow(rows[rank][cc], -1, L)
        rows[rank] = [u * iv % L for u in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][cc] % L:
                f = rows[i][cc]; rows[i] = [(u - f * v) % L for u, v in zip(rows[i], rows[rank])]
        piv.append(cc); rank += 1
    if any(r[d] % L for r in rows[rank:]): return None
    return [rows[i][d] for i in range(d)]


t0 = time.time()
for L in (47, 59, 61, 73, 97):
    txt = open(os.path.join(HERE, 'partitions_flint_family_steps_%d.txt' % L), encoding='utf-8').read()
    mA = re.search(r"the constants' polynomial \(lowest first\) on lane ℓ: \[[^\]]*\];  on lane ℓ\^(\d+): \[([^\]]*)\]", txt)
    mL = re.search(r"their polynomial \(lowest first\) on lane ℓ: \[[^\]]*\];  on lane ℓ\^(\d+): \[([^\]]*)\]", txt)
    K = int(mA.group(1)); m = L ** K
    CA = [int(u) for u in mA.group(2).split(',')]; CL = [int(u) for u in mL.group(2).split(',')]; d = len(CA) - 1
    # the base class on lane ℓ, and R with T₅ = R(T_ℓ)
    k = L - 3; n0 = k + 1
    Kz, dd = F.class_basis_exact(k, leg(2, L), leg(3, L)); assert dd == d
    Xc, Rc = F.normalized_basis(nmod_mat(n0, d, [int(Kz[i, j]) % L for i in range(n0) for j in range(d)], L), d, L)
    Hm = F.hecke_matrices(k, L, [('T', 5), ('T', L)])
    T5 = F.restrict(Hm[('T', 5)], Xc, Rc, d); TL = F.restrict(Hm[('T', L)], Xc, Rc, d)
    Rr = solve_R(TL, T5, d, L)
    say(); say('ℓ = %d: %d forms; constants and λ(5) on lane ℓ^%d; T₅ = R(T_ℓ) on lane ℓ with R = %s' % (L, d, K, Rr))
    # the groups: the factors of the constants' polynomial, lifted to lane ℓ^K
    fac = nmod_poly([u % L for u in CA], L).factor()[1]; rest = [u % m for u in CA]; groups = []
    for fp, ex in fac:
        assert ex == 1
        if len(fac) == 1: groups.append(rest); break
        hb = nmod_poly([u % L for u in rest], L) // fp
        Gp, Hh = F.hensel_factor(rest, [int(u) for u in fp.coeffs()], [int(u) for u in hb.coeffs()], L, K)
        groups.append(Gp); rest = Hh
        if len(groups) == len(fac) - 1: groups.append(rest); break
    for A in sorted(groups, key=len, reverse=True):
        O = Ring(A, m, L); e = O.e; x = O.x()
        lam = O.scal(leg(3, 5), O.ev([u % m for u in Rr], O.scal(leg(3, L), x)))      # λ₁ on lane ℓ, from the pairing
        paired = O.val(O.ev(CL, lam)) >= 1
        dCL = deriv(CL)
        for _ in range(K + 1): lam = O.sub(lam, O.mul(O.ev(CL, lam), O.inv(O.ev(dCL, lam))))
        assert O.val(O.ev(CL, lam)) >= K
        # h = CL/(y − λ₁)
        q = [None] * d; q[d - 1] = O.c(CL[d])
        for i in range(d - 1, 0, -1): q[i - 1] = O.add(O.c(CL[i]), O.mul(lam, q[i]))
        assert O.val(O.add(O.c(CL[0]), O.mul(lam, q[0]))) >= K
        KAP = {}
        for e5 in (1, -1):
            kap = e5 * leg(-3, 5) * pow(5, -2, m) % m; q3 = pow(5, -3, m); out = []
            for j in range(d):
                y = [1 if jj == j else 0 for jj in range(d)]; s = []
                for i in range(d):
                    s.append(y[0] % m)
                    if i < d - 1: y = [(y[1] + kap * y[0]) % m] + [(y[t + 1] + q3 * y[t - 1]) % m for t in range(1, len(y) - 1)]
                kj = O.c(0)
                for i in range(d): kj = O.add(kj, O.scal(s[i], q[i]))
                out.append(kj)
            KAP[e5] = out
        conj = []
        if e > 1:
            z = x
            for _ in range(e - 1):
                z = O.pw(z, L)
                for _ in range(K + 1): z = O.sub(z, O.mul(O.ev(A, z), O.inv(O.ev(deriv(A), z))))
                conj.append(z)
        say('   a group of %d form%s (constants\' factor on lane ℓ %s): λ₁(5) paired on lane ℓ: %s;  its conjugates: %d' % (
            e, '' if e == 1 else 's', [u % L for u in A], 'yes' if paired else 'NO', len(conj)))

        def G(M):
            Kj = KAP[leg(M, 5)]; r = O.c(0)
            for j in range(d): r = O.add(r, O.scal(pn((25 ** j * M + 1) // 24) % m, Kj[j]))
            return r

        for b in (1, 3, 5, 7):
            kk = (-pow(L ** b, -1, 24)) % 24; t = time.time(); own = []; oth = []; skipped = 0; n = 0; used = 0
            while used < 6:
                M = L ** b * (24 * n + kk); n += 1
                if M % 5 == 0: continue
                if (25 ** (d - 1) * L * L * M + 1) // 24 > NMAX: break
                used += 1; g0 = G(M)
                if not O.unit(g0): skipped += 1; continue
                ratio = O.mul(G(L * L * M), O.inv(g0))
                own.append(min(O.val(O.sub(ratio, x)), K)); oth.append([min(O.val(O.sub(ratio, z)), K) for z in conj])
            if not used: say('      levels %d → %d: past N = %.1e, not computed' % (b, b + 2, NMAX)); break
            say('      levels %d → %d (%d values of M): digits agreeing with its own constant: %s;  with its conjugates: %s%s  (%.0f s)' % (
                b, b + 2, used, own, oth if conj else '(none: one form)', ('; %d skipped (G(M) ≡ 0 on lane ℓ)' % skipped) if skipped else '',
                time.time() - t))
say(); say('the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
