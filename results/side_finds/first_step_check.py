# first_step_check.py, 14 Nov 2026 (README "Why exactly 13, 17, 19 and 23 settle fast"; a relayed review: "The slope-1
# side relies on Theorem 1.3, which needs ρ̄ irreducible ... check irreducibility for the first-step families").
# For each prime 29 ≤ ℓ ≤ LMAX: the first-step forms, the newforms g of level 6 and weight w = ℓ − 5 with
# U_2 = (2|ℓ)·2^((w−2)/2) and U_3 = (3|ℓ)·3^((w−2)/2). Twisted by ω, these are the class's forms added at the first weight
# step (README). Each is split into its systems on lane ℓ (T_q for every prime 5 ≤ q ≤ 43 other than ℓ), and for each
# system:
#   ordinary   : T_ℓ invertible mod ℓ on it. T_ℓ ≡ U_ℓ mod ℓ, read from the coefficients a_(ℓp) at n ≤ w positions p. For
#                these weights (≤ ℓ + 1), ordinary means ρ̄ is reducible at ℓ, as Theorem 1.3 needs. Its Serre weight is
#                then Sym^(ℓ−7) ⊗ det, with a = ℓ − 7, inside 2 … ℓ − 5 for ℓ ≥ 11 ("very generic").
#   irreducible: if some T_q has an eigenvalue outside F_ℓ, ρ̄ is not Eisenstein. Otherwise, for each of the ℓ − 1 pairs
#                ρ̄ ≅ ω^a ⊕ ω^b with a + b ≡ w − 1 (mod ℓ − 1), a prime q with a_q ≢ q^a + q^b (mod ℓ). That is a proof.
# A prime is covered when at least one of its first-step systems is ordinary and irreducible.
# COUNTED (the numbers), DERIVED (the verdicts).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/first_step_check.py [LMAX] [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat, nmod_poly
import forms6_flint as F

flint.ctx.threads = 2
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; LMAX = nums[0] if nums else 613; KEEP = '--keep' in sys.argv
LMIN = nums[1] if len(nums) > 1 else 29
NP = 8
OUT = os.path.join(HERE, 'first_step_check_%d-%d.txt' % (LMIN, LMAX))
if KEEP and os.path.exists(OUT): sys.exit('stopped before computing anything: %s exists' % OUT)
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def is_prime(n): return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def split_pieces(Mq, d, Nn, qs, L):
    m1 = L ** Nn; pieces = [((), F.identity(d, m1), list(range(d)), d)]
    for q in qs:
        newp = []
        for lab, Yp, Rp, rp in pieces:
            Mqr = F.restrict(Mq[q], Yp, Rp, rp)
            for f, ee in F.reduce_mat(Mqr, L).charpoly().factor()[1]:
                co = [int(c) for c in f.coeffs()]
                Bq = Mqr - F.identity(rp, m1) * ((-co[0]) % L) if len(co) == 2 else F.mat_poly(co, Mqr)
                Zm, Rz, r2 = F.kernel_part(Bq, L, Nn)
                if r2: newp.append((lab + (((-co[0]) % L) if len(co) == 2 else tuple(co),), Yp * Zm, [Rp[j] for j in Rz], r2))
        pieces = newp
    return pieces


def check(L):
    NP = max(3, min(8, int(61 / L.bit_length())))            # keep ℓ^N below 2^62 (FLINT's word-size matrices)
    w = L - 5; n = w + 1; m = L ** NP; s2, s3 = leg(2, L), leg(3, L)
    QX = [q for q in range(5, 44) if is_prime(q) and q != L]
    H = F.hecke_matrices(w, m, [('T', q) for q in QX] + [('U', 2), ('U', 3)])
    c2 = s2 * 2 ** ((w - 2) // 2) % m; c3 = s3 * 3 ** ((w - 2) // 2) % m
    Y, R, r = F.kernel_part(H[('U', 2)] - F.identity(n, m) * c2, L, NP)
    if r == 0: return 0, []
    Zm, Rz, r2 = F.kernel_part(F.restrict(H[('U', 3)], Y, R, r) - F.identity(r, m) * c3, L, NP)
    if r2 == 0: return 0, []
    Y, R = Y * Zm, [R[j] for j in Rz]; r = r2
    U2W = F.restrict(H[('U', 2)], Y, R, r); U3W = F.restrict(H[('U', 3)], Y, R, r)
    Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)] + \
         [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(r)] for i in range(r)]
    ker, N1 = F.padic_kernel(Bm, 2 * r, r, L, NP); d = len(ker)
    if d == 0: return 0, []
    assert N1 >= 1
    m1 = L ** N1
    Xm = F.mat_mod(r, d, [ker[c][i] % m1 for i in range(r) for c in range(d)], m1)
    Rn = F.pivot_rows(nmod_mat(r, d, [ker[c][i] % L for i in range(r) for c in range(d)], L), L)
    XR = F.inverse_mod(F.mat_mod(d, d, [ker[b][Rn[a]] % m1 for a in range(d) for b in range(d)], m1), L, N1)
    def on_new(T):
        TX = F.mat_mod(r, r, [int(v) % m1 for v in T.entries()], m1) * Xm
        return XR * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m1)
    Mq = {q: on_new(F.restrict(H[('T', q)], Y, R, r)) for q in QX}
    # the class basis (the columns of Xm) in full coordinates (coefficients of q^0 … q^w), mod ℓ
    Ym = F.mat_mod(n, r, [int(v) % m1 for v in Y.entries()], m1)
    V = F.reduce_mat(Ym * Xm, L)                              # n × d
    P = F.pivot_rows(V, L)                                    # d coefficient positions where the basis is invertible
    VP = nmod_mat(d, d, [int(V[P[i], j]) for i in range(d) for j in range(d)], L)
    # T_ℓ on the class mod ℓ: (T_ℓ f)(p) ≡ f(ℓp), so the coefficients a_(ℓp) of each basis form are needed
    Lt = L * max(P) + 1
    E1, E2, E3, E6 = (F.eta_series(dd, Lt, L) for dd in (1, 2, 3, 6))
    mul = lambda a, b: a.mul_low(b, Lt)
    u = mul(mul(E1.pow_trunc(5, Lt), E3), mul(E2, E6.pow_trunc(5, Lt)).inverse_series_trunc(Lt))
    g1 = mul(mul(E1, E6.pow_trunc(6, Lt)), mul(E2.pow_trunc(2, Lt), E3.pow_trunc(3, Lt)).inverse_series_trunc(Lt))
    psi = g1.pow_trunc(w, Lt)
    coefH = {}                                                # H_mm's coefficients at q^0 … q^w and at q^(ℓp)
    want = sorted(set(range(n)) | {L * p for p in P})
    for i in range(w + 1):
        mm = w - i
        coefH[mm] = {t: (int(psi[t - mm]) if t - mm >= 0 else 0) for t in want}
        if i < w: psi = mul(psi, u)
    Tl = [[0] * d for _ in range(d)]
    for j in range(d):
        v = [int(V[t, j]) for t in range(n)]
        c = [0] * n                                           # v = Σ c_mm H_mm (H_mm starts with q^mm)
        for mm in range(n):
            c[mm] = (v[mm] - sum(c[jj] * coefH[jj][mm] for jj in range(mm))) % L
        for i_, p in enumerate(P):
            Tl[i_][j] = sum(c[mm] * coefH[mm][L * p] for mm in range(n)) % L
    Ml = VP.inv() * nmod_mat(d, d, [Tl[i][j] for i in range(d) for j in range(d)], L)   # T_ℓ in the class basis
    out = []
    for lab, Yp, Rp, rp in split_pieces(Mq, d, N1, QX, L):
        Ypl = F.reduce_mat(Yp, L)
        Mlp = F.restrict(Ml, Ypl, Rp, rp)
        ordinary = int(Mlp.det()) % L != 0
        rational = all(isinstance(v, int) for v in lab)
        if not rational: irr, how = True, 'eigenvalues outside F_ℓ'
        else:
            ev = dict(zip(QX, lab)); missing = []
            for a in range(L - 1):
                b = (w - 1 - a) % (L - 1)
                if not any(ev[q] != (pow(q, a, L) + pow(q, b, L)) % L for q in QX): missing.append((a, b))
            irr = not missing
            how = 'every pair ω^a ⊕ ω^b missed at some q ≤ 43' if irr else 'Eisenstein-like for (a, b) = %s on every q ≤ 43' % missing
        out.append((rp, rational, ordinary, irr, how, int(Mlp.det()) % L))
    return d, out


say('ℓ = %d … %d: the first-step forms (level 6, weight ℓ − 5, signs ((2|ℓ), (3|ℓ))), by system: ordinary? irreducible?' % (LMIN, LMAX))
t0 = time.time(); covered = 0; primes = 0; bad = []
for L in range(LMIN, LMAX + 1):
    if not is_prime(L): continue
    t1 = time.time(); d, out = check(L); primes += 1
    good = [o for o in out if o[2] and o[3]]
    covered += bool(good)
    if not good or len(good) < len(out): bad.append(L)
    say('   ℓ = %3d: %d forms (⌊ℓ/24⌋ = %d) in %d systems (sizes %s); ordinary and irreducible: %d system%s%s  (%.0f s)' % (
        L, d, L // 24, len(out), [o[0] for o in out], len(good), '' if len(good) == 1 else 's',
        '' if len(good) == len(out) else '; the others: %s' % [(o[0], 'ordinary' if o[2] else 'NOT ordinary', o[4]) for o in out if not (o[2] and o[3])],
        time.time() - t1))
say(); say('primes covered (at least one first-step system ordinary and irreducible): %d of %d; primes where some first-step '
           'system is not: %s  (%.0f s)' % (covered, primes, bad or 'none', time.time() - t0))
