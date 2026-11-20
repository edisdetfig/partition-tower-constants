# eisenstein_proof.py, 14 Nov 2026 (README "Do the families meet the conditions of the proofs?"; a relayed review:
# "three primes is evidence, not proof, that the 7 system is Eisenstein mod 7").
# A proof, for each seed whose numbers looked Eisenstein (seed_types.py: T_q ≡ 1 + q^(w−1) mod ℓ for q ≤ 97), that its
# mod-ℓ representation really is reducible. Let f be the seed (weight w = ℓ − 1 at level 6, normalized a_1 = 1), and
# θ = q·d/dq. The claim is θ^(ℓ−1)·f ≡ θ^(ℓ−2)·G mod ℓ, where G is the weight-(ℓ + 1) Eisenstein series of level 6
# whose U_2 and U_3 numbers match:
#     a_n(G) = Π σ_ℓ(p^e) over p^e ∥ n with p ∤ 6, times e_2^(v_2(n))·e_3^(v_3(n)), e_p ∈ {1, p^ℓ}.
# θ^(ℓ−1)·f is f with its coefficients at multiples of ℓ removed; the seed is ordinary, so these are not 0.
# - Mod ℓ both sides are modular forms on Γ0(6) of weight ℓ² + ℓ − 2. (Katz: θ raises the weight by ℓ + 1, and
#   E_(ℓ−1) ≡ 1 pads the right side by ℓ − 1.)
# - Sturm's bound for Γ0(6) (index 12): two such forms are equal once their coefficients agree for n ≤ ℓ² + ℓ − 2.
# - Then a_q(f) ≡ q^(ℓ−2)·(1 + q^ℓ) ≡ 1 + q^(ℓ−2) for every prime q ∤ 6ℓ. So ρ̄_f ≅ 1 ⊕ ω^(ℓ−2) (Chebotarev,
#   Brauer–Nesbitt), and the family's ρ̄ = ρ̄_f ⊗ ω^m ≅ ω^m ⊕ ω^(m−1).
# For contrast the same comparison is run on every other seed, which must fail.
# COUNTED (the coefficients), DERIVED (the conclusion, from the standard theorems named).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/eisenstein_proof.py ℓ [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat, nmod_poly
import forms6_flint as F

flint.ctx.threads = 2
L = int([a for a in sys.argv[1:] if a.isdigit()][0]); KEEP = '--keep' in sys.argv
QS = [q for q in (5, 7, 11, 13) if q != L]; LM = L - 1; NPW = 10; LEN = L * L + L - 1      # n = 0 … ℓ² + ℓ − 2
OUT = os.path.join(HERE, 'eisenstein_proof_%d.txt' % L)
if KEEP and os.path.exists(OUT): sys.exit('stopped before computing anything: %s exists' % OUT)
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def sign(p, k): return leg(p, L) ** (1 + (k - (L - 3)) // (L - 1))


def split_pieces(Mq, d, Nn):
    m1 = L ** Nn; pieces = [((), F.identity(d, m1), list(range(d)), d)]
    for q in QS:
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


def class_vectors(w, s2, s3, N):
    """The class at weight w (new at 2 and 3 with these signs exactly), split by system: for each 1-dimensional piece,
    its label and its vector mod ℓ in forms6_flint's coordinates (the coefficients of q^0 … q^w)."""
    n = w + 1; m = L ** N
    H = F.hecke_matrices(w, m, [('T', q) for q in QS] + [('U', 2), ('U', 3)])
    c2 = s2 * 2 ** ((w - 2) // 2) % m; c3 = s3 * 3 ** ((w - 2) // 2) % m
    Y, R, r = F.kernel_part(H[('U', 2)] - F.identity(n, m) * c2, L, N)
    if r == 0: return []
    Zm, Rz, r2 = F.kernel_part(F.restrict(H[('U', 3)], Y, R, r) - F.identity(r, m) * c3, L, N)
    if r2 == 0: return []
    Y, R = Y * Zm, [R[j] for j in Rz]; r = r2
    U2W = F.restrict(H[('U', 2)], Y, R, r); U3W = F.restrict(H[('U', 3)], Y, R, r)
    Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)] + \
         [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(r)] for i in range(r)]
    ker, N1 = F.padic_kernel(Bm, 2 * r, r, L, N); d = len(ker)
    if d == 0: return []
    m1 = L ** N1
    Xm = F.mat_mod(r, d, [ker[c][i] % m1 for i in range(r) for c in range(d)], m1)
    Rn = F.pivot_rows(nmod_mat(r, d, [ker[c][i] % L for i in range(r) for c in range(d)], L), L)
    XR = F.inverse_mod(F.mat_mod(d, d, [ker[b][Rn[a]] % m1 for a in range(d) for b in range(d)], m1), L, N1)
    def on_new(T):
        TX = F.mat_mod(r, r, [int(v) % m1 for v in T.entries()], m1) * Xm
        return XR * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m1)
    out = []
    Ym = F.mat_mod(n, r, [int(v) % m1 for v in Y.entries()], m1)
    for lab, Yp, Rp, rp in split_pieces({q: on_new(F.restrict(H[('T', q)], Y, R, r)) for q in QS}, d, N1):
        if rp != 1: continue
        v = Ym * Xm * Yp                                         # n × 1, the piece's form in full coordinates
        out.append((lab, [int(v[i, 0]) % L for i in range(n)]))
    return out


def basis_series(w):
    """forms6_flint's basis H_m = q^m·φ·u^(w−m) of M_w(Γ0(6)), mod ℓ, to length LEN."""
    E1, E2, E3, E6 = (F.eta_series(d, LEN, L) for d in (1, 2, 3, 6))
    mul = lambda a, b: a.mul_low(b, LEN)
    u = mul(mul(E1.pow_trunc(5, LEN), E3), mul(E2, E6.pow_trunc(5, LEN)).inverse_series_trunc(LEN))
    g1 = mul(mul(E1, E6.pow_trunc(6, LEN)), mul(E2.pow_trunc(2, LEN), E3.pow_trunc(3, LEN)).inverse_series_trunc(LEN))
    psi = g1.pow_trunc(w, LEN); out = {}
    for i in range(w + 1):
        mm = w - i
        c = [int(x) for x in psi.coeffs()[:LEN]]; c += [0] * (LEN - len(c))
        out[mm] = [0] * mm + c[:LEN - mm]
        if i < w: psi = mul(psi, u)
    return out


def g_coeff(nn, e2, e3):
    """a_n of the weight-(ℓ+1) Eisenstein series with U_2 = e2, U_3 = e3 (normalized a_1 = 1), mod ℓ."""
    v = 1; x = nn
    for p, e in ((2, e2), (3, e3)):
        while x % p == 0: x //= p; v = v * e % L
    p = 5
    while x > 1:
        if p * p > x: p = x
        if x % p == 0:
            s = 1; t = 1
            while x % p == 0: x //= p; t = t * pow(p, L, L) % L; s = (s + t) % L
            v = v * s % L
        p += 1 if p == 2 else 2
    return v


c2 = sign(2, L - 3) * pow(2, (L - 5) // 2, L) % L; c3 = sign(3, L - 3) * pow(3, (L - 5) // 2, L) % L
say('ℓ = %d: is θ^(ℓ−1)·f ≡ θ^(ℓ−2)·G mod ℓ for each seed f, G the weight-(ℓ+1) Eisenstein series of level 6 with '
    'matching U_2, U_3? Coefficients compared for n = 0 … ℓ² + ℓ − 2 = %d (Sturm\'s bound at that weight on Γ0(6))' % (
        L, LEN - 1))
t0 = time.time(); proven = 0; tested = 0
for r in range(2, L, 2):
    w = r + 2; Bser = None
    for m in [m for m in range(LM) if (r + 2 * m - (L - 5)) % LM == 0]:
        u2 = pow(2, -m, L) * c2 % L; u3 = pow(3, -m, L) * c3 % L
        s2 = 1 if u2 * pow(2, -(r // 2), L) % L == 1 else -1; s3 = 1 if u3 * pow(3, -(r // 2), L) % L == 1 else -1
        for lab, v in class_vectors(w, s2, s3, NPW):
            if Bser is None: Bser = basis_series(w)
            # the coordinates c with Σ c_m H_m = the form: B is lower unitriangular (H_m starts with q^m)
            c = [0] * (w + 1)
            for mm in range(w + 1):
                c[mm] = (v[mm] - sum(c[j] * Bser[j][mm] for j in range(mm))) % L
            f = [sum(c[j] * Bser[j][nn] for j in range(w + 1)) % L for nn in range(LEN)]
            a1 = f[1]
            if a1 == 0: say('   %r (weight %d·ω^%d): a_1 ≡ 0, skipped' % (lab, w, m)); continue
            inv = pow(a1, -1, L); f = [x * inv % L for x in f]
            e2 = f[2] * pow(pow(2, L - 2, L), -1, L) % L; e3 = f[3] * pow(pow(3, L - 2, L), -1, L) % L
            fits = e2 in (1, 2 % L) and e3 in (1, 3 % L)
            first_bad = None
            if fits:
                for nn in range(LEN):
                    if nn % L == 0: continue                     # θ^(ℓ−1) kills these on both sides (and n = 0)
                    if pow(nn, L - 2, L) * g_coeff(nn, e2, e3) % L != f[nn]: first_bad = nn; break
            ok = fits and first_bad is None; tested += 1; proven += ok
            fam = tuple(x * pow(q, m, L) % L if isinstance(x, int) else x for q, x in zip(QS, lab))
            say('   family %r, seed weight %d·ω^%d: %s' % (
                fam, w, m, 'θ^(ℓ−1)·f ≡ θ^(ℓ−2)·G with (U_2, U_3) of G = (%d, %d): all %d coefficients agree. PROVEN: ρ̄ ≅ ω^%d ⊕ ω^%d'
                % (e2, e3, LEN, m % LM, (m - 1) % LM) if ok else 'not Eisenstein of this kind (%s)' % (
                    'its a_2, a_3 fit no such G' if not fits else 'first differs at n = %d' % first_bad)))
say(); say('proven Eisenstein: %d of %d seeds  (%.0f s)' % (proven, tested, time.time() - t0))
