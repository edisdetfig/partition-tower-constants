# seed_types.py, 14 Nov 2026 (README "Do the families meet the conditions of the proofs?").
# For every family of the class at ℓ, its seed (serre_weights.py's step 2, repeated here) and the facts that the
# theorems' hypotheses ask for:
#   ordinary : the seed's a_ℓ ≢ 0 mod ℓ. For weights 2 … ℓ + 1 that is the same as ρ̄ being reducible at ℓ (a
#              non-ordinary form of such a weight has ρ̄ irreducible at ℓ: Deligne, Fontaine). It is also Bergdall and
#              Pollack's regularity test (their Theorem 2.5) applied to the seed.
#   Eisenstein: whether T_q ≡ q^a + q^b mod ℓ on the seed for every prime q ≤ 97 other than 2, 3, ℓ (ρ̄ globally
#              reducible), with the exponents found.
#   a, b     : the Serre weight Sym^a ⊗ det^b in Liu–Truong–Xiao–Zhao's notation: a = r and b = m for a seed of weight
#              r + 2 twisted by ω^m with r ≤ ℓ − 2; the weight ℓ + 1 seed is their "a = 0, très ramifié" (Sym^(ℓ−1)).
# Then the verdict for each family:
#   LTXZ Theorem 1.3 (as read 14 Nov): ℓ ≥ 11, ρ̄ absolutely irreducible, reducible at ℓ, very generic 2 ≤ a ≤ ℓ − 5;
#   Wang 2026 (local ghost theorem, très ramifié): ℓ ≥ 11 and the Steinberg weight;
#   Bergdall–Pollack Conjecture 3.1 (regular ρ̄, including ω^j ⊕ ω^(j+1)): the seed is ordinary.
# COUNTED (the seeds' numbers), DERIVED (the verdicts, from the hypotheses as quoted in README).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/seed_types.py ℓ [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat
import forms6_flint as F

flint.ctx.threads = 2
L = int([a for a in sys.argv[1:] if a.isdigit()][0]); KEEP = '--keep' in sys.argv
QS = [q for q in (5, 7, 11, 13) if q != L]; LM = L - 1; NPW = 10
QX = [q for q in range(5, 98) if all(q % d for d in range(2, int(q ** 0.5) + 1)) and q != L]   # the Eisenstein test
OUT = os.path.join(HERE, 'seed_types_%d.txt' % L)
if KEEP and os.path.exists(OUT): sys.exit('stopped before computing anything: %s exists' % OUT)
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def sign(p, k): return leg(p, L) ** (1 + (k - (L - 3)) // (L - 1))


def split_pieces(Mq, d, Nn, qs):
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


def seed_pieces(w, s2, s3, N):
    """The class at weight w on Γ0(6) (new at 2 and 3 with these signs exactly), split by the systems of QS; for each
    piece, the characteristic polynomials mod ℓ of T_ℓ and of T_q for every q in QX."""
    n = w + 1; m = L ** N
    H = F.hecke_matrices(w, m, [('T', q) for q in QX] + [('T', L), ('U', 2), ('U', 3)])
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
    M = {q: on_new(F.restrict(H[('T', q)], Y, R, r)) for q in QX + [L]}
    out = []
    for lab, Yp, Rp, rp in split_pieces(M, d, N1, QS):
        cp = {q: [int(c) for c in F.reduce_mat(F.restrict(M[q], Yp, Rp, rp), L).charpoly().coeffs()] for q in QX + [L]}
        out.append((lab, rp, cp))
    return out


def twist(lab, m):
    out = []
    for q, v in zip(QS, lab):
        f = pow(q, m, L)
        out.append(v * f % L if isinstance(v, int) else tuple(v[i] * pow(f, len(v) - 1 - i, L) % L for i in range(len(v))))
    return tuple(out)


c2 = sign(2, L - 3) * pow(2, (L - 5) // 2, L) % L; c3 = sign(3, L - 3) * pow(3, (L - 5) // 2, L) % L
say('ℓ = %d: each family\'s seed, ordinary or not, Eisenstein or not (T_q for the %d primes q = 5 … 97 other than ℓ), '
    'and the theorems\' hypotheses' % (L, len(QX)))
rows = []; t0 = time.time()
for r in range(2, L, 2):
    for m in [m for m in range(LM) if (r + 2 * m - (L - 5)) % LM == 0]:
        u2 = pow(2, -m, L) * c2 % L; u3 = pow(3, -m, L) * c3 % L
        s2 = 1 if u2 * pow(2, -(r // 2), L) % L == 1 else -1; s3 = 1 if u3 * pow(3, -(r // 2), L) % L == 1 else -1
        for lab, rp, cp in seed_pieces(r + 2, s2, s3, NPW):
            fam = twist(lab, m)
            # ordinary: T_ℓ invertible mod ℓ on the piece (constant term of its characteristic polynomial ≠ 0)
            ordinary = cp[L][0] % L != 0
            lin = all(len(cp[q]) == 2 for q in QX)                     # all eigenvalues in F_ℓ (a 1-dimensional piece)
            eis = []
            if lin:
                ev = {q: (-cp[q][0]) % L for q in QX}
                eis = [(a_, b_) for a_ in range(LM) for b_ in range(a_, LM)
                       if all(ev[q] == (pow(q, a_, L) + pow(q, b_, L)) % L for q in QX)]
            rows.append((fam, r, m, rp, ordinary, eis, (-cp[L][0]) % L if len(cp[L]) == 2 else cp[L]))
say('   (%.0f s)' % (time.time() - t0))
say()
say('| family (T_q mod ℓ, q = %s) | seed | forms | a_ℓ of the seed mod ℓ | ordinary | Eisenstein (ρ̄ ≅ ω^a ⊕ ω^b, twisted) | '
    'Serre weight (a, b) | LTXZ Thm 1.3 | Wang (très ramifié) | BP Conj 3.1 (regular) |' % QS)
say('|---|---|---|---|---|---|---|---|---|---|')
for fam, r, m, rp, ordinary, eis, al in sorted(rows, key=lambda t: repr(t[0])):
    a, b = (0, m) if r == L - 1 else (r, m)
    tr = r == L - 1
    eis_t = [((a_ + m) % LM, (b_ + m) % LM) for a_, b_ in eis]           # ρ̄ = ρ̄_seed ⊗ ω^m
    irr = not eis
    ltxz = L >= 11 and irr and ordinary and not tr and 2 <= a <= L - 5
    why = [] if ltxz else [w for w, c in (('ℓ < 11', L < 11), ('globally reducible', not irr), ('not ordinary', not ordinary),
                                         ('a = 0 très ramifié', tr), ('a = %d not in 2 … ℓ − 5' % a, not tr and not 2 <= a <= L - 5)) if c]
    wang = L >= 11 and tr and ordinary
    say('| %r | weight %d·ω^%d | %d | %s | %s | %s | (%d, %d)%s | %s | %s | %s |' % (
        fam, r + 2, m, rp, al, 'yes' if ordinary else 'NO', ', '.join('ω^%d ⊕ ω^%d' % e for e in eis_t) or 'no', a, b,
        ' très ramifié' if tr else '', 'covered' if ltxz else 'no: ' + '; '.join(why), 'local theorem' if wang else '—',
        'regular' if ordinary else 'NOT regular'))
n = len(rows); cov = sum(1 for fam, r, m, rp, o, e, al in rows if L >= 11 and not e and o and r != L - 1 and 2 <= r <= L - 5)
say(); say('families: %d; seeds ordinary: %d; globally reducible: %d; covered by LTXZ Theorem 1.3: %d; Steinberg (Wang, '
           'local): %d' % (n, sum(1 for t in rows if t[4]), sum(1 for t in rows if t[5]), cov,
                           sum(1 for t in rows if t[1] == L - 1)))
