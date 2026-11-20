# A dig, 13 Nov 2026: IS 97'S SINGLE NOTE TOO QUIET? (README "Pictures of the chords", the test of 97's note; predicted
# there first.) Up to 10⁶ the single form g of 97's class was silent at 2,326 primes against 2,427.8 expected (z = −2.07).
# Here, for every prime Q < 10⁷: λ(Q) = (3|Q)·a_Q(f) on lane ℓ for 97's g and, as a control, for 37's twins A and B;
#   - its spread over the ℓ residues against Chebotarev's weights: value c with chance (ℓ + χ(c² − 4Q⁻³))/(ℓ² − 1), the
#     forms' mod-ℓ Galois images taken as the full ones; χ² on ℓ − 1 degrees of freedom;
#   - the silences (the three law kinds), counted against expected, up to 10⁶ and up to 10⁷.
# a_Q(f) = (Σ v_j·S_j(Q)) / (Σ v_j·S_j(1)), v the T_5-eigenvector of the class on lane ℓ (chords_page_data.py).
# COUNTED (the spreads and silences), DERIVED (the expected numbers).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/note97_check.py [--keep]
import os, sys, time, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import numpy as np
import flint
import forms6_flint as F
from flint import nmod_poly, nmod_mat
flint.ctx.threads = 4
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'note97_check.txt'), 'x', encoding='utf-8') if KEEP else None
QMAX = 10 ** 7


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def class_series(L, Nl):
    """The class basis S_j modulo ℓ to Nl terms (as in partitions_flint_lucky_dual.py), with Xc, Rc."""
    k = L - 3; n0 = k + 1
    Kz, d = F.class_basis_exact(k, leg(2, L), leg(3, L))
    Xc, Rc = F.normalized_basis(nmod_mat(n0, d, [int(Kz[i, j]) % L for i in range(n0) for j in range(d)], L), d, L)
    mul = lambda a, b: a.mul_low(b, Nl)
    E1, E2, E3, E6 = (F.eta_series(dd, Nl, L) for dd in (1, 2, 3, 6))
    E1_5 = E1.pow_trunc(5, Nl)
    w = mul(mul(E2, E6.pow_trunc(5, Nl)), mul(E1_5, E3).inverse_series_trunc(Nl)).left_shift(1).truncate(Nl)
    Fk = mul(mul(mul(E1_5, E1), E6), mul(E2.pow_trunc(3, Nl), E3.pow_trunc(2, Nl)).inverse_series_trunc(Nl)).pow_trunc(k, Nl)
    at = lambda p, i: int(p[i]) if i < p.length() else 0
    Bs = [[0] * n0 for _ in range(n0)]; h = Fk.truncate(n0); ws = w.truncate(n0)
    for mm in range(n0):
        for nn in range(n0): Bs[nn][mm] = at(h, nn)
        h = h.mul_low(ws, n0)
    W = nmod_mat(n0, n0, [x for row in Bs for x in row], L).inv() * Xc
    S = []
    for j in range(d):
        acc = nmod_poly([int(W[k, j])], L)
        for m in range(k - 1, -1, -1): acc = mul(acc, w) + nmod_poly([int(W[m, j])], L)
        Sj = mul(Fk, acc); assert all(at(Sj, nn) == int(Xc[nn, j]) for nn in range(n0))
        S.append(Sj)
    return S, Xc, Rc, d


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]


t0 = time.time()
for L, notes in ((97, (('g', 58),)), (37, (('A', 36), ('B', 5)))):
    t = time.time(); S, Xc, Rc, d = class_series(L, QMAX + 1)
    P = [int(q) for q in primes_upto(QMAX) if q >= 5 and q != L]; Pa = np.array(P, dtype=np.int64)
    k = L - 3; H = F.hecke_matrices(k, L, [('T', 5)]); T5 = F.restrict(H[('T', 5)], Xc, Rc, d)
    vals = [np.array([int(S[j][q]) for q in P], dtype=np.int64) for j in range(d)]
    s3 = np.where((Pa % 12 == 1) | (Pa % 12 == 11), 1, -1); sm3 = np.where(Pa % 3 == 1, 1, -1)
    r = Pa % L
    inv2 = np.array([0] + [pow(x, -2, L) for x in range(1, L)], dtype=np.int64)[r]
    D3 = np.array([0] + [pow(x, -3, L) for x in range(1, L)], dtype=np.int64)[r]
    legt = np.array([leg(x, L) for x in range(L)], dtype=np.int64)
    # Chebotarev weights for every value c, summed over the primes (by residue class of Q)
    cls = np.bincount(r, minlength=L)
    for name, a5 in notes:
        Kv, nul = (T5 - F.identity(d, L) * a5).nullspace(); assert nul == 1
        v = [int(Kv[i, 0]) for i in range(d)]
        a1 = sum(v[j] * int(S[j][1]) for j in range(d)) % L; inv = pow(a1, -1, L)
        aQ = (sum(v[j] * vals[j] for j in range(d)) % L) * inv % L
        lam = (s3 * aQ) % L
        for bound in (10 ** 6, QMAX):
            sel = Pa <= bound; lam_b = lam[sel]; r_b = r[sel]; cls_b = np.bincount(r_b, minlength=L)
            obs = np.bincount(lam_b, minlength=L).astype(float)
            exp = np.zeros(L)
            for rr in range(1, L):
                if cls_b[rr] == 0: continue
                Dr = pow(rr, -3, L)
                x = (np.arange(L) ** 2 - 4 * Dr) % L
                ch = np.where(x == 0, 0, legt[x])
                exp += cls_b[rr] * (L + ch) / (L * L - 1)
            chi2 = float(((obs - exp) ** 2 / exp).sum())
            tgt = (sm3[sel] * inv2[sel]) % L
            sil = {'class ε=+1': int((lam_b == tgt).sum()), 'class ε=−1': int((lam_b == (-tgt) % L).sum()), 'Ono': int((lam_b == 0).sum())}
            esil = {}
            for kind, cc in (('class ε=+1', tgt), ('class ε=−1', (-tgt) % L), ('Ono', np.zeros_like(tgt))):
                x = (cc * cc - 4 * D3[sel]) % L; ch = np.where(x == 0, 0, legt[x])
                esil[kind] = float(((L + ch) / (L * L - 1)).sum())
            ts, te = sum(sil.values()), sum(esil.values())
            say('ℓ = %d, %s, primes Q < %.0e (%d): spread over the %d residues χ² = %.1f on %d degrees of freedom (z = %+.2f);  '
                'silences %d against %.1f expected (z = %+.2f; by kind %s)' % (
                    L, name, bound, int(sel.sum()), L, chi2, L - 1, (chi2 - (L - 1)) / math.sqrt(2 * (L - 1)), ts, te,
                    (ts - te) / math.sqrt(te), ', '.join('%+.2f' % ((sil[kk] - esil[kk]) / math.sqrt(esil[kk])) for kk in sil)))
    say('   (%.0f s)' % (time.time() - t))
say('seconds: %.0f' % (time.time() - t0))
