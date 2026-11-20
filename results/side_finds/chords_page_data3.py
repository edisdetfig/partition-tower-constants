# A dig, 13 Nov 2026: 97'S BRAIDED TRIPLE, ITS SILENCES (for the page "Towers as Chords", the silences split apart; README
# "Pictures of the chords"; predicted there first). The triple's eigenform f₁ has its numbers in the field with 97³
# elements: a_Q(f₁) = (Σ w_j·S_j(Q)) / (Σ w_j·S_j(1)), w the T₅-eigenvector of the class on lane 97 for a root μ of the
# cubic factor of T₅'s characteristic polynomial (the class basis S_j as in chords_page_data.py). The triple is silent at Q
# when (3|Q)·a_Q(f₁) is the law's number c in F_97 (then its conjugates are too). PREDICTED: about 0.26 such primes below
# 10⁶ (Chebotarev in the field with q = 97³ elements: (q + χ(c² − 4Q⁻³))/(q² − 1) per kind), most likely none; so the
# whole tower of 97 (all four notes silent together) is silent at none. Writes chords_page_data3.json. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/chords_page_data3.py [--keep]
import os, sys, json, time, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly, nmod_mat
import forms6_flint as F
from sympy import primerange
flint.ctx.threads = 4
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'chords_page_data3.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def class_series(L, Nl):
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


t0 = time.time(); L = 97; QMAX = 10 ** 6
S, Xc, Rc, d = class_series(L, QMAX + 1)
k = L - 3; H = F.hecke_matrices(k, L, [('T', 5)]); T5 = F.restrict(H[('T', 5)], Xc, Rc, d)
cp = T5.charpoly(); cub = next(f for f, e in cp.factor()[1] if f.degree() == 3)
A = [int(u) for u in cub.coeffs()]                            # F_97[y]/(A), μ = y a root


def fmul(a, b):
    pr = [0] * 5
    for i in range(3):
        if a[i]:
            for j in range(3): pr[i + j] += a[i] * b[j]
    for t in (4, 3):
        cc = pr[t] % L
        if cc:
            for i in range(3): pr[t - 3 + i] -= cc * A[i]
        pr[t] = 0
    return tuple(pr[i] % L for i in range(3))


def fpow(a, e):
    r = (1, 0, 0)
    while e:
        if e & 1: r = fmul(r, a)
        a = fmul(a, a); e >>= 1
    return r


fadd = lambda a, b: tuple((u + v) % L for u, v in zip(a, b)); fsub = lambda a, b: tuple((u - v) % L for u, v in zip(a, b))
fsc = lambda c, a: tuple((c * u) % L for u in a); finv = lambda a: fpow(a, L ** 3 - 2); Z = (0, 0, 0)
mu = (0, 1, 0)
Mx = [[fsub((int(T5[i, j]) % L, 0, 0), mu if i == j else Z) for j in range(d)] for i in range(d)]   # T5 − μ
# kernel of Mx over the field (Gaussian elimination), one-dimensional
rows = [r[:] for r in Mx]; piv = []; rk = 0
for cc in range(d):
    pr = next((i for i in range(rk, d) if rows[i][cc] != Z), None)
    if pr is None: continue
    rows[rk], rows[pr] = rows[pr], rows[rk]; iv = finv(rows[rk][cc]); rows[rk] = [fmul(iv, x) for x in rows[rk]]
    for i in range(d):
        if i != rk and rows[i][cc] != Z:
            f = rows[i][cc]; rows[i] = [fsub(x, fmul(f, y)) for x, y in zip(rows[i], rows[rk])]
    piv.append(cc); rk += 1
free = [cc for cc in range(d) if cc not in piv]; assert len(free) == 1
wv = [Z] * d; wv[free[0]] = (1, 0, 0)
for i, cc in enumerate(piv): wv[cc] = fsub(Z, rows[i][free[0]])
a1 = Z
for j in range(d): a1 = fadd(a1, fsc(int(S[j][1]), wv[j]))
ia1 = finv(a1)
primes = [int(Q) for Q in primerange(5, QMAX + 1) if Q != L]
sil = {'p': [], 'm': [], 'o': []}; E = 0.0; q = L ** 3
for Q in primes:
    aQ = Z
    for j in range(d):
        c = int(S[j][Q]) if Q < S[j].length() else 0
        if c: aQ = fadd(aQ, fsc(c, wv[j]))
    lam = fsc(leg(3, Q), fmul(aQ, ia1)); t = leg(-3, Q) * pow(Q, -2, L) % L; D3 = pow(Q, -3, L)
    for kind, cc in (('p', t), ('m', (-t) % L), ('o', 0)):
        x = (cc * cc - 4 * D3) % L; ch = 0 if x == 0 else leg(x, L)
        E += (q + ch) / (q * q - 1)
        if lam == (cc, 0, 0): sil[kind].append(Q)
n = sum(len(v) for v in sil.values())
say('97\'s braided triple, primes Q from 5 to %d (%d): silent at %d (class +1 %s, class −1 %s, Ono %s); expected %.2f (derived)  (%.0f s)' % (
    QMAX, len(primes), n, sil['p'], sil['m'], sil['o'], E, time.time() - t0))
json.dump({'97': {'triple': sil, 'expected': E, 'primes': len(primes)}}, open(os.path.join(HERE, 'chords_page_data3.json'), 'w'))
say('wrote chords_page_data3.json')
