# A dig, 13 Nov 2026: THE LUCKY PRIMES BY A LIGHTER TEST, and the scan for 73 (README "Lucky primes for 73"; predicted
# there first). A law at a prime Q needs (3|Q)·T_Q ≡ c on ℓ's class modulo ℓ, with c = ε·(−3|Q)·Q⁻² for a class law
# (ε = ±1) or c = 0 for Ono's kind (partitions_flint_lucky.py, which tested the whole matrix of T_Q).
# THE LIGHTER TEST. a_n(f) = a_1(T_n·f) for every form f of the class and every n, so an operator t of the Hecke algebra
# with a_1(t·f) ≡ 0 for every f of the class has every coefficient of t·f ≡ 0: t is 0 on the class modulo ℓ. With
# t = (3|Q)·T_Q − c the law needs only (3|Q)·a_Q(S_j) ≡ c·a_1(S_j) for the class's basis S_j: q-expansions to Q_max terms.
# THE SERIES, modulo ℓ. forms6_flint.py's basis is H_m = q^m·φ·u^(k−m) = F^k·w^m, with F = η(z)⁶η(6z)/(η(2z)³η(3z)²)
# (weight 1, constant term 1) and w = q·E2·E6⁵/(E1⁵·E3) (= 1/t, starts with q). So S_j = F^k·P_j(w), with P_j of degree
# ≤ k whose coefficients are S_j's coordinates in that basis; P_j(w) is evaluated by baby steps and giant steps
# (Paterson–Stockmeyer), about 4 times fewer long products than Horner's rule. Checked: S_j's first k + 1 coefficients
# are the class's.
# THE EXPECTED NUMBERS (derived, Chebotarev, as in partitions_flint_lucky.py): a piece over the field of q = ℓ^f elements
# has the number c with chance (q + χ_q(c² − 4D))/(q² − 1), D = Q⁻³ its determinant; the pieces independent.
# --compare: the lists and counts must equal those of partitions_flint_lucky_37_43_41_47_53_59_61_13_17_19_23_29_31.txt.
# COUNTED (the laws), DERIVED (the expected numbers).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_lucky_dual.py ℓ …
#          [--max=1000000] [--workers=2] [--threads=3] [--compare] [--keep]
import os, sys, time, math, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OLD = os.path.join(HERE, 'partitions_flint_lucky_37_43_41_47_53_59_61_13_17_19_23_29_31.txt')
TYPES = ('class ε=+1', 'class ε=−1', 'Ono')


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def primes_upto(n):
    import numpy as np
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]


def one(job):
    L, QMAX, TH = job
    import numpy as np
    import flint
    from flint import nmod_poly, nmod_mat
    import forms6_flint as F
    flint.ctx.threads = TH
    try:
        import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
    except Exception: pass
    t0 = time.time(); k = L - 3; n0 = k + 1; out = []; say = out.append
    # ---- the class, exactly, and its pieces on lane ℓ ----
    Kz, d = F.class_basis_exact(k, leg(2, L), leg(3, L))
    X0 = nmod_mat(n0, d, [int(Kz[i, j]) % L for i in range(n0) for j in range(d)], L)
    if X0.rank() != d: say('ℓ = %d: the exact class does not stay of full rank on lane ℓ; stopped' % L); return L, out, None
    Xc, Rc = F.normalized_basis(X0, d, L)
    H = F.hecke_matrices(k, L, [('T', q) for q in (5, 7, 11, 13)])
    A = None
    for c, q in zip((1, 3, 7, 12), (5, 7, 11, 13)):
        Tq = F.restrict(H[('T', q)], Xc, Rc, d) * c
        A = Tq if A is None else A + Tq
    fac = A.charpoly().factor()[1]; pieces = sorted(f.degree() for f, e in fac for _ in range(e))
    rep = any(e > 1 for f, e in fac)
    say('ℓ = %d (weight %d): the class has %d form%s; on lane ℓ its pieces have degrees %s%s; pivot coefficients %s' % (
        L, k, d, '' if d == 1 else 's', pieces, ' (a REPEATED piece: the model below does not apply as stated)' if rep else '', Rc))
    if 1 not in Rc: say('   the coefficient q¹ is not a pivot; stopped'); return L, out, None
    # ---- the series, modulo ℓ, to Q_max + 1 terms ----
    Nl = QMAX + 1; mul = lambda a, b: a.mul_low(b, Nl)
    E1, E2, E3, E6 = (F.eta_series(dd, Nl, L) for dd in (1, 2, 3, 6))
    E1_5 = E1.pow_trunc(5, Nl)
    w = mul(mul(E2, E6.pow_trunc(5, Nl)), mul(E1_5, E3).inverse_series_trunc(Nl)).left_shift(1).truncate(Nl)
    Fz = mul(mul(mul(E1_5, E1), E6), mul(E2.pow_trunc(3, Nl), E3.pow_trunc(2, Nl)).inverse_series_trunc(Nl))
    del E1, E2, E3, E6, E1_5
    Fk = Fz.pow_trunc(k, Nl); del Fz
    at = lambda p, i: int(p[i]) if i < p.length() else 0
    Bs = [[0] * n0 for _ in range(n0)]; h = Fk.truncate(n0); ws = w.truncate(n0)
    for mm in range(n0):                                       # B[n][m] = coefficient of q^n in H_m = F^k·w^m
        for nn in range(n0): Bs[nn][mm] = at(h, nn)
        h = h.mul_low(ws, n0)
    W = nmod_mat(n0, n0, [x for row in Bs for x in row], L).inv() * Xc     # S_j's coordinates in the basis H_m
    s = 8; baby = [nmod_poly([1], L), w]
    for b in range(2, s): baby.append(mul(baby[-1], w))
    g = mul(baby[-1], w); nb = (n0 + s - 1) // s
    P = [int(q) for q in primes_upto(QMAX) if q >= 5 and q != L]
    vals = []
    for j in range(d):
        acc = None
        for a in range(nb - 1, -1, -1):
            blk = nmod_poly([], L)
            for b in range(s):
                m = a * s + b
                if m < n0 and int(W[m, j]): blk = blk + baby[b] * int(W[m, j])
            acc = blk if acc is None else mul(acc, g) + blk
        Sj = mul(Fk, acc); del acc
        if any(at(Sj, nn) != int(Xc[nn, j]) for nn in range(n0)):
            say('   the long series do not match the class; stopped'); return L, out, None
        vals.append(np.array([at(Sj, q) for q in P], dtype=np.int64)); del Sj
    del baby, g, w, Fk
    say('   q-expansions of the %d basis forms to %d terms (%.0f s)' % (d, Nl, time.time() - t0))
    # ---- every prime Q: the lighter test, and the expected numbers ----
    t = time.time(); Pa = np.array(P, dtype=np.int64); V = np.stack(vals); r = Pa % L
    s3 = np.where((Pa % 12 == 1) | (Pa % 12 == 11), 1, -1); sm3 = np.where(Pa % 3 == 1, 1, -1)
    inv2 = np.array([0] + [pow(x, -2, L) for x in range(1, L)], dtype=np.int64)[r]
    D = np.array([0] + [pow(x, -3, L) for x in range(1, L)], dtype=np.int64)[r]
    legt = np.array([leg(x, L) for x in range(L)], dtype=np.int64)
    a1 = np.array([int(Xc[1, j]) for j in range(d)], dtype=np.int64)
    lhs = (s3[None, :] * V) % L
    res = {}
    for ty, cc in zip(TYPES, ((sm3 * inv2) % L, (-sm3 * inv2) % L, np.zeros_like(Pa))):
        ok = np.all((lhs - cc[None, :] * a1[:, None]) % L == 0, axis=0)
        x = (cc * cc - 4 * D) % L; pr = np.ones(len(Pa)); chis = []
        for f in pieces:
            qf = float(L) ** f; ch = np.where(x == 0, 0, 1 if f % 2 == 0 else legt[x]); chis.append(ch)
            pr *= (qf + ch) / (qf * qf - 1)
        Qs = [int(q) for q in Pa[ok]]; E = float(pr.sum())
        res[ty] = {'all': Qs, 'count': len(Qs), 'expected': E}
        z = (len(Qs) - E) / math.sqrt(E) if E > 0 else float('nan')
        say('   %-11s expected %9.1f (derived)   counted %6d   (counted − expected)/√expected = %+.2f' % (ty, E, len(Qs), z))
        small = [q for q in Qs if q < 10000]
        say('      below 10,000: %s%s' % (small[:60], ' …' if len(small) > 60 else ''))
        if d >= 2 and len(Qs) <= 400: say('      all: %s' % Qs)
        if d == 1:
            for sg in (1, -1):
                sel = chis[0] == sg; cnt = int(ok[sel].sum()); nq = int(sel.sum())
                say('      where χ = %+d: %6d of %6d primes = 1/%.2f  (predicted 1/%d)' % (sg, cnt, nq, nq / cnt if cnt else float('inf'), L - sg))
    say('   %d primes Q from 5 to %d (ℓ left out)  (%.0f s);  seconds in all: %.0f' % (len(P), QMAX, time.time() - t, time.time() - t0))
    return L, out, res


def old_lists():
    res = {}; L = ty = None
    for line in open(OLD, encoding='utf-8'):
        m = re.match(r'ℓ = (\d+) \(weight', line)
        if m: L = int(m.group(1)); res[L] = {}; ty = None; continue
        m = re.match(r'\s+(class ε=\+1|class ε=−1|Ono)\s+expected\s+([\d.]+).*counted\s+(\d+)', line)
        if m and L: ty = m.group(1); res[L][ty] = {'count': int(m.group(3)), 'expected': float(m.group(2))}; continue
        m = re.match(r'\s+all: \[(.*)\]', line)
        if m and L and ty: res[L][ty]['all'] = [int(x) for x in m.group(1).split(',') if x.strip()]
    return res


if __name__ == '__main__':
    import multiprocessing as mp
    PS = [int(a) for a in sys.argv[1:] if a.isdigit()]
    QMAX = next((int(a[6:]) for a in sys.argv[1:] if a.startswith('--max=')), 10 ** 6)
    WK = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--workers=')), 2)
    TH = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--threads=')), 3)
    KEEP = '--keep' in sys.argv; COMPARE = '--compare' in sys.argv
    name = 'partitions_flint_lucky_dual_%s_max%d.txt' % ('_'.join(map(str, PS)), QMAX)
    fh = open(os.path.join(HERE, name), 'x', encoding='utf-8') if KEEP else None

    def say(*a):
        line = ' '.join(str(x) for x in a); print(line, flush=True)
        if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())

    t0 = time.time(); say('primes ℓ: %s; the lighter test at every prime Q up to %d' % (PS, QMAX))
    old = old_lists() if COMPARE else {}
    same = []
    with mp.Pool(WK) as pool:
        for L, lines, res in pool.imap(one, [(L, QMAX, TH) for L in PS]):
            say()
            for line in lines: say(line)
            if COMPARE and res is not None and L in old:
                for ty in TYPES:
                    o = old[L].get(ty, {}); ok_c = o.get('count') == res[ty]['count'] if QMAX == 10 ** 6 else None
                    ok_l = None if 'all' not in o else (o['all'] == [q for q in res[ty]['all'] if q <= 10 ** 6])
                    same.append(ok_c is not False and ok_l is not False)
                    say('   against the whole-matrix test of 12 Nov, %-11s: count %s, list %s' % (
                        ty, {True: 'the same', False: 'DIFFERENT', None: '—'}[ok_c], {True: 'the same', False: 'DIFFERENT', None: '(not printed then)'}[ok_l]))
    if COMPARE: say(); say('against the whole-matrix test: %d of %d comparisons the same' % (sum(same), len(same)))
    say(); say('seconds: %.0f' % (time.time() - t0))
