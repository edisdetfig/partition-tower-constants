# A dig, 12 Nov 2026: THE LUCKY PRIMES, COUNTED TO A MILLION. A law at a prime Q (p(N) ≡ 0 on lane ℓ for every N with
# 24N − 1 = Q²·M, ℓ | M, M in one class mod Q; or for every N with 24N − 1 = ℓ·Q³·m, Ono's kind) needs every form of ℓ's
# class to take one particular number on lane ℓ (README, the eleventh time): with λ_i = (3|Q)·a_Q(f_i) and
# δ = −(−3|Q)·Q⁻², the step c_i(Q²M) = (λ_i + εδ)·c_i(M) (ε = (M|Q)) and c_i(Q³m) = λ_i·c_i(Qm) give
#     a class law for ε = ±1  ⇔  λ_i ≡ −εδ for every i;      Ono's kind  ⇔  λ_i ≡ 0 for every i.
# PREDICTED before the run (derived from Chebotarev's density theorem, each piece of the class having the largest possible
# mod-ℓ Galois image and the pieces independent): a piece over the field with q = ℓ^f elements has the number c with chance
#     (q + χ_q(c² − 4D)) / (q² − 1),   D = Q^(ℓ−4) ≡ Q⁻³ its determinant, χ_q the quadratic character of that field
# (for f even every element of F_ℓ is a square there), and a law's chance is the product over the pieces: about 3/ℓ in all
# for 13 to 31 (one form), about 3/ℓ² for 37 to 59, about 3/ℓ³ for 61. The expected counts are summed prime by prime and
# printed before the counting; the counts should fall within a few √E of them. For one form the finer prediction: the
# share is 1/(ℓ − 1) where χ = +1 and 1/(ℓ + 1) where χ = −1. The laws known from partition numbers must be among those
# found: 37: 6599, 7541, 9547 (class) and 2999 (Ono's); 43: 4049, 8941 (class); none below 10,000 for 41, 47
# (partitions_split.txt, partitions_flint_arbiter.txt).
# Method (forms6_flint.py, FLINT only): the class at weight ℓ − 3 exactly (over the integers), its basis S_j with the
# identity at its pivot coefficients n_i, each S_j's q-expansion modulo ℓ to Q_max·max(n_i) terms; for each prime Q the
# matrix of T_Q on the class is M[i][j] = S_j(Q·n_i) + Q^(k−1)·S_j(n_i/Q), and a law needs (3|Q)·M ≡ c·I exactly.
# COUNTED (the laws), DERIVED (the expected numbers).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_lucky.py ℓ … [--max=1000000]
#          [--workers=3] [--keep]
import os, sys, time, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


KNOWN = {37: {'class': [6599, 7541, 9547], 'ono': [2999]}, 43: {'class': [4049, 8941], 'ono': []},
         41: {'class': [], 'ono': []}, 47: {'class': [], 'ono': []}}


def one(job):
    L, QMAX = job
    import flint
    from flint import nmod_poly, nmod_mat
    from sympy import primerange
    import forms6_flint as F
    flint.ctx.threads = 2
    try:
        import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
    except Exception: pass
    t0 = time.time(); k = L - 3; n0 = k + 1; out = []
    say = out.append
    # the class, exactly, and its pieces on lane ℓ
    Kz, d = F.class_basis_exact(k, leg(2, L), leg(3, L))
    Xc, Rc = F.normalized_basis(nmod_mat(n0, d, [int(Kz[i, j]) % L for i in range(n0) for j in range(d)], L), d, L)
    H = F.hecke_matrices(k, L, [('T', q) for q in (5, 7, 11, 13)])
    A = None
    for c, q in zip((1, 3, 7, 12), (5, 7, 11, 13)):
        Tq = F.restrict(H[('T', q)], Xc, Rc, d) * c
        A = Tq if A is None else A + Tq
    fac = A.charpoly().factor()[1]; pieces = sorted(f.degree() for f, e in fac for _ in range(e))
    rep = any(e > 1 for f, e in fac)
    say('ℓ = %d (weight %d): the class has %d form%s; on lane ℓ its pieces have degrees %s%s; pivot coefficients %s' % (
        L, k, d, '' if d == 1 else 's', pieces, ' (a REPEATED piece: the model below does not apply as stated)' if rep else '', Rc))
    # q-expansions of the class basis, modulo ℓ, to Lq terms (H_m = q^m·φ·u^(k−m); Horner in u)
    Lq = QMAX * max(Rc) + 1
    E1, E2, E3, E6 = (F.eta_series(dd, Lq, L) for dd in (1, 2, 3, 6))
    mul = lambda a, b: a.mul_low(b, Lq)
    u = mul(mul(E1.pow_trunc(5, Lq), E3), mul(E2, E6.pow_trunc(5, Lq)).inverse_series_trunc(Lq))
    phi = mul(mul(E1, E6.pow_trunc(6, Lq)), mul(E2.pow_trunc(2, Lq), E3.pow_trunc(3, Lq)).inverse_series_trunc(Lq)).pow_trunc(k, Lq)
    short = lambda p: nmod_poly([int(p[i]) for i in range(min(n0, p.length()))], L)
    Bs = [[0] * n0 for _ in range(n0)]; psi = short(phi); us = short(u)
    for i in range(n0):                                       # the basis's own first k + 1 coefficients (B lower unitriangular)
        mm = k - i
        for nn in range(mm, n0): Bs[nn][mm] = int(psi[nn - mm]) if nn - mm < psi.length() else 0
        if i < n0 - 1: psi = psi.mul_low(us, n0)
    W = nmod_mat(n0, n0, [x for row in Bs for x in row], L).inv() * Xc     # coordinates of S_j in the basis H_m
    S = []
    for j in range(d):
        R = nmod_poly([int(W[0, j])], L)
        for jj in range(1, n0):
            R = mul(R, u); c = int(W[jj, j])
            if c: R[jj] = (int(R[jj]) + c) % L if jj < R.length() else c
        S.append(mul(R, phi))
    del u, phi, E1, E2, E3, E6
    for j in range(d):                                        # check: S_j has the class's first coefficients
        if any(int(S[j][nn]) != int(Xc[nn, j]) for nn in range(n0)): say('   the long series do not match the class; stopped'); return out
    co = [S[j].coeffs() for j in range(d)]; del S
    get = lambda j, nn: int(co[j][nn]) if nn < len(co[j]) else 0
    say('   q-expansions to %d terms (%.0f s)' % (Lq, time.time() - t0))
    # every prime Q
    t = time.time(); pk = lambda Q: pow(Q, k - 1, L)
    TYPES = ('class ε=+1', 'class ε=−1', 'Ono')
    E = {ty: 0.0 for ty in TYPES}; C = {ty: [] for ty in TYPES}; Echi = {}; Cchi = {}; hit1 = 0; Ehit1 = 0.0; NQ = 0
    for Q in primerange(5, QMAX + 1):
        if Q == L: continue
        NQ += 1; s3 = leg(3, Q); D = pow(Q, -3, L); iQ2 = pow(Q * Q, -1, L)
        M = [[(s3 * (get(j, Q * Rc[i]) + (pk(Q) * get(j, Rc[i] // Q) if Rc[i] % Q == 0 else 0))) % L for j in range(d)] for i in range(d)]
        for ty, cc in zip(TYPES, (leg(-3, Q) * iQ2 % L, (-leg(-3, Q) * iQ2) % L, 0)):
            x = (cc * cc - 4 * D) % L; p = 1.0; chis = []
            for f in pieces:
                qf = L ** f; ch = 0 if x == 0 else (1 if f % 2 == 0 else leg(x, L)); chis.append(ch)
                p *= (qf + ch) / (qf * qf - 1)
            E[ty] += p
            key = (ty, tuple(chis)); Echi[key] = Echi.get(key, 0.0) + p; Cchi.setdefault(key, [0, 0]); Cchi[key][1] += 1
            if all(M[i][j] == (cc if i == j else 0) for i in range(d) for j in range(d)):
                C[ty].append(Q); Cchi[key][0] += 1
        A0 = nmod_mat(d, d, [M[i][j] for i in range(d) for j in range(d)], L)   # the marginal: SOME piece has the number 0
        if A0.det() == 0: hit1 += 1
        x0 = (-4 * D) % L; pr = 1.0
        for f in pieces: qf = L ** f; pr *= 1 - (qf + (1 if f % 2 == 0 else leg(x0, L))) / (qf * qf - 1)
        Ehit1 += 1 - pr
    say('   %d primes Q from 5 to %d (ℓ left out)  (%.0f s)' % (NQ, QMAX, time.time() - t))
    for ty in TYPES:
        z = (len(C[ty]) - E[ty]) / math.sqrt(E[ty]) if E[ty] > 0 else float('nan')
        say('   %-11s expected %9.1f (derived)   counted %6d   (counted − expected)/√expected = %+.2f' % (ty, E[ty], len(C[ty]), z))
        small = [Q for Q in C[ty] if Q < 10000]
        say('      below 10,000: %s%s' % (small[:60], ' …' if len(small) > 60 else ''))
        if d >= 2 and len(C[ty]) <= 400: say('      all: %s' % C[ty])
    if d == 1:
        for (ty, chis), (cnt, nq) in sorted(Cchi.items()):
            say('   %-11s where χ = %+d: %6d of %6d primes = 1/%.2f  (predicted 1/%d)' % (ty, chis[0], cnt, nq, nq / cnt if cnt else float('inf'), L - chis[0]))
    say('   some piece with the number 0 (one form is enough): counted %d, expected %.1f (derived)' % (hit1, Ehit1))
    if L in KNOWN:
        for kind, lst in KNOWN[L].items():
            have = C['Ono'] if kind == 'ono' else C['class ε=+1'] + C['class ε=−1']
            for Q in lst: say('   known from partition numbers, %s law at %d: %s' % (kind, Q, 'found' if Q in have else 'NOT FOUND'))
        lo = sorted(Q for ty in TYPES for Q in C[ty] if Q < 10000)
        say('   all laws below 10,000 here: %s  (partitions_split.txt had: %s)' % (lo, sorted(KNOWN[L]['class'] + KNOWN[L]['ono'])))
    say('   seconds: %.0f' % (time.time() - t0))
    return out


if __name__ == '__main__':
    import multiprocessing as mp
    PS = [int(a) for a in sys.argv[1:] if a.isdigit()]
    QMAX = next((int(a[6:]) for a in sys.argv[1:] if a.startswith('--max=')), 10 ** 6)
    WK = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--workers=')), 3); KEEP = '--keep' in sys.argv
    fh = open(os.path.join(HERE, 'partitions_flint_lucky_%s.txt' % '_'.join(map(str, PS))), 'x', encoding='utf-8') if KEEP else None

    def say(*a):
        line = ' '.join(str(x) for x in a); print(line, flush=True)
        if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())

    t0 = time.time(); say('primes ℓ: %s; laws counted at every prime Q up to %d' % (PS, QMAX))
    with mp.Pool(WK) as pool:
        for lines in pool.imap(one, [(L, QMAX) for L in PS]):
            say()
            for line in lines: say(line)
    say(); say('seconds: %.0f' % (time.time() - t0))
