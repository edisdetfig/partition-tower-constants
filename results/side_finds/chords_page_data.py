# A dig, 13 Nov 2026: THE DATA BEHIND THE PAGE "Towers as Chords" (README "Pictures of the chords"; predicted there first).
#  A  37's tower on 37's cycle: the whole tower H_b(n) = c(M), M = 37^b·(24n + 11), and each twin's tower alone
#     G_i(M) = c(25M) + (κ − λ_j(5))·c(M) (partitions_twins_split.py), levels b = 1, 3, 5, 7, 9, n = 0 … 23 (5 ∤ M), as far
#     as p(N) with N ≤ 1.2·10¹⁶. The angle of a value is its power of the primitive root 2 (the page "Every prime at zero").
#  B  the multipliers' digits (modulo 37⁵) between levels b and b + 2: each twin's, and the whole tower's.
#  C  73's braided strands: G₁(M) in the field with 73³ elements (partitions_braid_split.py, on lane 73), levels 1 and 3;
#     strands 2 and 3 are its Frobenius images; angles are discrete logarithms in that field's cycle of 389,016.
#  D  silences: every prime Q < 10⁶ where one twin of 37 alone, or 97's single form alone, takes a law's number
#     ((3|Q)·a_Q(f) ≡ ε·(−3|Q)·Q⁻², or 0); the whole tower's laws are the Q where all forms agree on one kind.
# Writes chords_page_data.json (for the page) and prints the checks. COUNTED (the checks), DERIVED (κ, λ, logarithms).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/chords_page_data.py [--keep]
import os, sys, re, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly, nmod_mat
import forms6_flint as F
flint.ctx.threads = 4
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'chords_page_data.txt'), 'x', encoding='utf-8') if KEEP else None
NMAX = 12 * 10 ** 15
OUT = {}


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


PC = {}


def c(M):
    N = (M + 1) // 24
    if N not in PC: PC[N] = int(flint.fmpz(N).partitions_p())
    return PC[N]


def digits(x, L, e): return [(x // L ** i) % L for i in range(e)]


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


def eigenform_coeffs(L, S, d, Xc, Rc, a5, primes):
    """a_Q(f) for the eigenform f of the class with a_5(f) ≡ a5 (a one-dimensional eigenspace of T_5 on lane ℓ)."""
    k = L - 3; H = F.hecke_matrices(k, L, [('T', 5)]); T5 = F.restrict(H[('T', 5)], Xc, Rc, d)
    Kv, nul = (T5 - F.identity(d, L) * a5).nullspace(); assert nul == 1
    v = [int(Kv[i, 0]) for i in range(d)]
    a1 = sum(v[j] * int(S[j][1]) for j in range(d)) % L; inv = pow(a1, -1, L)
    return {Q: sum(v[j] * (int(S[j][Q]) if Q < S[j].length() else 0) for j in range(d)) * inv % L for Q in primes}


def silences(L, aQ, Q):
    """Which law's number (3|Q)·a_Q(f) takes at Q: 'p' (class ε = +1), 'm' (class ε = −1), 'o' (Ono's kind), or None."""
    lam = leg(3, Q) * aQ % L; t = leg(-3, Q) * pow(Q, -2, L) % L
    return 'p' if lam == t else 'm' if lam == (-t) % L else 'o' if lam == 0 else None


t0 = time.time()
from sympy import primerange
# ---------------- A and B: 37 ----------------
L = 37; g37 = 2; m5 = L ** 5
dlog = {pow(g37, i, L): i for i in range(L - 1)}
steps = open(os.path.join(HERE, 'partitions_flint_family_steps_37.txt'), encoding='utf-8').read()
alphas = [[int(x) for x in g.split(',')] for g in re.findall(r'a constant on lane ℓ, digits lowest first: \[([^\]]*)\]', steps)]
lams = [[int(x) for x in g.split(',')] for g in re.findall(r'a root on lane ℓ, digits lowest first: \[([^\]]*)\]', steps)]
num = lambda dg, mod: sum(x * L ** i for i, x in enumerate(dg)) % mod
# pairing from partitions_twins_split.txt: λ(5) 32… with the constant 34…, λ(5) 1… with 15…
twins = [{'name': 'A', 'alpha': next(a for a in alphas if a[0] == 15), 'lam': next(l for l in lams if l[0] == 1)},
         {'name': 'B', 'alpha': next(a for a in alphas if a[0] == 34), 'lam': next(l for l in lams if l[0] == 32)}]
ns = [n for n in range(24) if (24 * n + 11) % 5]
levels = {}; mults = {}
for b in (1, 3, 5, 7, 9):
    t = time.time(); row = {'n': [], 'chord': [], 'A': [], 'B': []}
    for n in ns:
        M = L ** b * (24 * n + 11)
        if (25 * M + 1) // 24 > NMAX: break
        cM, c25 = c(M), c(25 * M)
        kap = leg(M, 5) * leg(-3, 5) * pow(5, -2, m5) % m5
        row['n'].append(n); row['chord'].append(cM % m5)
        for i, tw in enumerate(twins):
            lj = num(twins[1 - i]['lam'], m5)
            row[tw['name']].append((c25 + (kap - lj) * cM) % m5)
    levels[b] = row
    say('37, level %d: %d values of n (%.0f s)' % (b, len(row['n']), time.time() - t))
ang = lambda v: None if v % L == 0 else dlog[v % L]
A_out = {'ell': 37, 'root': g37, 'levels': {}, 'alpha': {tw['name']: tw['alpha'][:5] for tw in twins}}
for b, row in levels.items():
    A_out['levels'][b] = {'n': row['n'], 'chord': [ang(v) for v in row['chord']], 'A': [ang(v) for v in row['A']], 'B': [ang(v) for v in row['B']]}
# the turns and the multipliers' digits
B_out = {}
for b in (1, 3, 5, 7):
    r0, r2 = levels[b], levels.get(b + 2)
    if not r2 or not r2['n']: continue
    common = [n for n in r2['n'] if n in r0['n']]
    pair = {'n': common, 'chord': [], 'A': [], 'B': []}; turns = {'chord': set(), 'A': set(), 'B': set()}
    for n in common:
        i0, i2 = r0['n'].index(n), r2['n'].index(n)
        for key in ('chord', 'A', 'B'):
            v0, v2 = r0[key][i0], r2[key][i2]
            if v0 % L == 0: pair[key].append(None); continue
            r = v2 * pow(v0, -1, m5) % m5; pair[key].append(digits(r, L, 5))
            if v2 % L: turns[key].add((dlog[v2 % L] - dlog[v0 % L]) % (L - 1))
    B_out[b] = pair
    say('37, levels %d → %d (%d values of n): different turns on the cycle: whole tower %d, twin A %d (%s), twin B %d (%s); '
        'the constants turn by %d and %d' % (b, b + 2, len(common), len(turns['chord']), len(turns['A']), sorted(turns['A']),
                                             len(turns['B']), sorted(turns['B']), dlog[15], dlog[34]))
OUT['A'] = A_out; OUT['B'] = B_out
# ---------------- C: 73's braid ----------------
L = 73; m = L
steps = open(os.path.join(HERE, 'partitions_flint_family_steps_73.txt'), encoding='utf-8').read()
CA = [int(u) % L for u in re.search(r"the constants' polynomial \(lowest first\) on lane ℓ: \[([^\]]*)\]", steps).group(1).split(',')]
CL = [int(u) % L for u in re.search(r"their polynomial \(lowest first\) on lane ℓ: \[([^\]]*)\]", steps).group(1).split(',')]


def fmul(a, b):                                            # the field with 73³ elements: F_73[x]/(CA)
    pr = [0] * 5
    for i in range(3):
        for j in range(3): pr[i + j] += a[i] * b[j]
    for tdeg in (4, 3):
        cc = pr[tdeg] % L
        if cc:
            for i in range(3): pr[tdeg - 3 + i] -= cc * CA[i]
        pr[tdeg] = 0
    return [pr[i] % L for i in range(3)]


def fpow(a, e):
    r = [1, 0, 0]
    while e:
        if e & 1: r = fmul(r, a)
        a = fmul(a, a); e >>= 1
    return r


def fadd(a, b): return [(u + v) % L for u, v in zip(a, b)]
def fsc(k, a): return [(k * u) % L for u in a]
ORD = L ** 3 - 1; X = [0, 1, 0]
k = L - 3; n0 = k + 1
Kz, d = F.class_basis_exact(k, leg(2, L), leg(3, L))
Xc, Rc = F.normalized_basis(nmod_mat(n0, d, [int(Kz[i, j]) % L for i in range(n0) for j in range(d)], L), d, L)
Hm = F.hecke_matrices(k, L, [('T', 5), ('T', L)])
T5 = F.restrict(Hm[('T', 5)], Xc, Rc, d); TL = F.restrict(Hm[('T', L)], Xc, Rc, d)
P = [F.identity(d, L), TL, TL * TL]
rows = [[int(P[i][a, b]) for i in range(3)] + [int(T5[a, b])] for a in range(3) for b in range(3)]
Rm = nmod_mat(9, 3, [r[i] for r in rows for i in range(3)], L); rhs = nmod_mat(9, 1, [r[3] for r in rows], L)
sol = (Rm.transpose() * Rm).solve(Rm.transpose() * rhs); Rr = [int(sol[i, 0]) for i in range(3)]
assert all(sum(Rr[i] * r[i] for i in range(3)) % L == r[3] for r in rows)
xa = fsc(leg(3, L), X); Rx = fadd(fadd([Rr[0], 0, 0], fsc(Rr[1], xa)), fsc(Rr[2], fmul(xa, xa))); lam1 = fsc(leg(3, 5), Rx)
def ev(poly, y):
    r = [0, 0, 0]
    for co in reversed(poly): r = fadd(fmul(r, y), [co, 0, 0])
    return r
assert ev(CL, lam1) == [0, 0, 0]
q = [None] * 3; q[2] = [1, 0, 0]
for i in (2, 1): q[i - 1] = fadd([CL[i], 0, 0], fmul(lam1, q[i]))
KAP = {}
for e5 in (1, -1):
    kap = e5 * leg(-3, 5) * pow(5, -2, L) % L; q3 = pow(5, -3, L); out = []
    for j in range(3):
        y = [1 if jj == j else 0 for jj in range(3)]; s = []
        for i in range(3):
            s.append(y[0] % L)
            if i < 2: y = [(y[1] + kap * y[0]) % L] + [(y[t + 1] + q3 * y[t - 1]) % L for t in range(1, len(y) - 1)]
        kj = [0, 0, 0]
        for i in range(3): kj = fadd(kj, fsc(s[i], q[i]))
        out.append(kj)
    KAP[e5] = out
# a generator of the field's cycle, and discrete logarithms (baby steps, giant steps)
fac = [2, 3, 1801]
gen = next(gg for gg in ([a0, a1, 1] for a1 in range(L) for a0 in range(L)) if all(fpow(gg, ORD // p) != [1, 0, 0] for p in fac))
Msq = 624; baby = {}; cur = [1, 0, 0]
for j in range(Msq): baby.setdefault(tuple(cur), j); cur = fmul(cur, gen)
giant = fpow(fpow(gen, Msq), ORD - 1)


def flog(a):
    cur = a
    for i in range(ORD // Msq + 2):
        if tuple(cur) in baby: return (i * Msq + baby[tuple(cur)]) % ORD
        cur = fmul(cur, giant)
    raise ValueError


C_out = {'ell': 73, 'order': ORD, 'levels': {}, 'theta': [flog(X), flog(fpow(X, L)), flog(fpow(X, L * L))]}
for b in (1, 3):
    t = time.time(); row = {'n': [], 's1': [], 's2': [], 's3': [], 'chord': []}
    for n in range(16):
        M = L ** b * (24 * n + 23)
        if M % 5 == 0: continue
        Kj = KAP[leg(M, 5)]; G1 = [0, 0, 0]
        for j in range(3): G1 = fadd(G1, fsc(c(25 ** j * M) % L, Kj[j]))
        row['n'].append(n); row['chord'].append(c(M) % L)
        if G1 == [0, 0, 0]: row['s1'].append(None); row['s2'].append(None); row['s3'].append(None); continue
        row['s1'].append(flog(G1)); row['s2'].append(flog(fpow(G1, L))); row['s3'].append(flog(fpow(G1, L * L)))
    C_out['levels'][b] = row
    say('73, level %d: %d values of n (%.0f s)' % (b, len(row['n']), time.time() - t))
ok = sum(1 for i, n in enumerate(C_out['levels'][3]['n']) if n in C_out['levels'][1]['n'] and C_out['levels'][3]['s1'][i] is not None
         and C_out['levels'][1]['s1'][C_out['levels'][1]['n'].index(n)] is not None
         and (C_out['levels'][3]['s1'][i] - C_out['levels'][1]['s1'][C_out['levels'][1]['n'].index(n)]) % ORD == C_out['theta'][0])
tot = sum(1 for i, n in enumerate(C_out['levels'][3]['n']) if C_out['levels'][3]['s1'][i] is not None)
say('73: strand 1 turns by θ = %d (of %d) from level 1 to 3 for %d of %d values of n; the strands\' angles θ, 73θ, 73²θ = %s'
    % (C_out['theta'][0], ORD, ok, tot, C_out['theta']))
g73 = next(gg for gg in range(2, 73) if all(pow(gg, 72 // p, 73) != 1 for p in (2, 3)))
dl73 = {pow(g73, i, 73): i for i in range(72)}; C_out['root_small'] = g73
C_out['chord_angle'] = {b: [None if v == 0 else dl73[v] for v in C_out['levels'][b]['chord']] for b in (1, 3)}
OUT['C'] = C_out
# ---------------- D: silences ----------------
QMAX = 10 ** 6; D_out = {}
for L, notes in ((37, (('A', 36), ('B', 5))), (97, (('g', 58),))):
    t = time.time(); S, Xc, Rc, d = class_series(L, QMAX + 1)
    primes = [int(Q) for Q in primerange(5, QMAX + 1) if Q != L]
    res = {}
    for name, a5 in notes:
        aQ = eigenform_coeffs(L, S, d, Xc, Rc, a5, primes)
        res[name] = {Q: silences(L, aQ[Q], Q) for Q in primes}
    D_out[L] = {name: {kind: [Q for Q in primes if res[name][Q] == kind] for kind in 'pmo'} for name in res}
    D_out[L]['primes'] = len(primes)
    if L == 37:
        chord = {kind: [Q for Q in primes if res['A'][Q] == kind and res['B'][Q] == kind] for kind in 'pmo'}
        old = open(os.path.join(HERE, 'partitions_flint_lucky_37_43_41_47_53_59_61_13_17_19_23_29_31.txt'), encoding='utf-8').read()
        blk = old.split('ℓ = 43')[0]
        oldl = [[int(x) for x in g.split(',') if x.strip()] for g in re.findall(r'      all: \[([^\]]*)\]', blk)]
        same = [chord['p'] == oldl[0], chord['m'] == oldl[1], chord['o'] == oldl[2]]
        D_out[37]['chord'] = chord
        say('37: twin A silent at %d, twin B at %d of %d primes (each about 3 in 37 PREDICTED: %.0f); both at once with one kind: '
            '%d, %d, %d (class +1, class −1, Ono) — the lists of the nineteenth time: %s' % (
                sum(len(D_out[37]['A'][k]) for k in 'pmo'), sum(len(D_out[37]['B'][k]) for k in 'pmo'), len(primes),
                3 * len(primes) / 37, len(chord['p']), len(chord['m']), len(chord['o']), 'the same, all three' if all(same) else 'DIFFERENT %s' % same))
    else:
        say('97: the single form g is silent at %d of %d primes (about 3 in 97 PREDICTED: %.0f); the first: %s' % (
            sum(len(D_out[97]['g'][k]) for k in 'pmo'), len(primes), 3 * len(primes) / 97,
            sorted(D_out[97]['g']['p'] + D_out[97]['g']['m'] + D_out[97]['g']['o'])[:8]))
    say('   (%.0f s)' % (time.time() - t))
OUT['D'] = D_out
json.dump(OUT, open(os.path.join(HERE, 'chords_page_data.json'), 'w'), separators=(',', ':'))
say('wrote chords_page_data.json; the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
