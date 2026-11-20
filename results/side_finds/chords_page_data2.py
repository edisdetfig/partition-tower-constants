# A dig, 13 Nov 2026: MORE DATA FOR THE PAGE "Towers as Chords": THE FULL FORMS. Exact partition numbers (FLINT), and the notes' constants on lane ℓ from the family runs. Writes
# chords_page_data2.json. COUNTED (values, orders, loops), DERIVED (angles, turns).
#  A2  37: levels 1, 3, 5 for n = 0 … 359 (5 ∤ M): the whole tower's and each twin's angle on 37's cycle (root 2), and the
#      multipliers' digits modulo 37⁵ between levels 1 → 3 and 3 → 5 (as in chords_page_data.py).
#  C2  73: levels 1, 3 for n = 0 … 299: each strand's angle on the cycle of the field with 73³ elements, and the whole tower's
#      on 73's cycle.
#  S   every note's constant on lane ℓ for ℓ = 13 … 97: a root x of each factor of the constants' polynomial on lane ℓ
#      (degree e = 1, 2 or 3), its turn θ = log(x) on the cycle of the field with ℓ^e elements and its order
#      (ℓ^e − 1)/gcd(θ, ℓ^e − 1): following the note two levels at a time traces a star with that many points. The
#      tower's loop on lane ℓ is the least common multiple of its notes' orders (checked against partitions_loops.txt).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/chords_page_data2.py [--keep]
import os, sys, re, json, time, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly, nmod_mat
import forms6_flint as F
flint.ctx.threads = 4
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'chords_page_data2.txt'), 'x', encoding='utf-8') if KEEP else None
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


class Fq:
    """The field F_ℓ[x]/(A), A monic irreducible on lane ℓ of degree e (coefficients lowest first)."""
    def __init__(self, A, L):
        self.A = [a % L for a in A]; self.e = len(A) - 1; self.L = L; self.N = L ** self.e - 1
    def mul(self, a, b):
        e, L, A = self.e, self.L, self.A; pr = [0] * (2 * e - 1)
        for i in range(e):
            if a[i]:
                for j in range(e): pr[i + j] += a[i] * b[j]
        for t in range(2 * e - 2, e - 1, -1):
            cc = pr[t] % L
            if cc:
                for i in range(e): pr[t - e + i] -= cc * A[i]
            pr[t] = 0
        return tuple(pr[i] % L for i in range(e))
    def one(self): return tuple([1] + [0] * (self.e - 1))
    def pw(self, a, n):
        r = self.one()
        while n:
            if n & 1: r = self.mul(r, a)
            a = self.mul(a, a); n >>= 1
        return r
    def gen(self):
        fac = [p for p in range(2, int(self.N ** 0.5) + 2) if self.N % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1))]
        rest = self.N
        for p in fac:
            while rest % p == 0: rest //= p
        if rest > 1: fac.append(rest)
        import itertools
        for co in itertools.product(range(self.L), repeat=self.e):
            g = tuple(co)
            if g == tuple([0] * self.e): continue
            if all(self.pw(g, self.N // p) != self.one() for p in fac): return g
    def log(self, a, g):
        m = int(self.N ** 0.5) + 1; baby = {}; cur = self.one()
        for j in range(m): baby.setdefault(cur, j); cur = self.mul(cur, g)
        giant = self.pw(self.pw(g, m), self.N - 1); cur = a
        for i in range(m + 1):
            if cur in baby: return (i * m + baby[cur]) % self.N
            cur = self.mul(cur, giant)
        raise ValueError('no logarithm')


t0 = time.time()
# ---------------- A2: 37, many n ----------------
L = 37; m5 = L ** 5; dlog = {pow(2, i, L): i for i in range(L - 1)}
steps = open(os.path.join(HERE, 'partitions_flint_family_steps_37.txt'), encoding='utf-8').read()
lams = [[int(x) for x in g.split(',')] for g in re.findall(r'a root on lane ℓ, digits lowest first: \[([^\]]*)\]', steps)]
num = lambda dg: sum(x * L ** i for i, x in enumerate(dg)) % m5
lamA = num(next(l for l in lams if l[0] == 1)); lamB = num(next(l for l in lams if l[0] == 32))   # twin A keeps λ(5) 1…, B 32…
ns = [n for n in range(360) if (24 * n + 11) % 5]
A2 = {'n': ns, 'levels': {}, 'mult': {}}; vals = {}
for b in (1, 3, 5):
    t = time.time(); row = {'chord': [], 'A': [], 'B': []}; raw = {'chord': [], 'A': [], 'B': []}
    for n in ns:
        M = L ** b * (24 * n + 11); cM, c25 = c(M) % m5, c(25 * M) % m5
        kap = leg(M, 5) * leg(-3, 5) * pow(5, -2, m5) % m5
        gA = (c25 + (kap - lamB) * cM) % m5; gB = (c25 + (kap - lamA) * cM) % m5
        for key, v in (('chord', cM), ('A', gA), ('B', gB)):
            raw[key].append(v); row[key].append(None if v % L == 0 else dlog[v % L])
    A2['levels'][b] = row; vals[b] = raw
    say('37, level %d: %d values of n (%.0f s)' % (b, len(ns), time.time() - t))
for b in (1, 3):
    mm = {}
    for key in ('chord', 'A', 'B'):
        out = []
        for v0, v2 in zip(vals[b][key], vals[b + 2][key]):
            if v0 % L == 0: out.append(None); continue
            r = v2 * pow(v0, -1, m5) % m5; out.append([(r // L ** i) % L for i in range(5)])
        mm[key] = out
    A2['mult'][b] = mm
turnsA = {(a2 - a0) % 36 for a0, a2 in zip(A2['levels'][1]['A'], A2['levels'][3]['A']) if a0 is not None and a2 is not None}
turnsB = {(a2 - a0) % 36 for a0, a2 in zip(A2['levels'][1]['B'], A2['levels'][3]['B']) if a0 is not None and a2 is not None}
turnsW = {(a2 - a0) % 36 for a0, a2 in zip(A2['levels'][1]['chord'], A2['levels'][3]['chord']) if a0 is not None and a2 is not None}
say('37, levels 1 → 3 over %d values of n: turns of twin A %s, twin B %s, the whole tower %d different' % (len(ns), sorted(turnsA), sorted(turnsB), len(turnsW)))
OUT['A2'] = A2
# ---------------- C2: 73's strands, many n ----------------
L = 73
steps = open(os.path.join(HERE, 'partitions_flint_family_steps_73.txt'), encoding='utf-8').read()
CA = [int(u) % L for u in re.search(r"the constants' polynomial \(lowest first\) on lane ℓ: \[([^\]]*)\]", steps).group(1).split(',')]
CL = [int(u) % L for u in re.search(r"their polynomial \(lowest first\) on lane ℓ: \[([^\]]*)\]", steps).group(1).split(',')]
K = Fq(CA, L); X = (0, 1, 0)
fadd = lambda a, b: tuple((u + v) % L for u, v in zip(a, b)); fsc = lambda k, a: tuple((k * u) % L for u in a)
k = L - 3; n0 = k + 1
Kz, d = F.class_basis_exact(k, leg(2, L), leg(3, L))
Xc, Rc = F.normalized_basis(nmod_mat(n0, d, [int(Kz[i, j]) % L for i in range(n0) for j in range(d)], L), d, L)
Hm = F.hecke_matrices(k, L, [('T', 5), ('T', L)])
T5 = F.restrict(Hm[('T', 5)], Xc, Rc, d); TL = F.restrict(Hm[('T', L)], Xc, Rc, d)
P = [F.identity(d, L), TL, TL * TL]
rows = [[int(P[i][a, b]) for i in range(3)] + [int(T5[a, b])] for a in range(3) for b in range(3)]
Rm = nmod_mat(9, 3, [r[i] for r in rows for i in range(3)], L); rhs = nmod_mat(9, 1, [r[3] for r in rows], L)
sol = (Rm.transpose() * Rm).solve(Rm.transpose() * rhs); Rr = [int(sol[i, 0]) for i in range(3)]
xa = fsc(leg(3, L), X); lam1 = fsc(leg(3, 5), fadd(fadd((Rr[0], 0, 0), fsc(Rr[1], xa)), fsc(Rr[2], K.mul(xa, xa))))
ev = lambda poly, y: __import__('functools').reduce(lambda r, co: fadd(K.mul(r, y), (co, 0, 0)), reversed(poly), (0, 0, 0))
assert ev(CL, lam1) == (0, 0, 0)
q = [None] * 3; q[2] = (1, 0, 0)
for i in (2, 1): q[i - 1] = fadd((CL[i], 0, 0), K.mul(lam1, q[i]))
KAP = {}
for e5 in (1, -1):
    kap = e5 * leg(-3, 5) * pow(5, -2, L) % L; q3 = pow(5, -3, L); out = []
    for j in range(3):
        y = [1 if jj == j else 0 for jj in range(3)]; s = []
        for i in range(3):
            s.append(y[0] % L)
            if i < 2: y = [(y[1] + kap * y[0]) % L] + [(y[t + 1] + q3 * y[t - 1]) % L for t in range(1, len(y) - 1)]
        kj = (0, 0, 0)
        for i in range(3): kj = fadd(kj, fsc(s[i], q[i]))
        out.append(kj)
    KAP[e5] = out
g73 = K.gen(); gsmall = next(gg for gg in range(2, L) if all(pow(gg, (L - 1) // p, L) != 1 for p in (2, 3)))
dl73 = {pow(gsmall, i, L): i for i in range(L - 1)}
C2 = {'n': [], 'levels': {}, 'theta': [K.log(X, g73), K.log(K.pw(X, L), g73), K.log(K.pw(X, L * L), g73)], 'order': K.N}
nsC = [n for n in range(300) if (L * (24 * n + 23)) % 5]
C2['n'] = nsC
for b in (1, 3):
    t = time.time(); row = {'s1': [], 's2': [], 's3': [], 'chord': []}
    for n in nsC:
        M = L ** b * (24 * n + 23); Kj = KAP[leg(M, 5)]; G1 = (0, 0, 0)
        for j in range(3): G1 = fadd(G1, fsc(c(25 ** j * M) % L, Kj[j]))
        cm = c(M) % L; row['chord'].append(None if cm == 0 else dl73[cm])
        if G1 == (0, 0, 0):
            for s_ in ('s1', 's2', 's3'): row[s_].append(None)
            continue
        row['s1'].append(K.log(G1, g73)); row['s2'].append(K.log(K.pw(G1, L), g73)); row['s3'].append(K.log(K.pw(G1, L * L), g73))
    C2['levels'][b] = row
    say('73, level %d: %d values of n (%.0f s)' % (b, len(nsC), time.time() - t))
good = sum(1 for a, b2 in zip(C2['levels'][1]['s1'], C2['levels'][3]['s1']) if a is not None and b2 is not None and (b2 - a) % K.N == C2['theta'][0])
tot = sum(1 for a, b2 in zip(C2['levels'][1]['s1'], C2['levels'][3]['s1']) if a is not None and b2 is not None)
say('73: strand 1 turns by θ = %d for %d of %d values of n' % (C2['theta'][0], good, tot))
OUT['C2'] = C2
# ---------------- S: every note's star ----------------
loops_known = {13: 12, 17: 16, 19: 18, 23: 22, 29: 28, 31: 2, 37: 36, 41: 10, 43: 21, 47: 2208, 53: 26, 59: 870, 61: 3720}
S = {}
for L in (13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 73, 97):
    txt = open(os.path.join(HERE, 'partitions_flint_family_%d.txt' % L), encoding='utf-8').read()
    CA = [int(u) % L for u in re.search(r"the constants' polynomial \(lowest first\) on lane ℓ: \[([^\]]*)\]", txt).group(1).split(',')]
    notes = []
    for fp, ex in nmod_poly(CA, L).factor()[1]:
        A = [int(u) for u in fp.coeffs()]; Kf = Fq(A, L); g = Kf.gen()
        x = tuple([(-A[0]) % L]) if Kf.e == 1 else tuple([0, 1] + [0] * (Kf.e - 2))
        th = Kf.log(x, g); order = Kf.N // math.gcd(th, Kf.N)
        notes.append({'e': Kf.e, 'N': Kf.N, 'theta': th, 'order': order, 'first': (-A[0]) % L if Kf.e == 1 else None, 'factor': A})
    loop = 1
    for nt in notes: loop = loop * nt['order'] // math.gcd(loop, nt['order'])
    S[L] = {'notes': notes, 'loop': loop}
    say('ℓ = %d: notes %s; the loop (lcm of the orders) %d%s' % (L, ['%s: turn %d of %d, %d points' % ('constant %d' % nt['first'] if nt['e'] == 1 else 'a group of %d in the field of %d^%d' % (nt['e'], L, nt['e']), nt['theta'], nt['N'], nt['order']) for nt in notes],
                                                      loop, ('; partitions_loops.txt: %d, %s' % (loops_known[L], 'the same' if loops_known[L] == loop else 'DIFFERENT')) if L in loops_known else ''))
OUT['S'] = S
json.dump(OUT, open(os.path.join(HERE, 'chords_page_data2.json'), 'w'), separators=(',', ':'))
say('wrote chords_page_data2.json; the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
