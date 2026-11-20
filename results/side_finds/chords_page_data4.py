# A dig, 13 Nov 2026: THE PAGES' DATA FOR 67, 71, 79, 83.
# From partitions_flint_family_steps_ℓ.txt (the constants' polynomial on lane ℓ), with the same methods as
# chords_page_data2.py (S: every note's turn and order; the loop), braid_angles.py (the strands' turns, their sum, the
# product of the constants) and braid_polar_data.py (size and direction steps; for pairs every number of the field).
# Writes chords_page_data4.json. COUNTED (orders, loops, identities), DERIVED (logarithms).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/chords_page_data4.py [--keep]
import os, sys, re, json, math, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
from flint import nmod_poly
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'chords_page_data4.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


class Fq:
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
    def primes_of(self, n):
        out, p = [], 2
        while p * p <= n:
            if n % p == 0:
                out.append(p)
                while n % p == 0: n //= p
            p += 1
        if n > 1: out.append(n)
        return out
    def gen(self, avoid=(), reverse=False):
        fac = self.primes_of(self.N)
        for co in itertools.product(range(self.L), repeat=self.e):
            g = tuple(reversed(co)) if reverse else tuple(co)
            if g == tuple([0] * self.e) or g in avoid: continue
            if all(self.pw(g, self.N // p) != self.one() for p in fac): return g
    def log(self, a, g):
        m = int(self.N ** 0.5) + 1; baby = {}; cur = self.one()
        for j in range(m): baby.setdefault(cur, j); cur = self.mul(cur, g)
        giant = self.pw(self.pw(g, m), self.N - 1); cur = a
        for i in range(m + 1):
            if cur in baby: return (i * m + baby[cur]) % self.N
            cur = self.mul(cur, giant)
        raise ValueError


OUT = {'S': {}, 'angles': {}, 'polar': {}}
for L in (67, 71, 79, 83):
    txt = open(os.path.join(HERE, 'partitions_flint_family_steps_%d.txt' % L), encoding='utf-8').read()
    CA = [int(u) % L for u in re.search(r"the constants' polynomial \(lowest first\) on lane ℓ: \[([^\]]*)\]", txt).group(1).split(',')]
    notes = []
    for fp, ex in nmod_poly(CA, L).factor()[1]:
        A = [int(u) for u in fp.coeffs()]; K = Fq(A, L)
        x = tuple([(-A[0]) % L]) if K.e == 1 else tuple([0, 1] + [0] * (K.e - 2))
        g = K.gen()                                            # as chords_page_data2.py (S)
        th = K.log(x, g); order = K.N // math.gcd(th, K.N)
        notes.append({'e': K.e, 'N': K.N, 'theta': th, 'order': order, 'first': (-A[0]) % L if K.e == 1 else None, 'factor': A})
        if K.e >= 2:
            # the strands' angles (as braid_angles.py)
            conj = [K.pw(x, L ** i) for i in range(K.e)]
            ga = K.gen(avoid=set(conj), reverse=True)
            tq = [K.log(c, ga) for c in conj]
            prod = K.one()
            for c in conj: prod = K.mul(prod, c)
            Nx = prod[0]; Ng = K.pw(ga, K.N // (L - 1))[0]
            logNx = next(i for i in range(L - 1) if pow(Ng, i, L) == Nx); Ssum = sum(tq) % K.N; step = K.N // (L - 1)
            ratio_ok = all((tq[(i + 1) % K.e] - L * tq[i]) % K.N == 0 for i in range(K.e))
            sum_ok = Ssum % step == 0 and Ssum // step == logNx
            OUT['angles'][L] = {'e': K.e, 'turns': tq, 'N': K.N, 'sum': Ssum, 'steps': Ssum // step, 'norm': Nx, 'square': leg(Nx, L) == 1}
            say('ℓ = %d, a braid of %d: turns %s of %s; ratio ×ℓ: %s; sum = %d steps of 360°/%d = the angle of the product %d (%s): %s' % (
                L, K.e, tq, format(K.N, ','), 'yes' if ratio_ok else 'NO', Ssum // step, L - 1, Nx, 'a square' if leg(Nx, L) == 1 else 'not a square', 'yes' if sum_ok else 'NO'))
            # size and direction (as braid_polar_data.py)
            Kdir = K.N // (L - 1); D = K.pw(x, L - 1)
            ordD = next(dd for dd in range(1, Kdir + 1) if Kdir % dd == 0 and K.pw(D, dd) == K.one())
            E = D if ordD == Kdir else K.pw(ga, L - 1)
            tab = {}; cur = K.one()
            for i in range(Kdir): tab[cur] = i; cur = K.mul(cur, E)
            size = lambda v: K.pw(v, Kdir)[0]; direc = lambda v: tab[K.pw(v, L - 1)]
            info = {'e': K.e, 'Kdir': Kdir, 'Na': size(x), 'step': direc(x), 'own': ordD == Kdir, 'orderD': ordD}
            if K.e == 2:
                info['all'] = [[size(v), direc(v)] for v in itertools.product(range(L), repeat=2) if v != (0, 0)]
            OUT['polar'][L] = info
            say('   size ×%d; direction %d of %s per two levels%s' % (info['Na'], info['step'], format(Kdir, ','), ' (its own clock)' if info['own'] else ''))
    loop = 1
    for nt in notes: loop = loop * nt['order'] // math.gcd(loop, nt['order'])
    OUT['S'][L] = {'notes': notes, 'loop': loop}
    say('ℓ = %d: notes %s; the loop %s two-level steps' % (L, ['constant %d (%d points)' % (n['first'], n['order']) if n['e'] == 1 else 'a braid of %d (%s points)' % (n['e'], format(n['order'], ',')) for n in notes], format(loop, ',')))
json.dump(OUT, open(os.path.join(HERE, 'chords_page_data4.json'), 'w'), separators=(',', ':'))
say('wrote chords_page_data4.json')
