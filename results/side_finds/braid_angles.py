# A dig, 13 Nov 2026: WHAT THE BRAIDS' ANGLES OBEY (# README "The braids' angles"; predicted there first). For each braided group (pairs at 47, 59, 61; triples at 73, 97):
# the field F = F_ℓ[x]/(A) of ℓ^e elements (A the group's factor of the constants' polynomial on lane ℓ,
# partitions_flint_family_ℓ.txt), a generator γ of its cycle, and each strand's turn θ_i = log_γ(x^(ℓ^i)), i = 0 … e − 1.
# PREDICTED (derived from Frobenius x ↦ x^ℓ and the norm N(x) = x^(1 + ℓ + … + ℓ^(e−1)) ∈ F_ℓ):
#   ratio:  θ_(i+1) ≡ ℓ·θ_i (mod ℓ^e − 1), every strand;
#   sum:    Σθ_i ≡ log_γ(N(x)), a multiple of (ℓ^e − 1)/(ℓ − 1): on the circle, a multiple of 360°/(ℓ − 1), and exactly the
#           angle of the product of the group's constants, N(x) = (−1)^e·A(0), on ℓ's own cycle walked by N(γ);
#   differences: θ_1 − θ_0 = log_γ(x^(ℓ−1)), a multiple of ℓ − 1 places: on the circle a multiple of 360°/((ℓ^e − 1)/(ℓ − 1)).
#   The sum is an odd number of steps exactly when N(x) is not a square mod ℓ (then the constants are not squares, and a
#   one-level law is impossible; README "The towers of 73 and 97").
# Also, for the twins of 37, 41, 43, 53 (two constants in F_ℓ): the exponent k with B ≡ A^k (mod ℓ), when the order of B
# divides that of A. COUNTED (the identities on the computed logarithms), DERIVED (the predictions).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/braid_angles.py [--keep]
import os, sys, re, math, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
from flint import nmod_poly
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'braid_angles.txt'), 'x', encoding='utf-8') if KEEP else None


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
    def gen(self, avoid=()):
        fac = self.primes_of(self.N)
        for co in itertools.product(range(self.L), repeat=self.e):
            g = tuple(reversed(co))
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


res = {}
say('The braids: each strand\'s turn θ_i = log_γ(x^(ℓ^i)) on the cycle of the field with ℓ^e elements (γ the first generator in the')
say('field\'s listing other than the group\'s own constants); angles in degrees of a full turn.')
for L in (47, 59, 61, 73, 97):
    txt = open(os.path.join(HERE, 'partitions_flint_family_%d.txt' % L), encoding='utf-8').read()
    CA = [int(u) % L for u in re.search(r"the constants' polynomial \(lowest first\) on lane ℓ: \[([^\]]*)\]", txt).group(1).split(',')]
    for fp, ex in nmod_poly(CA, L).factor()[1]:
        A = [int(u) for u in fp.coeffs()]
        if len(A) - 1 < 2: continue
        K = Fq(A, L); e = K.e; N = K.N; x = tuple([0, 1] + [0] * (e - 2))
        conj = [K.pw(x, L ** i) for i in range(e)]
        g = K.gen(avoid=set(conj))
        th = [K.log(cj, g) for cj in conj]
        ratio_ok = all((th[(i + 1) % e] - L * th[i]) % N == 0 for i in range(e))
        prod = K.one()
        for cj in conj: prod = K.mul(prod, cj)
        assert all(c == 0 for c in prod[1:]); Nx = prod[0]
        Ng = K.pw(g, N // (L - 1))[0]                           # N(γ) = γ^(1 + ℓ + …): a primitive root mod ℓ
        logNx = next(i for i in range(L - 1) if pow(Ng, i, L) == Nx)
        S = sum(th) % N; step = N // (L - 1)
        sum_ok = S % step == 0 and S // step == logNx
        diff = (th[1] - th[0]) % N; dstep = L - 1
        deg = [360.0 * t / N for t in th]
        say(); say('ℓ = %d, a braid of %d (field of %d^%d = %s elements; constants\' factor on lane ℓ %s)' % (L, e, L, e, format(N + 1, ','), A))
        say('   turns: %s of %s places = %s°' % (', '.join(format(t, ',') for t in th), format(N, ','), ', '.join('%.3f' % d for d in deg)))
        say('   ratio: each strand\'s turn ≡ %d × the one before (mod %s): %s' % (L, format(N, ','), 'as PREDICTED' if ratio_ok else 'NOT as predicted'))
        say('   sum: %s places = %d × %s, i.e. %.3f° = %d × %.4f° (a multiple of 360°/%d);  the product of the constants N(x) = %d (mod %d) = %d^%d on %d\'s cycle walked by N(γ) = %d: %s' % (
            format(S, ','), S // step if S % step == 0 else -1, format(step, ','), 360.0 * S / N, S // step if S % step == 0 else -1, 360.0 / (L - 1), L - 1,
            Nx, L, Ng, logNx, L, Ng, 'as PREDICTED' if sum_ok else 'NOT as predicted'))
        say('   N(x) = %d is %s mod %d; the sum is an %s number of steps (%d): %s' % (
            Nx, 'a square' if leg(Nx, L) == 1 else 'not a square', L, 'even' if logNx % 2 == 0 else 'odd', logNx,
            'as PREDICTED' if (leg(Nx, L) == 1) == (logNx % 2 == 0) else 'NOT as predicted'))
        say('   difference θ_1 − θ_0 = %s places = %d × %d: on the circle a multiple of 360°/%s = %.5f°: %s' % (
            format(diff, ','), diff // dstep if diff % dstep == 0 else -1, dstep, format(N // dstep, ','), 360.0 * dstep / N, 'as PREDICTED' if diff % dstep == 0 else 'NOT as predicted'))
        res[L] = {'e': e, 'theta': th, 'N': N, 'norm': Nx, 'logNorm': logNx, 'Ng': Ng}
# the twins: B as a power of A
say(); say('The twins (two constants on lane ℓ): B ≡ A^k (mod ℓ)')
for L in (37, 41, 43, 53):
    txt = open(os.path.join(HERE, 'partitions_flint_family_%d.txt' % L), encoding='utf-8').read()
    al = [int(g.split(',')[0]) for g in re.findall(r'a constant on lane ℓ, digits lowest first: \[([^\]]*)\]', txt)]
    a, b = al[0], al[1]
    oa = next(k for k in range(1, L) if pow(a, k, L) == 1); ob = next(k for k in range(1, L) if pow(b, k, L) == 1)
    if oa % ob == 0:
        k = next(k for k in range(oa) if pow(a, k, L) == b)
        say('   ℓ = %d: constants %d (order %d) and %d (order %d): %d ≡ %d^%d (mod %d)' % (L, a, oa, b, ob, b, a, k, L))
    elif ob % oa == 0:
        k = next(k for k in range(ob) if pow(b, k, L) == a)
        say('   ℓ = %d: constants %d (order %d) and %d (order %d): %d ≡ %d^%d (mod %d)' % (L, b, ob, a, oa, a, b, k, L))
    else:
        say('   ℓ = %d: constants %d (order %d) and %d (order %d): neither is a power of the other' % (L, a, oa, b, ob))
