# A dig, 13 Nov 2026: 607, THE EXCEPTION, AS A PICTURE (README "607 on a
# page"). The class of 607 at weight 604 (as partitions_flint_607_pair.py: forms6_flint.class_padic modulo 607^3), with
# T_607 and T₅ on it. Written:
#   - the constants' polynomial, (3|ℓ)·T_ℓ on lane 607, and its factors: the notes. Each note's turn and order as in
#     chords_page_data4.py. For the braid of 22: its order, its size step N(α) and the order of its direction α^606, from
#     the factors of 607^22 − 1 (its logarithm in the field with 607^22 elements is not computed: the largest prime factor
#     is about 4·10^15);
#   - the pair whose T_607 numbers meet on lane 607, in size and direction (field with 607² elements; size: the norm, a
#     number mod 607; direction: v^606 on the circle of 608 places, counted from E = G^606, G the first generator found):
#     at 5 (T₅), at 607 on the first digit (T_607 ≡ 367) and one digit deeper (C = (T_607 − 367)/607). Each form is
#     followed by its own eigenvector, so form A at 5 and form A one digit deeper are the same form.
# Writes chords_page_data607.json. COUNTED (factors, orders, eigenvalues), DERIVED (logarithms on the circles).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/chords_page_data607.py [--keep]
import os, sys, json, math, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat, nmod_poly, fmpz
import forms6_flint as F
flint.ctx.threads = 4
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'chords_page_data607.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


class F2:
    """The field with ℓ² elements as F_ℓ[t]/(t² + c1·t + c0); elements (a, b) = a + b·t."""
    def __init__(self, c0, c1, L): self.c0, self.c1, self.L = c0 % L, c1 % L, L
    def mul(self, x, y):
        L = self.L; a, b = x; c, d = y; bd = b * d                       # t² = −c1·t − c0
        return ((a * c - bd * self.c0) % L, (a * d + b * c - bd * self.c1) % L)
    def pw(self, x, n):
        r = (1, 0)
        while n:
            if n & 1: r = self.mul(r, x)
            x = self.mul(x, x); n >>= 1
        return r
    def inv(self, x): return self.pw(x, self.L ** 2 - 2)


t0 = time.time(); L = 607; k = L - 3; N = 3; m = L ** N; s3 = leg(3, L)
H = F.hecke_matrices(k, m, [('U', 2), ('U', 3), ('T', L), ('T', 5)])
W, RW, r, Wn, Rn, d, Np = F.class_padic(H, k, L, N, leg(2, L), leg(3, L))
mp = L ** Np
TL = F.restrict(F.mat_mod(k + 1, k + 1, [int(v) % mp for v in H[('T', L)].entries()], mp), Wn, Rn, d)
T5 = F.restrict(F.mat_mod(k + 1, k + 1, [int(v) % mp for v in H[('T', 5)].entries()], mp), Wn, Rn, d)
say('ℓ = 607, weight 604: the class has %d forms, known to 607^%d; (3|607) = %+d  (%.0f s)' % (d, Np, s3, time.time() - t0))

# ---- the notes on lane 607: the constants' polynomial (3|ℓ)·T_ℓ ----
TLb = F.reduce_mat(TL, L)
CA = (TLb * s3).charpoly()
fac = CA.factor()[1]
g0 = next(a for a in range(2, L) if all(pow(a, (L - 1) // q, L) != 1 for q in (2, 3, 101)))   # 3
logt = {}; x_ = 1
for i in range(L - 1): logt[x_] = i; x_ = x_ * g0 % L
ordp = lambda a: (L - 1) // math.gcd(logt[a % L], L - 1)
# the factors of 607^22 − 1 = 606 · 608 · Φ11(607) · Φ22(607)
parts = [L - 1, L + 1, (L ** 11 - 1) // (L - 1), (L ** 11 + 1) // (L + 1)]
P22 = {}
for q in parts:
    for p, e in fmpz(q).factor(): P22[int(p)] = P22.get(int(p), 0) + int(e)
assert math.prod(p ** e for p, e in P22.items()) == L ** 22 - 1 and all(fmpz(p).is_prime() for p in P22)
say('607^22 − 1 = %s' % ' · '.join('%d^%d' % (p, e) if e > 1 else str(p) for p, e in sorted(P22.items())))


def order_in(x, f, n, primes):
    o = n; one = nmod_poly([1], L)
    for p in primes:
        while o % p == 0 and x.pow_mod(o // p, f) == one: o //= p
    return o


notes = []; braid = None
for f, ex in fac:
    A = [int(u) for u in f.coeffs()]; e = len(A) - 1
    for _ in range(ex):
        if e == 1:
            c = (-A[0]) % L
            notes.append({'e': 1, 'N': L - 1, 'theta': logt[c], 'order': ordp(c), 'first': c, 'factor': A, 'mult': ex})
        else:
            n22 = L ** e - 1; prs = [p for p in P22] if e == 22 else None
            if e != 22: raise SystemExit('unexpected braid degree %d' % e)
            fx = nmod_poly(A, L); x = nmod_poly([0, 1], L)
            o = order_in(x, fx, n22, prs)
            Kdir = n22 // (L - 1); Dx = x.pow_mod(L - 1, fx)
            oD = order_in(Dx, fx, Kdir, [p for p in prs if Kdir % p == 0])
            Na = int(x.pow_mod(Kdir, fx).coeffs()[0]) % L
            assert Na == A[0] % L                                  # the norm is the product of the 22 constants
            notes.append({'e': e, 'N': str(n22), 'theta': None, 'order': str(o), 'first': None, 'factor': A, 'mult': ex})
            braid = {'e': e, 'N': str(n22), 'order': str(o), 'Kdir': str(Kdir), 'Na': Na, 'Na_order': ordp(Na), 'Na_log': logt[Na],
                     'orderD': str(oD), 'own': oD == Kdir, 'square': leg(Na, L) == 1}
for nt in notes:
    say('   note: %s' % ('constant %d (turn %d of 606, %d points)%s' % (nt['first'], nt['theta'], nt['order'], '  ×%d (the pair)' % nt['mult'] if nt['mult'] > 1 else '')
                         if nt['e'] == 1 else 'a braid of %d: order %s of %s' % (nt['e'], format(int(nt['order']), ','), format(int(nt['N']), ','))))
say('   the braid of 22: size ×%d per two levels (the product of its 22 constants; %s; %d points on the 607-hour clock); '
    'direction: order %s of the circle of %s places (%s)' % (braid['Na'], 'a square' if braid['square'] else 'not a square', braid['Na_order'],
                                                               format(int(braid['orderD']), ','), format(int(braid['Kdir']), ','), 'it generates it: its own clock' if braid['own'] else 'a smaller circle'))
loop = 1
for nt in notes: o = int(nt['order']); loop = loop * o // math.gcd(loop, o)
say('   the whole tower of 607 comes back after %s two-level steps on lane 607 (the least common multiple)' % format(loop, ','))

# ---- the pair: the piece where T_607 ≡ 367 twice ----
rep = [(f, e) for f, e in TLb.charpoly().factor()[1] if e > 1]
assert len(rep) == 1 and rep[0][0].degree() == 1
a = (-int(rep[0][0].coeffs()[0])) % L
Y, R, s = F.kernel_part(TL - F.identity(d, mp) * a, L, Np)
B = F.restrict(TL, Y, R, s); B5 = F.restrict(T5, Y, R, s)
Bl = [[int(B[i, j]) for j in range(s)] for i in range(s)]; B5l = [[int(B5[i, j]) % L for j in range(s)] for i in range(s)]
assert s == 2 and all((Bl[i][j] - (a if i == j else 0)) % L == 0 for i in range(2) for j in range(2))
C = [[((Bl[i][j] - (a if i == j else 0)) // L) % L for j in range(2)] for i in range(2)]
chi5 = [(B5l[0][0] * B5l[1][1] - B5l[0][1] * B5l[1][0]) % L, (-(B5l[0][0] + B5l[1][1])) % L, 1]
chiC = [(C[0][0] * C[1][1] - C[0][1] * C[1][0]) % L, (-(C[0][0] + C[1][1])) % L, 1]
say('   the pair: T_607 ≡ %d (×2) on lane 607, so both constants are (3|607)·%d = %d; T₅ there: %s (irreducible: %s); C: %s (irreducible: %s)' % (
    a, a, s3 * a % L, chi5, nmod_poly(chi5, L).factor()[1][0][0].degree() == 2, chiC, nmod_poly(chiC, L).factor()[1][0][0].degree() == 2))
K2 = F2(chi5[0], chi5[1], L)                                       # F_(607²) = F_607[t]/(T₅'s polynomial): γ_A = t
gamA = (0, 1); gamB = K2.pw(gamA, L)
# form A's eigenvector of T₅ for t: v = (b12, t − b11) (b12 ≠ 0), then C·v = δ_A·v
assert B5l[0][1] != 0
v = ((B5l[0][1], 0), ((-B5l[0][0]) % L, 1))
Cv0 = ((C[0][0] * v[0][0] + C[0][1] * v[1][0]) % L, (C[0][0] * v[0][1] + C[0][1] * v[1][1]) % L)
Cv1 = ((C[1][0] * v[0][0] + C[1][1] * v[1][0]) % L, (C[1][0] * v[0][1] + C[1][1] * v[1][1]) % L)
delA = K2.mul(Cv0, K2.inv(v[0]))
assert K2.mul(delA, v[1]) == Cv1                                   # v is an eigenvector of C too
delB = K2.pw(delA, L)
# a generator of F_(607²)^* and the direction circle of 608 places
NN = L * L - 1; prs2 = [2, 3, 19, 101]
G = next(g for g in itertools.product(range(L), repeat=2) if g != (0, 0) and all(K2.pw(g, NN // p) != (1, 0) for p in prs2))
E = K2.pw(G, L - 1); dtab = {}; cur = (1, 0)
for i in range(L + 1): dtab[cur] = i; cur = K2.mul(cur, E)
size = lambda u: K2.pw(u, L + 1)[0]
direc = lambda u: dtab[K2.pw(u, L - 1)]
pair = {}
for name, (uA, uB), what in (('at5', (gamA, gamB), 'T₅'), ('at607', ((a, 0), (a, 0)), 'T_607, first digit'), ('deeper', (delA, delB), 'C, one digit deeper')):
    pair[name] = {'A': [size(uA), direc(uA)], 'B': [size(uB), direc(uB)], 'A_el': list(uA), 'B_el': list(uB)}
    say('   %s: form A size %d, direction %d of 608; form B size %d, direction %d of 608' % (what, size(uA), direc(uA), size(uB), direc(uB)))
OUT = {'L': L, 'k': k, 'forms': d, 'sign3': s3, 'notes': notes, 'braid': braid, 'loop': str(loop), 'pair_T': a, 'pair_const': s3 * a % L,
       'chi5': chi5, 'chiC': chiC, 'pair': pair, 'generator': g0, 'P22': {str(p): e for p, e in P22.items()}}
if KEEP: json.dump(OUT, open(os.path.join(HERE, 'chords_page_data607.json'), 'x'), separators=(',', ':'))
say('wrote chords_page_data607.json' if KEEP else '(not kept)'); say('seconds: %.0f' % (time.time() - t0))
