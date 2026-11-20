# A dig, 12 Nov 2026: THE TOWER'S LOOPS. Two levels up multiply the tower's limit by its constants (STATEMENT.md, A1), so on
# lane ℓ^m the tower comes back to itself after c_m two-level steps, c_m the order of the constants (of the matrix they make,
# the companion matrix of their polynomial) modulo ℓ^m. Boylan–Webb printed c₃₇ = 36 and c₄₁ = 10 for m = 1
# (partitions_bw_check.txt agreed). Here every prime from 5 to 61, at every depth the constants reach.
# Write C^(c₁) = I + ℓ^v·E (E not 0 on lane ℓ): then c_m = c₁ for m ≤ v and c₁·ℓ^(m−v) beyond (ℓ odd; derived). v is the
# loop's "return depth": the tower comes back to itself to v digits after c₁ steps.
# PREDICTED before the run (derived, the constants treated as random ℓ-adic units):
#   - v = 1 as a rule; v ≥ 2 with chance 1/ℓ per constant (about one case in all expected among 5 to 61);
#   - one constant (ℓ = 13 to 31) is a primitive root on lane ℓ (c₁ = ℓ − 1) with chance φ(ℓ − 1)/(ℓ − 1).
# Part 2 checks the shortest loops on exact partition numbers (FLINT): 5 and 31, whose constants are −1 on lane ℓ (so two
# levels flip the sign and four come back), 11 (c₁ = 5: ten levels) and 7 (c₁ = 6: twelve levels). H′_b = p(ℓᵇn + δ_b),
# divided by the power of ℓ of Watson's and Atkin's laws for 5, 7, 11 (5ᵇ, 7^(⌊b/2⌋+1), 11ᵇ). PREDICTED: H′_(b+2c₁) ≡ H′_b
# and H′_(b+2) ≡ α·H′_b on lane ℓ, for b = 1, 2 and n = 0..9. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_loops.py [--keep]
import os, sys, re, glob, math, time
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
import flint
from flint import nmod_poly
flint.ctx.threads = 4
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_loops.txt'), 'x', encoding='utf-8') if KEEP else None
DIG = {5: [4, 2, 3, 1, 0, 2, 0, 3, 0, 1, 1, 0, 1, 4, 0, 1, 3, 2, 3, 0, 3, 1, 0, 4, 3, 0, 4, 1, 2, 4, 2, 3, 0],
       7: [5, 6, 1, 1, 2, 1, 3, 1, 0, 6, 5, 3, 1, 1], 11: [9, 3, 2, 0]}      # settled from partition numbers (partitions_constants.txt)


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def poly_of(L):
    if L in DIG:
        K = len(DIG[L]); a = sum(d * L ** j for j, d in enumerate(DIG[L])); return K, [(-a) % L ** K, 1], 'partitions_constants.txt'
    best = None
    for f in glob.glob(os.path.join(HERE, 'partitions_flint_family_%d*.txt' % L)):
        mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\]', open(f, encoding='utf-8').read())
        if mt and (best is None or int(mt.group(1)) > best[0]): best = (int(mt.group(1)), [int(v) for v in mt.group(2).split(',')], os.path.basename(f))
    return best


def matmul(A, B, m):
    return [[sum(A[i][l] * B[l][j] for l in range(len(B))) % m for j in range(len(B[0]))] for i in range(len(A))]


def matpow(A, e, m):
    d = len(A); R = [[int(i == j) for j in range(d)] for i in range(d)]
    while e:
        if e & 1: R = matmul(R, A, m)
        A = matmul(A, A, m); e >>= 1
    return R


def factorint(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f


t0 = time.time()
say('1  THE LOOPS, from the constants (counted; c₁ in two-level steps, so 2·c₁ levels)')
rows = []
for L in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61):
    K, poly, name = poly_of(L); d = len(poly) - 1; m = L ** K
    C = [[(-poly[d - 1 - j]) % m if i == 0 else int(j == i - 1) for j in range(d)] for i in range(d)]   # companion matrix
    fac = nmod_poly([c % L for c in poly], L).factor()[1]; degs = sorted(f.degree() for f, e in fac)
    E = 1
    for f, e in fac: E = math.lcm(E, L ** f.degree() - 1)
    if any(e > 1 for f, e in fac): E *= L
    I = [[int(i == j) for j in range(d)] for i in range(d)]
    c1 = E
    for p in factorint(E):
        while c1 % p == 0 and matpow(C, c1 // p, L) == I: c1 //= p
    D = matpow(C, c1, m); diff = [(D[i][j] - I[i][j]) % m for i in range(d) for j in range(d)]
    v = min((next((e for e in range(K) if x % L ** (e + 1)), K) for x in diff), default=K)
    lane = [c % L for c in poly]; most = max(L ** f - 1 for f in degs)
    rows.append((L, c1, v))
    say('   ℓ = %2d: %d constant%s (pieces of degree %s on lane ℓ; polynomial on lane ℓ %s; %s, %d digits): loop on lane ℓ c₁ = %d '
        '(largest possible %d%s); return depth v = %s%s' % (
            L, d, '' if d == 1 else 's', degs, lane, name, K, c1, most, ', a primitive root' if d == 1 and c1 == L - 1 else '',
            v if v < K else '≥ %d' % K, '' if v == 1 else '  ← the tower comes back to %d digits' % v))
say('   predicted share of primitive roots for one constant: 13: %.2f, 17: %.2f, 19: %.2f, 23: %.2f, 29: %.2f, 31: %.2f' % tuple(
    sum(1 for a in range(1, L) if math.gcd(a, L - 1) == 1) / (L - 1) for L in (13, 17, 19, 23, 29, 31)))

say(); say('2  THE SHORTEST LOOPS ON EXACT PARTITION NUMBERS (FLINT; counted)')
LAW = {5: lambda b: b, 7: lambda b: b // 2 + 1, 11: lambda b: b, 31: lambda b: 0}
for L, c1, alpha in ((5, 2, 4), (31, 2, 30), (11, 5, 9), (7, 6, 5)):
    t = time.time(); top = 2 + 2 * c1
    H = {}
    for b in range(1, top + 1):
        db = pow(24, -1, L ** b)
        H[b] = [int(flint.fmpz(L ** b * n + db).partitions_p()) // L ** LAW[L](b) % L for n in range(10)]
    for b in (1, 2):
        loop = H[b + 2 * c1] == H[b]; step = all((alpha * x - y) % L == 0 for x, y in zip(H[b], H[b + 2]))
        say('   ℓ = %2d, b = %d: H′_%d on lane ℓ %s;  H′_%d %s  →  back after %d levels: %s;  H′_%d ≡ %d·H′_%d: %s  (%.0f s)' % (
            L, b, b, H[b], b + 2 * c1, H[b + 2 * c1], 2 * c1, 'yes' if loop else 'NO', b + 2, alpha, b, 'yes' if step else 'NO', time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
