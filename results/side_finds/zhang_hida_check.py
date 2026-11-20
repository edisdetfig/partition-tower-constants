# A dig, 12 Nov 2026: AT 19, IS ZHANG'S LEFTOVER THE HIDA FAMILY? (README "Zhang's conjectures…", predicted there first.)
# From the kept towers (zhang_check_j02_N8000000_….json): the law L₄·(X − p³) = (X − u²)(X − p)(X − ū²)(X − p³u⁻⁴)(X − p³)
# at its one complete level, 5, for n = 1, 2, 3 (19⁵·n < 8·10⁶). Predicted: the leading digit of the residual runs like
# n³·a_n(ΔE₄) mod 19 (θ³ of the weight −2 member of the Hida family through ΔE₄), not like σ₃(n). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/zhang_hida_check.py [--keep]
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'zhang_hida_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def val(x, p, K):
    x %= p ** K
    if x == 0: return K
    v = 0
    while x % p == 0: x //= p; v += 1
    return v


def pmul(A, B, M):
    o = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        for j, b in enumerate(B): o[i + j] = (o[i + j] + a * b) % M
    return o


rec = json.load(open(os.path.join(HERE, 'zhang_check_j02_N8000000_5_7_11_13_17_19_23_29_31_37_41_43_307.json'), encoding='utf-8'))['primes']['19']
p, K, a = 19, rec['K'], rec['a_p']; M = p ** K
u = a % p
for _ in range(100): u = (u - (u * u - a * u + p) * pow((2 * u - a) % M, -1, M)) % M
roots = [u * u % M, p, p * p * pow(u, -2, M) % M, p ** 3 * pow(u, -4, M) % M, p ** 3]
Q = [1]
for r_ in roots: Q = pmul(Q, [(-r_) % M, 1], M)
# a_n of ΔE₄ (weight 16) and σ₃(n), exact, for the comparison
N = 8
P = [1] + [0] * (N - 1)
for m in range(1, N):
    for _ in range(24): P = [P[k] - (P[k - m] if k >= m else 0) for k in range(N)]
D = [0] + P[:N - 1]; E4 = [1] + [240 * sum(d ** 3 for d in range(1, k + 1) if k % d == 0) for k in range(1, N)]
f16 = [sum(D[i] * E4[k - i] for i in range(k + 1)) for k in range(N)]
s3 = lambda n: sum(d ** 3 for d in range(1, n + 1) if n % d == 0)
say('p = 19, lane 19^%d, a_p(C) = %d; the law L₄·(X − p³) at level 5, n = 1, 2, 3 (the only n with 19⁵·n < 8·10⁶)' % (K, a))
say('   ΔE₄ = %s…;  predicted ratios n³·a_n(ΔE₄)/a_1 mod 19: n = 2: %d, n = 3: %d;  the control σ₃(n): %d, %d' % (
    f16[1:5], 8 * f16[2] % p, 27 * f16[3] % p, s3(2) % p, s3(3) % p))
for r in ('1', '2', '3'):
    T = {int(m): v for m, v in rec['towers'][r].items()}
    R = [sum(Q[s] * T[s][n - 1] for s in range(6)) % M for n in (1, 2, 3)]
    vs = [val(x, p, K) for x in R]; d = min(vs)
    lead = [(x // p ** d) % p if v == d else None for x, v in zip(R, vs)]
    if lead[0] is None: say('   r = %s: depths %s; the n = 1 residual is deeper than the others, so no ratio' % (r, vs)); continue
    inv = pow(lead[0], -1, p)
    ratios = ['%d' % (l * inv % p) if l is not None else 'deeper' for l in lead[1:]]
    say('   r = %s: depths of R(1), R(2), R(3): %s;  leading digits %s;  R(2)/R(1), R(3)/R(1) mod 19: %s' % (r, vs, lead, ratios))
