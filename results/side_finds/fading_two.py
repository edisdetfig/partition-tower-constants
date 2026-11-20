# A dig, 13 Nov 2026: WHAT FADES AT 19 AND 23 — TWO QUIET NOTES (README "What fades at 19 and 23: two quiet notes at once";
# predicted there first). As in tower_fading.py: E_b(n) = H(ℓ²M) − α·H(M), M = ℓ^b·(24n + k), with α from
# partitions_flint_family_ℓ.txt to 30 digits, so E is known modulo ℓ^30. If what fades is exactly two notes, then
# E_(b+4)(n) = s·E_(b+2)(n) − p·E_b(n) for every n, with one s and one p (s = β₁ + β₂, p = β₁·β₂).
# For each level triple (b, b + 2, b + 4) and each pair of values of n, the two equations are solved for s and p, each to
# the digits its valuations allow. Printed: the digits of s and p that every pair knowing them shares (and how many pairs
# know the last one); the rule tested on every n with those values; and the two notes β/ℓ², the roots of
# z² − (s/ℓ²)·z + p/ℓ⁴ (on lane ℓ, or a conjugate pair when the discriminant is not a square there).
# FLINT's exact partition numbers, at idle priority. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/fading_two.py [--keep]
import os, sys, re, time
HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
flint.ctx.threads = 6
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'fading_two.txt'), 'x', encoding='utf-8') if KEEP else None
NMAX = 3 * 10 ** 15; D = 30; NS = 8


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


PC = {}


def c(M):
    N = (M + 1) // 24
    if N not in PC: PC[N] = int(flint.fmpz(N).partitions_p())
    return PC[N]


def val(x, L, cap):
    x %= L ** cap
    if x == 0: return cap
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def solve(a, b, L):
    """a, b: (E_b, E_(b+2), E_(b+4)) at two values of n, each known modulo ℓ^D. Returns (s, p, prec): s and p modulo
    ℓ^prec, or None when the two equations do not determine them, or 'no fit' when s or p would not be integral."""
    vmin = min(val(x, L, D) for x in a + b); P = D + vmin; mP = L ** P   # products of two E's are known modulo ℓ^P
    det = (a[0] * b[1] - a[1] * b[0]) % mP
    ns = (a[0] * b[2] - a[2] * b[0]) % mP
    np_ = (a[1] * b[2] - a[2] * b[1]) % mP
    vd = val(det, L, P)
    if vd >= P: return None
    if val(ns, L, P) < vd or val(np_, L, P) < vd: return 'no fit'
    prec = P - vd; mq = L ** prec; ui = pow((det // L ** vd) % mq, -1, mq)
    return (ns // L ** vd) * ui % mq, (np_ // L ** vd) * ui % mq, prec


def consensus(vals, L):
    """vals: (x, prec) pairs. The lowest digits that every value knowing them shares: (x mod ℓ^K, K, pairs knowing digit K)."""
    K = 0
    while True:
        know = [x for x, pr in vals if pr > K]
        if len(know) < 2 or len({x % L ** (K + 1) for x in know}) != 1: break
        K += 1
    know = [x for x, pr in vals if pr >= K]
    return (know[0] % L ** K if know else 0), K, len([1 for x, pr in vals if pr >= K])


def digits(x, L, lo, hi):
    return [(x // L ** i) % L for i in range(lo, hi)]


t0 = time.time()
for L in (19, 23):
    txt = open(os.path.join(HERE, 'partitions_flint_family_%d.txt' % L), encoding='utf-8').read()
    al = [int(x) for x in re.search(r'a constant on lane ℓ, digits lowest first: \[([^\]]*)\]', txt).group(1).split(',')]
    m = L ** D; aInt = sum(x * L ** i for i, x in enumerate(al[:D]))
    say(); say('ℓ = %d: α to %d digits; (3|ℓ) = %+d' % (L, D, leg(3, L)))
    E = {}
    for b in range(1, 40, 2):
        kk = (-pow(L ** b, -1, 24)) % 24; row = {}
        for n in range(NS):
            M = L ** b * (24 * n + kk)
            if (L * L * M + 1) // 24 > NMAX: break
            row[n] = (c(L * L * M) - aInt * c(M)) % m
        if not row: break
        E[b] = row
        say('   level %2d: v(E) over n = %s  (%.0f s)' % (b, [val(row[n], L, D) for n in sorted(row)], time.time() - t0))
    for b in sorted(E):
        if b + 2 not in E or b + 4 not in E: continue
        ns_ = [n for n in range(NS) if n in E[b] and n in E[b + 2] and n in E[b + 4]]
        res = []; nofit = 0; undet = 0
        for i, n1 in enumerate(ns_):
            for n2 in ns_[i + 1:]:
                r = solve((E[b][n1], E[b + 2][n1], E[b + 4][n1]), (E[b][n2], E[b + 2][n2], E[b + 4][n2]), L)
                if r is None: undet += 1
                elif r == 'no fit': nofit += 1
                else: res.append(r)
        say('   levels (%d, %d, %d), n = %s: %d pairs solved, %d not determined, %d with no integral fit' % (
            b, b + 2, b + 4, ns_, len(res), undet, nofit))
        if len(res) < 2: continue
        s, Ks, cs = consensus([(r[0], r[2]) for r in res], L)
        p, Kp, cp = consensus([(r[1], r[2]) for r in res], L)
        vs_, vp_ = val(s, L, Ks), val(p, L, Kp)
        say('      s: %d digits shared by every pair that knows them (%d pairs know the last); v(s) %s %d; s/ℓ², digits lowest first: %s' % (
            Ks, cs, '=' if vs_ < Ks else '≥', vs_, digits(s, L, 2, Ks) if Ks > 2 else '—'))
        say('      p: %d digits shared (%d pairs know the last); v(p) %s %d; p/ℓ⁴, digits lowest first: %s' % (
            Kp, cp, '=' if vp_ < Kp else '≥', vp_, digits(p, L, 4, Kp) if Kp > 4 else '—'))
        tests = []
        for n in ns_:
            e0, e2, e4 = E[b][n], E[b + 2][n], E[b + 4][n]
            bound = min(D, Ks + val(e2, L, D), Kp + val(e0, L, D))   # the residual is known to ℓ^bound
            R = (e4 - s * e2 + p * e0) % L ** bound
            tests.append('n=%d: v(E_%d) = %d, rule holds to ℓ^%d of ℓ^%d' % (n, b + 4, val(e4, L, D), val(R, L, bound), bound))
        say('      the rule on every n: ' + '; '.join(tests))
        if Ks > 2 and Kp > 4 and vs_ >= 2 and vp_ == 4:
            S = s // L ** 2; Pn = p // L ** 4; K = min(Ks - 2, Kp - 4); mK = L ** K
            dl = (S * S - 4 * Pn) % L
            if dl == 0: say('      the two notes β/ℓ²: they agree on lane ℓ (the discriminant is 0 there)')
            elif leg(dl, L) == -1: say('      the two notes β/ℓ²: a conjugate pair (the discriminant is not a square on lane ℓ); sum %s, product %s' % (
                digits(S, L, 0, K), digits(Pn, L, 0, K)))
            else:
                roots = []
                for z0 in range(L):
                    if (z0 * z0 - S * z0 + Pn) % L: continue
                    z = z0
                    for _ in range(K + 2): z = (z - (z * z - S * z + Pn) * pow((2 * z - S) % mK, -1, mK)) % mK
                    roots.append(digits(z, L, 0, K))
                say('      the two notes β/ℓ² (to %d digits), digits lowest first: %s' % (K, roots))
say(); say('the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
