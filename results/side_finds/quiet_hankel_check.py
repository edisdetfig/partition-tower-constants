# A dig, 13 Nov 2026: THE QUIET NOTES WITHOUT THE LOUD NOTE'S CONSTANT (README section of that name; predicted there
# first). For ℓ = 13 and 17: H_b(n) = p(ℓ^b·n + δ_b) for b = 1, 3, …, 11 and n = 0 … 7 (exact partition numbers from FLINT,
# kept modulo ℓ^30), then D_b(n) = H_(b+4)·H_b − H_(b+2)², and the ratio D_(b+2)(n)/D_b(n): its valuation and the digits of
# its unit part that all n share. No constant is used to form them. Compared with α·β/ℓ², α from
# partitions_flint_family_ℓ.txt (13 from the N30 file) and β/ℓ² from tower_fading.txt. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/quiet_hankel_check.py [--keep]
import os, sys, re, time
HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
flint.ctx.threads = 2
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'quiet_hankel_check.txt'), 'x', encoding='utf-8') if KEEP else None
D, NMAX = 30, 3 * 10 ** 14
BETA = {13: [10, 11, 6, 2, 10, 10], 17: [3, 9, 1, 13, 4, 3]}


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


t0 = time.time()
for L in (13, 17):
    fn = 'partitions_flint_family_13_N30.txt' if L == 13 else 'partitions_flint_family_%d.txt' % L
    al = [int(x) for x in re.search(r'a constant on lane ℓ, digits lowest first: \[([^\]]*)\]', open(os.path.join(HERE, fn), encoding='utf-8').read()).group(1).split(',')]
    m = L ** D; alpha = sum(x * L ** i for i, x in enumerate(al[:D])); bu = sum(x * L ** i for i, x in enumerate(BETA[L]))
    ab = alpha * bu % L ** len(BETA[L]); abd = [(ab // L ** i) % L for i in range(len(BETA[L]))]
    H = {}
    for b in range(1, 12, 2):
        kk = (-pow(L ** b, -1, 24)) % 24; row = {}
        for n in range(8):
            M = L ** b * (24 * n + kk)
            if (M + 1) // 24 > NMAX: break
            row[n] = c(M) % m
        H[b] = row
    say(); say('ℓ = %d: α·β/ℓ² from the known constants, digits lowest first: %s' % (L, abd))
    Dv = {}
    for b in range(1, 8, 2):
        if b + 4 not in H: break
        Dv[b] = {n: (H[b + 4][n] * H[b][n] - H[b + 2][n] ** 2) % m for n in H[b + 4] if n in H[b] and n in H[b + 2]}
    for b in sorted(Dv):
        vs = [val(Dv[b][n], L, D) for n in sorted(Dv[b])]
        line = '   levels (%d, %d, %d): v(D) over n = %s' % (b, b + 2, b + 4, vs)
        if b + 2 in Dv:
            rats = []
            for n in Dv[b + 2]:
                if n not in Dv[b]: continue
                d0, d2 = Dv[b][n], Dv[b + 2][n]; v0, v2 = val(d0, L, D), val(d2, L, D)
                mp = L ** max(1, min(D - v0, D - v2) - v0)                 # products of H's: known to ℓ^(D + v(H)) only
                u0, u2 = d0 // L ** v0, d2 // L ** v2
                rats.append((v2 - v0, u2 * pow(u0 % mp, -1, mp) % mp, mp))
            k = 0; mpmin = min(r[2] for r in rats)
            while L ** (k + 1) <= mpmin and len({r[1] % L ** (k + 1) for r in rats}) == 1: k += 1
            dg = [(rats[0][1] // L ** i) % L for i in range(min(6, k))]
            agree = next((i for i in range(min(k, len(abd))) if dg[i] != abd[i]), min(k, len(abd)))
            line += ';  D(next)/D: valuation %s, unit part digits shared by all %d values of n: %s (%d digits);  against α·β/ℓ²: %d agree' % (
                sorted({r[0] for r in rats}), len(rats), dg, k, agree)
        say(line)
say(); say('the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
