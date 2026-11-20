# A dig, 13 Nov 2026: WHAT FADES — THE SPEED'S OWN NOTE (README "What fades"; predicted there first). For ℓ = 13, 17, 19,
# 23, 29, 31 (one-form classes, so the whole tower is one note plus what fades): E_b(n) = H(ℓ²M) − α·H(M), M = ℓ^b·(24n + k),
# modulo ℓ^D (D = 30; α from partitions_flint_family_ℓ.txt, 13 from the N30 file). Printed for each level b: the
# valuations v(E_b(n)) over n = 0 … 7, and the ratio E_(b+2)(n)/E_b(n): its valuation and the lowest digits of its unit
# part, and how many of those digits all n share. FLINT's exact partition numbers, at idle priority. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/tower_fading.py [--keep]
import os, sys, re, time
HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
flint.ctx.threads = 6
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'tower_fading.txt'), 'x', encoding='utf-8') if KEEP else None
NMAX = 3 * 10 ** 15; D = 30


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


PC = {}


def c(M):
    N = (M + 1) // 24
    if N not in PC: PC[N] = int(flint.fmpz(N).partitions_p())
    return PC[N]


def val(x, L, cap):
    if x == 0: return cap
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


t0 = time.time()
for L in (13, 17, 19, 23, 29, 31):
    fn = 'partitions_flint_family_13_N30.txt' if L == 13 else 'partitions_flint_family_%d.txt' % L
    txt = open(os.path.join(HERE, fn), encoding='utf-8').read()
    al = [int(x) for x in re.search(r'a constant on lane ℓ, digits lowest first: \[([^\]]*)\]', txt).group(1).split(',')]
    Dl = min(D, len(al)); m = L ** Dl; aInt = sum(x * L ** i for i, x in enumerate(al[:Dl]))
    say(); say('ℓ = %d: α to %d digits' % (L, Dl))
    E = {}
    for b in range(1, 40, 2):
        kk = (-pow(L ** b, -1, 24)) % 24; row = {}
        for n in range(8):
            M = L ** b * (24 * n + kk)
            if (L * L * M + 1) // 24 > NMAX: break
            row[n] = (c(L * L * M) - aInt * c(M)) % m
        if not row: break
        E[b] = row
    for b in sorted(E):
        vs = [val(E[b][n], L, Dl) for n in sorted(E[b])]
        line = '   level %2d: v(E) over n = %s' % (b, vs)
        if b + 2 in E:
            common = [n for n in E[b + 2] if n in E[b]]
            rats = []
            for n in common:
                e0, e2 = E[b][n], E[b + 2][n]; v0, v2 = val(e0, L, Dl), val(e2, L, Dl)
                prec = Dl - v2                                   # digits of the unit parts that are known
                u0, u2 = e0 // L ** v0, e2 // L ** v2
                mp = L ** max(1, min(Dl - v0, Dl - v2))
                ratio_unit = u2 * pow(u0 % mp, -1, mp) % mp
                rats.append((v2 - v0, ratio_unit, mp))
            vr = sorted({r[0] for r in rats})
            # digits shared by the unit parts of all ratios
            k = 0; mpmin = min(r[2] for r in rats)
            while L ** (k + 1) <= mpmin and len({r[1] % L ** (k + 1) for r in rats}) == 1: k += 1
            dg = [(rats[0][1] // L ** i) % L for i in range(min(6, k + 1))]
            line += ';  E(b + 2)/E(b): valuation %s, unit part digits %s… shared by all %d values of n: %d digit%s' % (
                vr, dg, len(rats), k, '' if k == 1 else 's')
        say(line)
say(); say('the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
