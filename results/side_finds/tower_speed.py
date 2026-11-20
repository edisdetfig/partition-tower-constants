# A dig, 13 Nov 2026: HOW FAST DOES EACH TOWER REACH ITS CONSTANT? (README "The speed of the towers"; predicted there
# first). For ℓ = 13, 17, 19, 23, 29, 31 the class is one form, so the whole tower is a single note: the multiplier
# H(ℓ²M)/H(M), M = ℓ^b·(24n + k), k ≡ −ℓ^(−b) mod 24, against the constant α (partitions_flint_family_ℓ.txt, 17 to 30
# digits). Counted: the digits the multiplier shares with α at levels b → b + 2, for n = 0 … 7 (smallest over n), as deep
# as p(N) with N ≤ 3·10¹⁵ allows. FLINT's exact partition numbers, at idle priority. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/tower_speed.py [--keep]
import os, sys, re, time
HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
flint.ctx.threads = 4
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'tower_speed.txt'), 'x', encoding='utf-8') if KEEP else None
NMAX = 3 * 10 ** 15; D = 14


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


PC = {}


def c(M):
    N = (M + 1) // 24
    if N not in PC: PC[N] = int(flint.fmpz(N).partitions_p())
    return PC[N]


t0 = time.time()
for L in (13, 17, 19, 23, 29, 31):
    txt = open(os.path.join(HERE, 'partitions_flint_family_%d.txt' % L), encoding='utf-8').read()
    al = [int(x) for x in re.search(r'a constant on lane ℓ, digits lowest first: \[([^\]]*)\]', txt).group(1).split(',')]
    Dl = min(D, len(al)); m = L ** Dl; aInt = sum(x * L ** i for i, x in enumerate(al[:Dl]))
    say(); say('ℓ = %d: the constant α = %s … (%d digits used)' % (L, ', '.join(map(str, al[:8])), Dl))
    for b in range(1, 40, 2):
        kk = (-pow(L ** b, -1, 24)) % 24; t = time.time(); got = []; skipped = 0
        for n in range(8):
            M = L ** b * (24 * n + kk)
            if (L * L * M + 1) // 24 > NMAX: break
            h0, h2 = c(M) % m, c(L * L * M) % m
            if h0 % L == 0: skipped += 1; continue
            r = h2 * pow(h0, -1, m) % m; diff = (r - aInt) % m; v = 0
            while v < Dl and diff % L == 0: diff //= L; v += 1
            got.append(v)
        if not got and not skipped: say('   levels %d → %d: past N = %.0e, not computed' % (b, b + 2, NMAX)); break
        say('   levels %2d → %2d: digits shared with α, n = 0 … %d: %s  (smallest %s)%s  (%.0f s)' % (
            b, b + 2, len(got) + skipped - 1, got, min(got) if got else '—', ('; %d skipped (H ≡ 0 on lane ℓ)' % skipped) if skipped else '',
            time.time() - t))
say(); say('the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
