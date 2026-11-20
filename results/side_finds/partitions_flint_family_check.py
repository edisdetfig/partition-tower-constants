# A dig, 12 Nov 2026: 53'S AND 59'S CONSTANTS, CHECKED ON PARTITION NUMBERS. The FLINT runs partitions_flint_family_53.txt
# and _59.txt (modular forms only, 20 weights, modulo ℓ^20) gave each prime's two constants: 53's both on lane 53, 59's a
# pair A ± B·√2 (no √2 on lane 59). If they are the tower's constants, the law Σ c_i·H_(b+2i) ≡ 0 (c: their polynomial)
# holds ever deeper, as it did for 37, 41, 43, 47 (partitions_flint_deeper.txt, partitions_hida_more_check.txt).
# PREDICTED before the run: lane ℓ at the levels (1,3,5) and (2,4,6), lane ℓ² at (3,5,7) and (4,6,8), lane ℓ³ at (5,7,9),
# for n = 0 and n = 1 (H_b = p(ℓᵇ·n + δ_b), δ_b = 1/24 on lane ℓᵇ). FLINT's exact partition numbers on 8 threads. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_family_check.py ℓ … [--keep]
import os, sys, time, re, glob
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
import flint
flint.ctx.threads = 8
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
PS = [int(a) for a in sys.argv[1:] if a.isdigit()]
fh = open(os.path.join(HERE, 'partitions_flint_family_check_%s.txt' % '_'.join(map(str, PS))), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); PC = {}


def P(n):
    if n not in PC: PC[n] = int(flint.fmpz(n).partitions_p())
    return PC[n]


TOP = 9
for L in PS:
    best = None
    for f in glob.glob(os.path.join(HERE, 'partitions_flint_family_%d*.txt' % L)):
        mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\]', open(f, encoding='utf-8').read())
        if mt and (best is None or int(mt.group(1)) > best[0]): best = (int(mt.group(1)), [int(v) for v in mt.group(2).split(',')], os.path.basename(f))
    K, poly, name = best; d = len(poly) - 1; MK = L ** K
    val = lambda x: next((e for e in range(K) if x % L ** (e + 1)), K)
    say(); say('ℓ = %d: the constants\' polynomial on lane ℓ^%d (%s), degree %d' % (L, K, name, d))
    for n in (0, 1):
        t = time.time()
        H = {b: P(L ** b * n + pow(24, -1, L ** b)) % MK for b in range(1, TOP + 1)}
        rows = []
        for b in range(1, TOP - 2 * d + 1):
            r = sum(poly[i] * H[b + 2 * i] for i in range(d + 1)) % MK
            rows.append('(%s): ℓ^%d' % (','.join(str(b + 2 * i) for i in range(d + 1)), val(r)))
        say('   n = %d: the law is divisible by  %s  (%.0f s)' % (n, ';  '.join(rows), time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
