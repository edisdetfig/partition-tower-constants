# A dig, 12 Nov 2026: 47'S AND 61'S CONSTANTS, CHECKED ON PARTITION NUMBERS. partitions_flint_family_47.txt and _61.txt give,
# from modular forms only (the FLINT engine forms6_flint.py, no PARI), the polynomial whose roots are each prime's
# constants (47: two, conjugate in the field with 47² elements; 61: three, one on lane 61 and a conjugate pair). If they are
# the tower's constants, the law Σ p_i·H_(b+2i) ≡ 0 holds ever deeper: PREDICTED before the run (as for 37, 41, 43):
# 47: lane 47 at levels (1,3,5) and (2,4,6), lane 47² at (3,5,7) and (4,6,8), lane 47³ at (5,7,9); 61: lane 61 at
# (1,3,5,7) and (2,4,6,8), lane 61² at (3,5,7,9). H_b = p(δ_b), δ_b = 1/24 on lane ℓᵇ (n = 0, and n = 1 for 47);
# FLINT's partition numbers on 8 threads. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_hida_more_check.py [--keep]
import os, sys, time, re
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
import flint
flint.ctx.threads = 8
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_hida_more_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); PC = {}


def P(n):
    if n not in PC: PC[n] = int(flint.fmpz(n).partitions_p())
    return PC[n]


for L, top in ((47, 9), (61, 9)):
    text = open(os.path.join(HERE, 'partitions_flint_family_%d.txt' % L), encoding='utf-8').read()
    mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\]', text); K = int(mt.group(1)); poly = [int(v) for v in mt.group(2).split(',')]
    d = len(poly) - 1; MK = L ** K; t = time.time()
    val = lambda x: next((e for e in range(K) if x % L ** (e + 1)), K)
    say(); say('ℓ = %d: the constants\' polynomial on lane ℓ^%d (partitions_flint_family_%d.txt): %s (degree %d)' % (L, K, L, poly, d))
    for n in ((0, 1) if L == 47 else (0,)):
        H = {b: P(L ** b * n + pow(24, -1, L ** b)) % MK for b in range(1, top + 1)}
        rows = []
        for b in range(1, top - 2 * d + 1):
            r = sum(poly[i] * H[b + 2 * i] for i in range(d + 1)) % MK
            rows.append('(%s): ℓ^%d' % (','.join(str(b + 2 * i) for i in range(d + 1)), val(r)))
        say('   n = %d: the law is divisible by  %s  (%.0f s)' % (n, ';  '.join(rows), time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
