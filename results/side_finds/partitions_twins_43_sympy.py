# A dig, 11 Nov 2026: THE 43 TWINS CHECK, WITH A SECOND PARTITION-NUMBER PROGRAM. partitions_twins_43.txt found the two-number
# law failing on lane 43 at levels 3, 5, 7, which the theory forbids (every odd level obeys the lane-43 law). Rechecked:
# PARI's numbpart(56628877314) on lane 43 is 23 (a scratch rerun gave the same), while sympy's partition() gives 36, the
# value the law needs. The smaller values agree in both programs. Here every level H_b(0) = p(δ_b), b = 1 to 7, is taken
# from sympy (Hardy–Ramanujan–Rademacher through mpmath, the program of partitions_constants.py), and the law is checked
# with the twins' constants read from partitions_twins_43.txt (forms only, PARI). PREDICTED (as for 37 and 41): lane 43 at
# levels 1, 3, 5 and 2, 4, 6; lane 43² at levels 3, 5, 7. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_twins_43_sympy.py [--keep]
import os, sys, time, re
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
from sympy.functions.combinatorial.numbers import partition
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_twins_43_sympy.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); L = 43
text = open(os.path.join(HERE, 'partitions_twins_43.txt'), encoding='utf-8').read()
digs = [[int(d) for d in s.split(',')] for s in re.findall(r'digits lowest first: \[([^\]]*)\]', text)]
K = min(len(d) for d in digs); MK = L ** K
a1, a2 = (sum(d[i] * L ** i for i in range(K)) for d in digs[:2]); T, N = (a1 + a2) % MK, a1 * a2 % MK
say('the twins\' constants (partitions_twins_43.txt): %s and %s; on lane 43: T = %d, N = %d' % (digs[0], digs[1], T % L, N % L))
val = lambda x: next((e for e in range(K) if x % L ** (e + 1)), K)
H = {}
for b in range(1, 8):
    t = time.time(); d = pow(24, -1, L ** b); v = int(partition(d)); H[b] = v % MK
    say('   p(δ_%d) = p(%d) by sympy: on lane 43: %d  (%.0f s)' % (b, d, v % L, time.time() - t))
for b in (1, 2, 3):
    r = (H[b + 4] - T * H[b + 2] + N * H[b]) % MK
    say('   levels %d, %d, %d: H_(b+4) − T·H_(b+2) + N·H_b is divisible by 43^%d' % (b, b + 2, b + 4, val(r)))
say(); say('seconds: %.0f' % (time.time() - t0))
