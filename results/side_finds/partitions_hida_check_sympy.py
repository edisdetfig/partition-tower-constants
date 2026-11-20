# A dig, 11 Nov 2026: THE 13 CHECK AGAIN, WITH A SECOND PARTITION-NUMBER PROGRAM. PARI's numbpart was found wrong once
# (p(56628877314) on lane 43: PARI 23, sympy 36; partitions_twins_43_sympy.py). partitions_hida_check.txt used PARI's values up
# to p(132114388022). Here the same ratios H_(b+2)(0)/H_b(0), b = 4 to 8, come from sympy's partition() (the program of
# partitions_constants.py) and are compared with the forms' digits of partitions_hida_digits.txt. Expected: the same
# agreement as with PARI (4 to 8 digits). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_hida_check_sympy.py [--keep]
import os, sys, time, re
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
from sympy.functions.combinatorial.numbers import partition
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_hida_check_sympy.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); L = 13; K = 11
text = open(os.path.join(HERE, 'partitions_hida_digits.txt'), encoding='utf-8').read()
block = text[text.index('ℓ = %d:' % L):]
forms = [int(d) for d in re.search(r'u\(−2\), digits lowest first: \[([^\]]*)\]', block).group(1).split(',')]
say('the constant of 13 from modular forms alone (partitions_hida_digits.txt): %s' % forms)
H = {}; M = L ** K
for b in range(4, 11):
    t = time.time(); d = pow(24, -1, L ** b); H[b] = int(partition(d)) % M
    say('   p(δ_%d) = p(%d) by sympy  (%.0f s)' % (b, d, time.time() - t))
for b in range(4, 9):
    r = H[b + 2] * pow(H[b], -1, M) % M; dg = [(r // L ** i) % L for i in range(K)]
    agree = next((i for i in range(K) if dg[i] != forms[i]), K)
    say('   levels %d and %d: the ratio\'s digits %s;  agrees with the forms\' digits in the first %d' % (b, b + 2, dg, agree))
say(); say('seconds: %.0f' % (time.time() - t0))
