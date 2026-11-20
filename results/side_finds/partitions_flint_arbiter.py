# A dig, 12 Nov 2026: FLINT AS THE ARBITER. PARI's numbpart and sympy's partition() disagreed on p(56628877314) on lane 43
# (PARI 23 and 62975 on lane 43³, computed twice; sympy 36 and 63332; partitions_twins_43_sympy.txt, README "a thirteenth
# time"). FLINT's partition function (python-flint 0.9.0, partitions_p, Arb's ball arithmetic: FLINT 3 documents the result
# as exact) settles it. Then every PARI value that a kept result rests on is recomputed with FLINT:
#   1  the disputed value;
#   2  the law members and controls of partitions_split.txt (37 at Q = 6599, 7541, 9547; 43 at 4049, 8941), read from that
#      file with PARI's residues, compared one by one;
#   3  the twins' law along ℓ (partitions_twins.txt for 37 and 41, partitions_twins_43.txt for 43; constants read from
#      them), with p(δ_b), b = 1 to 7, from FLINT;
#   4  13's deep check (partitions_hida_check.txt), the ratio digits of the levels two apart, against the forms' digits
#      read from partitions_hida_digits.txt;
#   5  the law at 2999 for 37 directly: p(540,584,055,355), past sympy's limit, PREDICTED a multiple of 37
#      (partitions_split_2999.txt checked its cause).
# COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_arbiter.py [--keep]
import os, sys, time, re
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
import flint
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_flint_arbiter.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); PC = {}
read = lambda name: open(os.path.join(HERE, name), encoding='utf-8').read()


def P(n):
    if n not in PC: PC[n] = int(flint.fmpz(n).partitions_p())
    return PC[n]


say('python-flint %s' % flint.__version__)

# ---- 1. the disputed value ----
t = time.time(); v = P(56628877314)
say('1  p(56628877314) by FLINT: %d on lane 43, %d on lane 43³  (%.0f s)' % (v % 43, v % 43 ** 3, time.time() - t))
say('   PARI gave 23 and 62975 (twice), sympy 36 and 63332  →  %s' % (
    'sympy is right, PARI is wrong' if (v % 43, v % 43 ** 3) == (36, 63332) else
    'PARI is right, sympy is wrong' if (v % 43, v % 43 ** 3) == (23, 62975) else 'neither'))

# ---- 2. the law members of partitions_split.txt ----
L = None; rows = []; agree = tot = 0
for line in read('partitions_split.txt').splitlines():
    m = re.match(r'ℓ = (\d+):', line)
    if m: L = int(m.group(1))
    if 'members p(' in line:
        for n, r in re.findall(r'p\((\d+)\) ≡ (\d+)', line):
            n, r = int(n), int(r); f = P(n) % L; tot += 1; agree += f == r
            if f != r: rows.append('ℓ = %d, p(%d): PARI %d, FLINT %d' % (L, n, r, f))
say('2  the law members and controls of partitions_split.txt: FLINT agrees with PARI at %d of %d%s' % (
    agree, tot, '' if not rows else ';  differ: ' + ';  '.join(rows)))

# ---- 3. the twins' law along ℓ ----
for L, name in ((37, 'partitions_twins.txt'), (41, 'partitions_twins.txt'), (43, 'partitions_twins_43.txt')):
    text = read(name); block = text[text.index('ℓ = %d, twin 1' % L):]
    digs = [[int(d) for d in s.split(',')] for s in re.findall(r'digits lowest first: \[([^\]]*)\]', block)[:2]]
    K = min(len(d) for d in digs); MK = L ** K
    a1, a2 = (sum(d[i] * L ** i for i in range(K)) for d in digs); T, N = (a1 + a2) % MK, a1 * a2 % MK
    H = {b: P(pow(24, -1, L ** b)) % MK for b in range(1, 8)}
    val = lambda x: next((e for e in range(K) if x % L ** (e + 1)), K)
    say('3  ℓ = %d (constants from %s): H_(b+4) − T·H_(b+2) + N·H_b divisible by ℓ^ %s for levels (1,3,5), (2,4,6), (3,5,7)' % (
        L, name, [val((H[b + 4] - T * H[b + 2] + N * H[b]) % MK) for b in (1, 2, 3)]))

# ---- 4. 13's deep check ----
text = read('partitions_hida_digits.txt'); block = text[text.index('ℓ = 13:'):]
forms = [int(d) for d in re.search(r'u\(−2\), digits lowest first: \[([^\]]*)\]', block).group(1).split(',')]
L, K = 13, 11; M = L ** K; H = {b: P(pow(24, -1, L ** b)) % M for b in range(4, 11)}
out = []
for b in range(4, 9):
    r = H[b + 2] * pow(H[b], -1, M) % M; dg = [(r // L ** i) % L for i in range(K)]
    out.append('(%d, %d): %d' % (b, b + 2, next((i for i in range(K) if dg[i] != forms[i]), K)))
say('4  13: digits shared with the forms by the ratio of levels two apart (FLINT):  %s' % ';  '.join(out))

# ---- 5. the law at 2999 ----
Q = 2999; m = next(m for m in range(1, 200) if (37 * Q ** 3 * m + 1) % 24 == 0 and m % Q)
N = (37 * Q ** 3 * m + 1) // 24; t = time.time(); v = P(N) % 37
say('5  p(%d) (24N − 1 = 37·2999³·%d) by FLINT: %d on lane 37  →  %s  (%.0f s)' % (
    N, m, v, 'a multiple of 37, as derived' if v == 0 else 'NOT a multiple', time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
