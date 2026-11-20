# A dig, 12 Nov 2026: 7'S AND 11'S NEW DIGITS AGAINST DEEPER PARTITION NUMBERS. The FLINT slope runs gave, from modular forms
# only, 20 digits for 7 and 12 for 11 (partitions_flint_slopes_7.txt, _11.txt). The partition numbers had settled 14 and 4
# (partitions_constants.txt, sympy, levels up to 13 and 10). For 11 the forms' 5th and 6th digits (1, 7) differ from what
# PARI's three-weight series printed (5, 8), which was past that series' reach.
# PREDICTED before the run: with FLINT's exact values one level pair deeper at a time, r_b = H′_(b+2)(n)/H′_b(n) (H′ = p(ℓᵇn + δ_b)
# divided by the law's power: 7^(⌊b/2⌋+1), 11ᵇ) shares ever more digits with the forms' constant: for 11 at least 6 by
# the deepest pair (so the 5th digit is 1, not 5), for 7 at least 17. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_lawprimes_check.py [--keep]
import os, sys, re, time
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
import flint
flint.ctx.threads = 8
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_flint_lawprimes_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


LAW = {7: lambda b: b // 2 + 1, 11: lambda b: b}
TOP = {7: 18, 11: 15}
t0 = time.time()
for L in (11, 7):
    text = open(os.path.join(HERE, 'partitions_flint_slopes_%d.txt' % L), encoding='utf-8').read()
    forms = [int(v) for v in re.search(r'digits lowest first[^:]*: \[([^\]]*)\]', text).group(1).split(',')]
    K = len(forms) + 2; M = L ** K
    say(); say('ℓ = %d: the forms\' constant (partitions_flint_slopes_%d.txt), %d digits: %s' % (L, L, len(forms), forms))
    H = {}
    for b in range(1, TOP[L] + 1):
        t = time.time(); db = pow(24, -1, L ** b); row = []
        for n in range(3):
            v = int(flint.fmpz(L ** b * n + db).partitions_p())
            if v % L ** LAW[L](b): say('   b = %d, n = %d: NOT divisible by the law\'s power' % (b, n))
            row.append((v // L ** LAW[L](b)) % M)
        H[b] = row
        if b >= 3:
            a = b - 2; n = next((n for n in range(3) if H[a][n] % L), None)
            if n is None: say('   b = %d: no unit among n = 0..2 at level %d' % (b, a)); continue
            r = H[b][n] * pow(H[a][n], -1, M) % M; dg = [(r // L ** j) % L for j in range(K)]
            same = next((i for i in range(len(forms)) if dg[i] != forms[i]), len(forms))
            say('   levels (%d, %d), n = %d: r shares %2d digits with the forms\' constant; r = %s…  (%.0f s)' % (a, b, n, same, dg[:same + 2], time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
