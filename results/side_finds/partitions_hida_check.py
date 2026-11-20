# A dig, 11 Nov 2026: THE NEW DIGITS OF THE 13 CONSTANT, CHECKED ON PARTITION NUMBERS. partitions_hida_digits.txt computes the
# constant of 13 from modular forms alone, to 11 digits; the partition numbers of partitions_constants.txt had settled 6.
# Here the ratio of the levels two apart, H_(b+2)(0)/H_b(0) with H_b(0) = p(δ_b), δ_b = 1/24 on lane 13ᵇ, is computed from
# exact partition numbers (PARI's numbpart) up to b = 10 (N = 132,114,388,022), and its digits are compared with the
# forms' digits, read from partitions_hida_digits.txt. PREDICTED before the run: each deeper pair agrees with the forms' digits
# in more places, and the pair (8, 10) reaches digits 7 and 8, which no partition number had reached. COUNTED.
#   usage (WSL): PYTHONPATH=<folder with cypari2> python3 -B -s results/side_finds/partitions_hida_check.py [--keep]
import os, sys, time, re
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(4 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_hida_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); L = 13; K = 11
text = open(os.path.join(HERE, 'partitions_hida_digits.txt'), encoding='utf-8').read()
block = text[text.index('ℓ = %d:' % L):]
forms = [int(d) for d in re.search(r'u\(−2\), digits lowest first: \[([^\]]*)\]', block).group(1).split(',')]
say('the constant of %d from modular forms alone (partitions_hida_digits.txt), digits lowest first: %s' % (L, forms))
H = {}
for b in range(4, 11):
    t = time.time(); d = pow(24, -1, L ** b); g = pari.numbpart(d); H[b] = int(g)
    say('   p(δ_%d) = p(%d): %d digits, exact  (%.0f s)' % (b, d, int(pari.logint(g, 10)) + 1, time.time() - t))
M = L ** K
for b in range(4, 9):
    r = H[b + 2] * pow(H[b], -1, M) % M; dg = [(r // L ** i) % L for i in range(K)]
    agree = next((i for i in range(K) if dg[i] != forms[i]), K)
    say('   levels %d and %d: the ratio\'s digits %s;  agrees with the forms\' digits in the first %d' % (b, b + 2, dg, agree))
say(); say('seconds: %.0f' % (time.time() - t0))
