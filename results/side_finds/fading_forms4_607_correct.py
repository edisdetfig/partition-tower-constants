# A check, 14 Nov 2026 (README "607's quiet braids", "PREDICTED before the 2,422 result"): undo the last step of
# fading_forms4.py on the first digit, for the braids it carried to weight -2.
# The tool's last step is Lagrange through x = 1, 2, 3 read at x = -1. On x^j it gives c_j = 6 - 8*2^j + 3*3^j instead of
# (-1)^j. At 29 and 31 the quiet note's first digit is (constant)*x at all 7 nearby weights (counted), so for a braid the
# coefficient of T^(d-j) is (constant)*x^j on the first digit, and multiplying the printed first digit by (-1)^j/c_j
# undoes the step. Reads the printed polynomials from fading_forms4_607_s1.txt; writes fading_forms4_607_correct.txt.
# COUNTED (the factorizations), DERIVED (the correction).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/fading_forms4_607_correct.py
import os, re, ast
from flint import nmod_poly
HERE = os.path.dirname(os.path.abspath(__file__)); L = 607
txt = open(os.path.join(HERE, 'fading_forms4_607_s1.txt'), encoding='utf-8').read()
fh = open(os.path.join(HERE, 'fading_forms4_607_correct.txt'), 'x', encoding='utf-8')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); fh.write(s + '\n')
c = lambda j: (6 - 8 * 2 ** j + 3 * 3 ** j) % L
fac = lambda p: sorted(f.degree() for f, e in p.factor()[1] for _ in range(e))
say('607: the braids carried to weight -2 by fading_forms4.py (plain Lagrange through x = 1, 2, 3), first digits')
for m in re.finditer(r'a group of (\d+) roots.*?digits lowest first: (\[\[.*?\]\])\n\s+its factors on lane ℓ: degrees (\[[^\]]*\])', txt, re.S):
    d = int(m.group(1)); co = [v[0] for v in ast.literal_eval(m.group(2))]; printed_fac = m.group(3)
    assert len(co) == d + 1 and co[-1] == 1
    corr = [(co[i] * (-1) ** (d - i) * pow(c(d - i), -1, L)) % L for i in range(d + 1)]
    prod = ((-1) ** d * corr[0]) % L; tr = (-corr[d - 1]) % L
    say('braid of %d:' % d)
    say('   printed, first digits lowest first:', co)
    say('   printed, factor degrees on lane 607:', fac(nmod_poly(co, L)), '(the run printed %s)' % printed_fac)
    say('   corrected by (-1)^j/c_j, lowest first:', corr)
    say('   corrected, factor degrees on lane 607:', fac(nmod_poly(corr, L)))
    say('   corrected: product of the roots %d, sum of the roots %d (derived from its coefficients)' % (prod, tr))
fh.close()
