# A try, 14 Nov 2026. Split 607's quiet
# braids at the tower's weight into their strands, on the first digit, with numbers from each braid's own field.
# Input: the polynomials fading_forms4.py printed for weight -2 (fading_forms4_607_s1.txt), first digits, as printed and
# with the run's three-weight curve undone ((-1)^j/c_j, as in fading_forms4_607_correct.py).
# For each factor f of degree k on lane 607, the field is F = F_607[t]/(f(t)). Its k roots are t, t^607, t^(607^2), ...:
# the strands. Each is printed as its coefficients c0, c1, ... of c0 + c1*t + ... (mod 607). Checked for every strand: it
# is a root of f; all k are distinct; the 607th power steps from each to the next and the k-th returns to the first.
# Also printed: the strands' product (the braid's size) and sum, which must be ordinary numbers mod 607.
# COUNTED (the field arithmetic); the input polynomials are DERIVED (README "607's quiet braids").
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/fading_forms4_607_split.py
import os, re, ast
from flint import nmod_poly
HERE = os.path.dirname(os.path.abspath(__file__)); L = 607
txt = open(os.path.join(HERE, 'fading_forms4_607_s1.txt'), encoding='utf-8').read()
fh = open(os.path.join(HERE, 'fading_forms4_607_split.txt'), 'x', encoding='utf-8')


def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); fh.write(s + '\n'); fh.flush()


P = lambda co: nmod_poly([int(v) % L for v in co], L)
c = lambda j: (6 - 8 * 2 ** j + 3 * 3 ** j) % L


def powmod(a, e, f):
    r, b = P([1]), a % f
    while e:
        if e & 1: r = (r * b) % f
        b = (b * b) % f; e >>= 1
    return r


def at(f, s):                     # f(s), computed in F_607[t]/(f)
    acc = P([0])
    for g in reversed([int(v) for v in f.coeffs()]): acc = (acc * s + P([g])) % f
    return acc


def vec(a, k):
    co = [int(v) for v in a.coeffs()]
    return co + [0] * (k - len(co))


def scalar(a, k):
    v = vec(a, k)
    return v[0] if all(x == 0 for x in v[1:]) else v


def split(name, co):
    G = P(co); d = len(co) - 1
    facs = [(f, e) for f, e in G.factor()[1]]
    say('%s (lowest first): %s' % (name, co))
    say('   factors on lane 607, of degrees %s' % sorted(f.degree() for f, e in facs for _ in range(e)))
    for f, e in sorted(facs, key=lambda fe: fe[0].degree()):
        k = f.degree(); fl = [int(v) for v in f.coeffs()]
        if k == 1:
            say('   - a lone root on lane 607 itself, an ordinary number: %d' % ((-fl[0]) % L)); continue
        s = [P([0, 1]) % f]
        for i in range(k): s.append(powmod(s[-1], L, f))
        roots = all(at(f, x).degree() < 0 for x in s[:k])
        distinct = len({tuple(vec(x, k)) for x in s[:k]}) == k
        closes = s[k] == s[0]
        prod, sm = P([1]), P([0])
        for x in s[:k]: prod = (prod * x) % f; sm = (sm + x) % f
        say('   - a braid of %d strands, in the field with 607^%d elements; t is a root of %s' % (k, k, fl))
        say('     every strand a root: %s; all %d distinct: %s; the 607th power steps through them and returns: %s' % (roots, k, distinct, closes))
        say('     the strands multiply to %s and add up to %s' % (scalar(prod, k), scalar(sm, k)))
        for i, x in enumerate(s[:k]):
            say('     strand %2d = t^(607^%d) = %s' % (i + 1, i, vec(x, k)))


say('607: the quiet braids at weight -2 (first digit), split in their own fields; input fading_forms4_607_s1.txt')
for m in re.finditer(r'a group of (\d+) roots.*?digits lowest first: (\[\[.*?\]\])', txt, re.S):
    d = int(m.group(1)); co = [v[0] for v in ast.literal_eval(m.group(2))]
    corr = [(co[i] * (-1) ** (d - i) * pow(c(d - i), -1, L)) % L for i in range(d + 1)]
    say(); say('=================== the braid of %d ===================' % d)
    split('as the run printed it', co)
    split('with the three-weight curve undone', corr)
fh.close()
