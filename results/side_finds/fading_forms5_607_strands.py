# The last step of the plan, 14 Nov 2026 (README "607's quiet braids"): each strand of 607's two quiet braids, refined
# from the first digit to the three digits the join gives (the first two confirmed, the third likely).
# Input: the scaled carry in fading_forms5_607_s1_fit.txt, the polynomial G of the constants (3|607)·u/607 at weight −2,
# modulo 607^3. Its strands are its d roots in the ring R = (Z/607^3)[t]/(G(t)), the unramified extension of degree d:
#   strand 1 is t itself; strand i+1 is the root congruent to (strand i)^607 on the first digit, found by Newton's
#   method in R (the roots are distinct on the first digit, so each lifts uniquely: Hensel).
# Each strand is printed as its coefficients c0, c1, …, c_(d−1) of c0 + c1·t + … in R, each as digits lowest first.
# Checked: each strand is a root of G modulo 607^3; the d strands are distinct; their product and sum are ordinary numbers
# (no t), equal to (−1)^d·G(0) and −(the coefficient of t^(d−1)). DERIVED from the join's polynomials.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/fading_forms5_607_strands.py
import os, re, ast
HERE = os.path.dirname(os.path.abspath(__file__)); L = 607; D = 3; M = L ** D
txt = open(os.path.join(HERE, 'fading_forms5_607_s1_fit.txt'), encoding='utf-8').read()
fh = open(os.path.join(HERE, 'fading_forms5_607_strands.txt'), 'x', encoding='utf-8')


def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); fh.write(s + '\n'); fh.flush()


def dig(v): return [(v // L ** i) % L for i in range(D)]


def ring(G):                                # arithmetic in (Z/M)[t]/(G), G monic, lowest first
    d = len(G) - 1
    def red(c):
        c = [x % M for x in c]
        for k in range(len(c) - 1, d - 1, -1):
            q = c[k]
            if q:
                for i in range(d + 1): c[k - d + i] = (c[k - d + i] - q * G[i]) % M
        return (c + [0] * d)[:d]
    def mul(a, b):
        c = [0] * (2 * d)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b): c[i + j] += x * y
        return red(c)
    def add(a, b): return [(x + y) % M for x, y in zip(a, b)]
    def const(v): return red([v % M])
    def pw(a, e):
        r, b = const(1), a
        while e:
            if e & 1: r = mul(r, b)
            b = mul(b, b); e >>= 1
        return r
    def at(P, s):                           # P(s) for a polynomial P (lowest first) with ordinary coefficients
        acc = const(0)
        for g in reversed(P): acc = add(mul(acc, s), const(g))
        return acc
    def inv(a):                             # a unit of R: inverse on the first digit by a^(L^d − 2), then Newton
        x = pw([v % L for v in a], L ** d - 2)          # correct modulo 607 (R/607 is the field with 607^d elements)
        for _ in range(3): x = add(mul(x, const(2)), [(-v) % M for v in mul(a, mul(x, x))])   # x·(2 − a·x)
        return x
    return d, mul, add, const, pw, at, inv


for m in re.finditer(r'a group of (\d+) roots.*?scaled carry: .*?\n\s+(\[\[.*?\]\])', txt, re.S):
    d = int(m.group(1)); co = ast.literal_eval(m.group(2))
    G = [sum(v * L ** i for i, v in enumerate(c)) % M for c in co]
    assert len(G) == d + 1 and G[-1] == 1
    d, mul, add, const, pw, at, inv = ring(G)
    dG = [(i * G[i]) % M for i in range(1, d + 1)]      # G'
    t = [0] * d; t[1 if d > 1 else 0] = 1
    strands = [t]
    for i in range(1, d):
        y = pw(strands[-1], L)                           # a root on the first digit
        for _ in range(3): y = add(y, [(-v) % M for v in mul(at(G, y), inv(at(dG, y)))])   # Newton: y − G(y)/G'(y)
        strands.append(y)
    ok_root = all(all(v == 0 for v in at(G, s)) for s in strands)
    distinct = len({tuple(s) for s in strands}) == d
    prod, sm = const(1), const(0)
    for s in strands: prod = mul(prod, s); sm = add(sm, s)
    scal = lambda a: a[0] if all(v == 0 for v in a[1:]) else None
    say(); say('the braid of %d, modulo 607^3: G = %s (coefficients lowest first, digits lowest first)' % (d, [dig(g) for g in G]))
    say('   every strand a root modulo 607^3: %s; all %d distinct: %s' % (ok_root, d, distinct))
    pv, sv = scal(prod), scal(sm)
    say('   the strands multiply to %s (expected %s) and add up to %s (expected %s); digits lowest first' % (
        dig(pv) if pv is not None else 'NOT an ordinary number', dig(((-1) ** d * G[0]) % M),
        dig(sv) if sv is not None else 'NOT an ordinary number', dig((-G[d - 1]) % M)))
    for i, s in enumerate(strands):
        say('   strand %d: %s' % (i + 1, [dig(v) for v in s]))
fh.close()
