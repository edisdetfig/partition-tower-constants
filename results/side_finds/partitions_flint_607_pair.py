# A dig, 12 Nov 2026: 607's TWO FORMS THAT SHARE THEIR NUMBER, ONE DIGIT DEEPER. At 607 the class (25 forms of weight
# 604 on Γ0(6)) has, on lane 607, a repeated root of T_607 with two independent eigenvectors, and T₅ tells the two forms
# apart (partitions_flint_607_span_*.txt). So they are different forms whose numbers at 607 agree on the first digit.
# PREDICTED before the run (derived, the generic case): on the two-dimensional piece T_607 = a·I + 607·C, and C has two
# DIFFERENT roots on lane 607: the two numbers separate at the second digit, and the pair is unramified (its discriminant
# has 607-adic valuation exactly 2), so the collision is an accident of one Hecke number, not a congruence of forms.
# Method (forms6_flint.py, FLINT only): the class at weight 604 on lane 607^4, separated exactly (class_padic); T_607 on it;
# the generalized eigenspace of the repeated root (kernel_part); C = (B − a·I)/607 there. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_607_pair.py [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat, nmod_poly
import forms6_flint as F
flint.ctx.threads = 8
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_flint_607_pair.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


t0 = time.time(); L = 607; k = L - 3; N = 4; m = L ** N
H = F.hecke_matrices(k, m, [('U', 2), ('U', 3), ('T', L), ('T', 5)])
W, RW, r, Wn, Rn, d, Np = F.class_padic(H, k, L, N, leg(2, L), leg(3, L))
mp = L ** Np
TL = F.restrict(F.mat_mod(k + 1, k + 1, [int(v) % mp for v in H[('T', L)].entries()], mp), Wn, Rn, d)
T5 = F.restrict(F.mat_mod(k + 1, k + 1, [int(v) % mp for v in H[('T', 5)].entries()], mp), Wn, Rn, d)
say('ℓ = 607, weight 604: the class has %d forms, known to 607^%d  (%.0f s)' % (d, Np, time.time() - t0))
fac = F.reduce_mat(TL, L).charpoly().factor()[1]
rep = [(f, e) for f, e in fac if e > 1]
say('   T_607 on lane 607: the repeated factors %s' % [([int(c) for c in f.coeffs()], e) for f, e in rep])
for f, e in rep:
    a = (-int(f.coeffs()[0])) % L                                  # the repeated root (a linear factor)
    Y, R, s = F.kernel_part(TL - F.identity(d, mp) * a, L, Np)   # the generalized eigenspace for roots ≡ a
    B = F.restrict(TL, Y, R, s); B5 = F.restrict(T5, Y, R, s)
    Bl = [[int(B[i, j]) for j in range(s)] for i in range(s)]
    semis = all((Bl[i][j] - (a if i == j else 0)) % L == 0 for i in range(s) for j in range(s))
    say('   the root %d: a piece of dimension %d; T_607 there is a·I on lane 607: %s' % (a, s, semis))
    if semis:
        mq = L ** (Np - 1)
        C = [[((Bl[i][j] - (a if i == j else 0)) // L) % mq for j in range(s)] for i in range(s)]
        chiC = F.berkowitz(C, mq)
        rootsC = nmod_poly([c % L for c in chiC], L).factor()[1]
        disc = (chiC[1] ** 2 - 4 * chiC[0]) % mq if s == 2 else None
        vd = next((e for e in range(Np - 1) if disc % L ** (e + 1)), Np - 1) if disc is not None else None
        say('   C = (T_607 − a)/607 there: characteristic polynomial on lane 607 %s, factors %s;  the pair\'s discriminant has '
            '607-adic valuation %s (2 + %s)' % ([c % L for c in chiC], [([int(x) for x in g.coeffs()], ee) for g, ee in rootsC],
                                               None if vd is None else 2 + vd, vd))
    say('   T₅ on the piece, on lane 607: characteristic polynomial %s' % [int(c) % L for c in F.berkowitz([[int(B5[i, j]) for j in range(s)] for i in range(s)], L)])
say('seconds: %.0f' % (time.time() - t0))
