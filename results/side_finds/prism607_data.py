# A dig, 13 Nov 2026: 607 THROUGH THE PRISM AT 5 (README "607 through the prism"). The class at weight 604 modulo 607 (forms6_flint.class_padic, as in
# chords_page_data607.py), split by T₆₀₇ into its pieces on lane 607 (the lone note, the pair, the braid of 22); on each piece
# T₅'s characteristic polynomial mod 607 and its factors; whether the pieces' T₅ polynomials are pairwise coprime (then the
# step at 5 tells the pieces apart on lane 607); whether T₅ on the whole class has 25 distinct roots (then it tells every
# form apart, in the bigger fields); and T₆₀₇ on the pair (≡ 367·I: the tower's own step cannot tell the pair apart on the
# first digit). Writes prism607_data.json. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/prism607_data.py [--keep]
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly
import forms6_flint as F
flint.ctx.threads = 4
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'prism607_data.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


t0 = time.time(); L = 607; k = L - 3; N = 2; m = L ** N; s3 = leg(3, L)
H = F.hecke_matrices(k, m, [('U', 2), ('U', 3), ('T', L), ('T', 5)])
W, RW, r, Wn, Rn, d, Np = F.class_padic(H, k, L, N, leg(2, L), leg(3, L))
mp = L ** Np
TL = F.restrict(F.mat_mod(k + 1, k + 1, [int(v) % mp for v in H[('T', L)].entries()], mp), Wn, Rn, d)
T5 = F.restrict(F.mat_mod(k + 1, k + 1, [int(v) % mp for v in H[('T', 5)].entries()], mp), Wn, Rn, d)
say('ℓ = 607, weight 604: the class has %d forms (to 607^%d)  (%.0f s)' % (d, Np, time.time() - t0))
TLb, T5b = F.reduce_mat(TL, L), F.reduce_mat(T5, L)
pieces = []
for f, e in TLb.charpoly().factor()[1]:
    co = [int(c) for c in f.coeffs()]
    B = TL - F.identity(d, mp) * ((-co[0]) % L) if len(co) == 2 else F.mat_poly(co, TL)
    Y, R, s = F.kernel_part(B, L, Np)
    T5p = F.reduce_mat(F.restrict(T5, Y, R, s), L); TLp = F.reduce_mat(F.restrict(TL, Y, R, s), L)
    chi5 = T5p.charpoly(); fac5 = chi5.factor()[1]
    kind = 'the lone note' if s == 1 else ('the pair' if s == 2 else 'the braid of %d' % s)
    scalar = all(int(TLp[i, j]) == ((-co[0]) % L if i == j else 0) for i in range(s) for j in range(s)) if len(co) == 2 else False
    const = s3 * (-co[0]) % L if len(co) == 2 else None
    pieces.append({'kind': kind, 'dim': s, 'T607_factor': co, 'const': const, 'T607_scalar': scalar,
                   'T5_poly': [int(c) for c in chi5.coeffs()], 'T5_factors': [([int(c) for c in g.coeffs()], ee) for g, ee in fac5], 'chi5': chi5})
    say('   %s: %d form%s; T₆₀₇ factor %s%s; constant %s; T₅ there: polynomial of degree %d, factors of degrees %s' % (
        kind, s, '' if s == 1 else 's', co if len(co) <= 3 else 'of degree %d' % (len(co) - 1), ' (T₆₀₇ = 367·I there on lane 607)' if scalar and s > 1 else '',
        const if const is not None else 'in the field with 607^%d elements' % (len(co) - 1), chi5.degree(), sorted(g.degree() for g, ee in fac5 for _ in range(ee))))
cop = all(pieces[i]['chi5'].gcd(pieces[j]['chi5']).degree() == 0 for i in range(len(pieces)) for j in range(i + 1, len(pieces)))
whole5 = T5b.charpoly(); sqf = whole5.gcd(whole5.derivative()).degree() == 0
say('the pieces\' T₅ polynomials pairwise coprime on lane 607 (the step at 5 tells the pieces apart): %s' % cop)
say('T₅ on the whole class: 25 distinct roots in the bigger fields (the step at 5 tells every form apart): %s' % sqf)
lone = next(p for p in pieces if p['dim'] == 1)
say('the lone note\'s number at 5: a₅ ≡ %d, so (3|5)·a₅ ≡ %d (its multiplier at 5)' % ((-lone['T5_poly'][0]) % L, leg(3, 5) * (-lone['T5_poly'][0]) % L))
OUT = {'pieces': [{kk: vv for kk, vv in p.items() if kk != 'chi5'} for p in pieces], 'coprime': cop, 'T5_squarefree': sqf,
       'lone_a5': (-lone['T5_poly'][0]) % L}
if KEEP: json.dump(OUT, open(os.path.join(HERE, 'prism607_data.json'), 'x'), separators=(',', ':')); say('wrote prism607_data.json')
say('seconds: %.0f' % (time.time() - t0))
