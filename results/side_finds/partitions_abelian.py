# A dig, 11 Nov 2026: WHAT REPLACES THE CURVE FOR 17 TO 31. partitions_nocurve.txt proved that no elliptic curve over Q
# carries the towers' rule for 17, 19, 23, 29, 31. The natural replacement is a weight-2 newform g with trivial character whose
# coefficients lie in a number field K (an abelian variety over Q of dimension [K:Q]): its traces are bounded only conjugate by
# conjugate, so Hasse's bound does not block it. Here, with PARI (2.15.4 through cypari2, run in WSL), for every prime ℓ from
# 7 to 31:
#   1  the tower's newform f of level 6 (weight and a₂, a₃ as found in partitions_newforms.txt), recomputed by PARI, and the
#      target u_Q = Q^(−(k−2)/2)·a_Q(f) mod ℓ (partitions_nocurve.txt: a curve would need trace u_Q at every Q);
#   2  every newform of weight 2, trivial character, at the levels ℓ, 2ℓ, 3ℓ, 6ℓ, ℓ², 2ℓ², 3ℓ², 6ℓ² (where Serre's
#      conjecture puts such a form, the prime-to-ℓ conductor dividing 6), reduced at every prime 𝔩 of K above ℓ, against
#      (d|Q)·u_Q for every twist d in ±1, ±2, ±3, ±6 and the same times ℓ, at every prime Q up to 300 not dividing 6ℓ.
# For 7, 11, 13 this must find the curves 42a, 726i/726d, 1014b/1014f (dimension 1). Counts COUNTED; the field data are PARI's.
#   usage (WSL, with the cypari2 wheel unpacked in a folder): PYTHONPATH=<folder> python3 -B -s results/side_finds/partitions_abelian.py
#   [primes …] [--keep]   (with primes given, only those are done, and the .txt is named partitions_abelian_<primes>.txt)
import os, re, sys, time
os.nice(19)
import cypari2
pari = cypari2.Pari(); pari.allocatemem(4 * 10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
ONLY = [a for a in sys.argv[1:] if a.isdigit()]                 # optional: the primes to do (then the .txt is named after them)
fh = open(os.path.join(HERE, 'partitions_abelian%s.txt' % (('_' + '_'.join(ONLY)) if ONLY else '')), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); QMAX = 300
only = [int(a) for a in sys.argv[1:] if a.isdigit()]
txt = open(os.path.join(HERE, 'partitions_newforms.txt'), encoding='utf-8').read()
WHICH = {int(l): (int(k), int(a2), int(a3)) for l, k, a2, a3 in re.findall(
    r'ℓ = +(\d+), level 1 \(\d+ primes\): weight (\d+), a₂ = ([+-]\d+), a₃ = ([+-]\d+)', txt)}
say('PARI %s;  the towers\' newforms (partitions_newforms.txt): %s' % ('.'.join(map(str, pari.version())), ', '.join(
    'ℓ = %d: weight %d, a₂ = %+d, a₃ = %+d' % (l, *WHICH[l]) for l in sorted(WHICH))))
PR = [int(q) for q in pari('primes([5, %d])' % QMAX)]


def tower_form(l):
    k, a2, a3 = WHICH[l]; mf = pari('mfinit([6, %d], 0)' % k)
    for F in pari.mfeigenbasis(mf):
        c = pari.mfcoefs(F, QMAX)
        if int(c[2]) == a2 and int(c[3]) == a3: return k, [int(x) for x in c]
    raise ValueError('no newform with these a2, a3')


for l in sorted(WHICH):
    if only and l not in only: continue
    k, f = tower_form(l); b0 = (k - 2) // 2; Qs = [Q for Q in PR if Q != l]
    u = {Q: f[Q] * pow(Q, -b0, l) % l for Q in Qs}
    say(); say('ℓ = %d: the tower\'s newform has weight %d (PARI: q %+d·q² %+d·q³ %+d·q⁴ %+d·q⁵ %+d·q⁷ …); target u_Q = Q^-%d·a_Q(f) mod %d' % (
        l, k, f[2], f[3], f[4], f[5], f[7], b0, l))
    found = []
    for N in (l, 2 * l, 3 * l, 6 * l, l * l, 2 * l * l, 3 * l * l, 6 * l * l):
        t = time.time(); mf = pari('mfinit([%d, 2], 0)' % N); dim = int(pari.mfdim(mf))
        if dim == 0: continue
        B = pari.mfeigenbasis(mf); flds = pari.mffields(mf); degs = [int(pari.poldegree(P)) for P in flds]
        for i, F in enumerate(B):
            P = flds[i]; c = pari.mfcoefs(F, QMAX)
            facs = pari.factormod(P, l)[0]                                  # the primes of K above ℓ (by their factors)
            for g in facs:
                res = {}
                ok_g = True
                for Q in Qs:
                    try:
                        r = pari('(lift(%s) * Mod(1, %d)) %% (%s)' % (c[Q], l, g)) if degs[i] > 1 else pari('Mod(%s, %d)' % (c[Q], l))
                    except cypari2.PariError:
                        ok_g = False; break
                    res[Q] = r
                if not ok_g: continue
                for d in (1, -1, 2, -2, 3, -3, 6, -6, l, -l, 2 * l, -2 * l, 3 * l, -3 * l, 6 * l, -6 * l):
                    good = True
                    for Q in Qs:
                        tQ = int(pari.kronecker(d, Q)) * u[Q] % l
                        if degs[i] > 1:
                            if pari('(%s) - Mod(%d, %d)' % (res[Q], tQ, l)) != 0: good = False; break
                        elif int(pari.lift(res[Q])) != tQ: good = False; break
                    if good: found.append((N, i + 1, degs[i], int(pari.poldegree(g)), d, str(P)))
        say('   level %5d: %3d newforms in %3d Galois orbits (degrees %s)  (%.0f s)' % (
            N, dim, len(degs), ' '.join(map(str, sorted(degs))) if len(degs) <= 25 else '%s … %s' % (' '.join(map(str, sorted(degs)[:12])), max(degs)), time.time() - t))
    if not found: say('   no newform of weight 2 at these levels has the tower\'s numbers on lane %d (twists included)' % l)
    for N, i, dg, rd, d, P in found:
        say('   MATCH: level %d, orbit %d, Hecke field of degree %d (%s), at a prime above %d of residue degree %d, twist d = %d' % (
            N, i, dg, P if len(P) < 80 else P[:77] + '…', l, rd, d))
say(); say('seconds: %.0f' % (time.time() - t0))
