# sign_check.py, 14 Nov 2026 (README "Do the families meet the conditions of the proofs?"; a relayed review: "Within one
# family, the signs at 2 and 3 are fixed by ρ̄ ... the class can be recovered by subtraction").
# The check behind that argument: every system (family) of forms new at 2 and 3 should occur in only ONE of the four
# sign classes (U_2, U_3 = ±2^((k−2)/2), ±3^((k−2)/2)), at level 6ℓ and at level 6. Reason (DERIVED): for a form new at
# 2, the Frobenius eigenvalues of ρ̄ at 2 are {α, 2α} with α fixed by the U_2 sign, and {α, 2α} = {β, 2β} with β ≠ α
# would need 3α ≡ 0, impossible for ℓ ≠ 3; the same at 3 with {α, 3α} (8α ≡ 0, impossible for ℓ ≠ 2).
# The code is level6l_counts2.py's own (its functions are read from that file and run unchanged), with every sign
# class split into systems instead of only the class. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/sign_check.py ℓ x1 x2 … [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); SRC = os.path.join(HERE, 'level6l_counts2.py')
args = [a for a in sys.argv[1:] if a.isdigit()]; KEEPF = '--keep' in sys.argv
L0, XS0 = int(args[0]), [int(a) for a in args[1:]]
code = open(SRC, encoding='utf-8').read().split('# ---- the run ----')[0]
sys.argv = [SRC, str(L0), str(max(XS0))]                    # no --keep: level6l_counts2.py's part writes no files
g = {'__file__': SRC, '__name__': 'sign_check'}
exec(compile(code, SRC, 'exec'), g)
L, F, QS = g['L'], g['F'], g['QS']
OUT = os.path.join(HERE, 'sign_check_%d.txt' % L)
if KEEPF and os.path.exists(OUT): sys.exit('stopped before computing anything: %s exists' % OUT)
fh = open(OUT, 'x', encoding='utf-8') if KEEPF else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


say('ℓ = %d: the systems of each of the four sign classes (forms new at 2 and 3), at level 6ℓ and level 6, x = %s; '
    'each system should occur in one class only' % (L, XS0))
bad_total = 0; t0 = time.time()
for x in XS0:
    k = -2 + (L - 1) * (x + 1)
    if k < 4: continue
    S6l, c6l = g['dim_S'](6 * L, k); n = S6l + c6l; B = k * (L + 1); Lx = max(QS + [2, 3]) * B + 1
    chosen, tried = g['select_Q'](k, n, B)
    assert len(chosen) == n
    Sr = g['Series'](k, Lx, L ** g['MS']); rows = [Sr.form(d) for d in chosen]; del Sr
    rows, rounds, total = g['saturate'](rows, B, Lx, L ** g['MS'])
    Np = min(g['NP'], g['MS'] - rounds)
    H, ok, comm = g['operators'](rows, B, k, Np, Lx); del rows
    assert ok and comm
    H6 = F.hecke_matrices(k, L ** g['NP6'], [('T', q) for q in QS] + [('U', 2), ('U', 3)])
    for lev, HH, nn, Nn in (('6ℓ', H, n, Np), ('6', H6, k + 1, g['NP6'])):
        where = {}
        for s2 in (1, -1):
            p2 = g['u2_part'](HH, nn, k, s2, Nn)
            for s3 in (1, -1):
                d, N1, sy, my, rj = g['the_class'](HH, k, p2, s3, Nn, systems=True)
                for lab, (cnt, _) in (sy or {}).items(): where.setdefault(lab, []).append(((s2, s3), cnt))
        multi = {lab: v for lab, v in where.items() if len(v) > 1}
        bad_total += len(multi)
        say('   x = %d (weight %d), level %s: %d systems over the four sign classes; in more than one class: %s' % (
            x, k, lev, len(where), ', '.join('%s %s' % (lab, v) for lab, v in multi.items()) or 'none'))
say(); say('systems found in more than one sign class: %d  (%.0f s)' % (bad_total, time.time() - t0))
