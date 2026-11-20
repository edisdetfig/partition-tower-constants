# share_spectrum.py, 14 Nov 2026 (README "Why depth one? Every family starts as a unit"; predicted there first).
# The tower's depths at levels b = 1 … 4 against the class's slope spectrum at weight −2, if every family starts as a unit:
#   a system of slope s sits at depth k·s at level 2k and k·s + ⌈s/2⌉ at level 2k + 1 (s = 0: the loud notes, units).
# The spectrum: every seed Sym^r ⊗ det^m on the row (r + 2m ≡ ℓ − 5 mod ℓ − 1), with c(r + 2; ε) newforms of level 6 (signs
# ε from the twisted class numbers; first_step_w610.jsonl), each one system; its slope at −2 from slope_rule.py (the formula
# min(m, ℓ − 4 − m) where it was checked, else its ghost count; the Eisenstein-type seeds found at 13 … 23 by
# serre_weights.py with their own count). Predicted multiset on an S × S grid: the S smallest depths, padded with the cap.
#   --saved: share_test's saved invariant factors (full grids, levels 1 … 3, modulo ℓ⁵: cap 5) at the primes LMIN … LMAX;
#   --grid=S[,shift]: fresh S × S grids (rows and columns shifted by `shift` if given), levels 1 … 4, modulo ℓ⁶ (cap 6).
# Partition numbers exact (FLINT), workers at idle priority. COUNTED (depths), DERIVED (the spectrum).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/share_spectrum.py LMIN LMAX (--saved |
#          --grid=S[,shift]) [--only=p,q,…] [--procs=N] [--keep]
import os, sys, io, json, glob, time, contextlib
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import share_test as ST

_src = open(os.path.join(HERE, 'slope_rule.py'), encoding='utf-8').read()
_ns = {'__file__': os.path.join(HERE, 'slope_rule.py')}
_argv = sys.argv; sys.argv = ['slope_rule_defs']
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(_src[:_src.index("say('slopes at the tower")], 'slope_rule_defs', 'exec'), _ns)
sys.argv = _argv
family_slope = _ns['family_slope']
C = {}
for _line in open(os.path.join(HERE, 'first_step_w610.jsonl'), encoding='utf-8'):
    _d = json.loads(_line); C[_d['w']] = {tuple(int(x) for x in k.split(',')): v for k, v in _d['c'].items()}
EIS = {13: [(10, 11)], 17: [(14, 15)], 19: [(16, 17)], 23: [(20, 10), (20, 21)]}   # serre_weights_ℓ.txt, Eisenstein type


def spectrum(L):
    LM = L - 1; out = []
    for r in range(2, L, 2):
        for m in range(LM):
            if (r + 2 * m - (L - 5)) % LM: continue
            j = (L - 5 - r - 2 * m) // LM
            e2, e3 = ST.jacobi(2, L) ** ((1 + j) % 2), ST.jacobi(3, L) ** ((1 + j) % 2)
            mu = C.get(r + 2, {}).get((e2, e3), 0)
            if mu == 0: continue
            if (r, m) in EIS.get(L, []):
                s_eis = family_slope(L, (r, m), (m - 1) % LM)[2]
                out.append((s_eis, r, m, 'Eisenstein type')); mu -= 1
            if mu <= 0: continue
            s = min(m, L - 4 - m) if (r <= L - 3 and m <= L - 5) else family_slope(L, (r, m), None)[2]
            out += [(s, r, m, '')] * mu
    return sorted(out, key=lambda x: (x[0] if x[0] is not None else 999))


def depth(s, b):
    if s is None: return 999
    k, odd = divmod(b, 2)
    return k * s + (-(-s // 2) if odd else 0)


def predicted(sp, b, S, cap):
    d = sorted(min(depth(s, b), cap) for s, r, m, kind in sp)[:S]
    return d + [cap] * (S - len(d))


def main():
    nums = [int(a) for a in sys.argv[1:] if a.isdigit()]
    opt = {a.split('=')[0]: a.split('=')[1] for a in sys.argv[1:] if a.startswith('--') and '=' in a}
    KEEP, SAVED = '--keep' in sys.argv, '--saved' in sys.argv; procs = int(opt.get('--procs', 6))
    LMIN, LMAX = nums[0], nums[1]
    primes = [int(x) for x in opt['--only'].split(',')] if '--only' in opt else [p for p in range(LMIN, LMAX + 1) if ST.is_prime(p) and p >= 13]
    if SAVED: tag = 'saved_%d-%d' % (LMIN, LMAX)
    else:
        g = opt['--grid'].split(','); S = int(g[0]); shift = int(g[1]) if len(g) > 1 else 0
        tag = 'grid%d%s_%s' % (S, ('s%d' % shift) if shift else '', ('only_' + '_'.join(map(str, primes))) if '--only' in opt else '%d-%d' % (LMIN, LMAX))
    out = os.path.join(HERE, 'share_spectrum_%s.txt' % tag)
    if KEEP and os.path.exists(out): sys.exit('stopped before computing anything: %s exists' % out)
    fh = open(out, 'x', encoding='utf-8') if KEEP else None

    def say(*a):
        line = ' '.join(str(x) for x in a); print(line, flush=True)
        if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())

    t0 = time.time(); good = tot = 0; misses = []
    if SAVED:
        say('THE SAVED SHARE TESTS (full grids, modulo ℓ⁵) against the slope spectrum with unit starts, levels 1, 2, 3')
        rows = []
        for f in glob.glob(os.path.join(HERE, 'share_test_*.jsonl')):
            for line in open(f, encoding='utf-8'): rows.append(json.loads(line))
        for d in sorted(rows, key=lambda d: d['l']):
            L, S = d['l'], d['S']
            if L not in primes or S != d['R'] + d['K'] + 2: continue
            sp = spectrum(L); res = []
            for b in (1, 2, 3):
                pr = predicted(sp, b, S, 5); me = d['vals'][str(b)]; res.append((b, me == pr, me, pr))
            ok = all(x[1] for x in res); good += ok; tot += 1
            if not ok: misses.append(L)
            say('ℓ = %3d (grid %d): %s' % (L, S, '; '.join('level %d %s%s' % (b, 'as predicted' if okb else 'MEASURED %s, PREDICTED %s' % (me, pr), '') for b, okb, me, pr in res)))
        say(); say('as predicted at %d of %d primes; not: %s  (%.0f s)' % (good, tot, misses or 'none', time.time() - t0)); return

    say('FRESH %d × %d GRIDS%s, levels 1 … 4, modulo ℓ⁶, against the slope spectrum with unit starts' % (S, S, (' (rows and columns shifted by %d)' % shift) if shift else ''))
    M = 6
    with Pool(procs, initializer=ST.init_worker) as pool:
        for L in primes:
            t1 = time.time(); P = L ** M; sp = spectrum(L); res = []; n_args = 0
            qs = []; q = 5; k = 0
            while len(qs) < S:
                if ST.is_prime(q) and q != L:
                    if k >= shift: qs.append(q)
                    k += 1
                q += 1
            for b in (1, 2, 3, 4):
                cls = (23 * pow(L ** b, -1, 24)) % 24; ts = []; t = cls; k = 0
                while len(ts) < S:
                    if t % L and ST.squarefree(t):
                        if k >= shift: ts.append(t)
                        k += 1
                    t += 24
                args = sorted({(L ** b * t * Q * Q + 1) // 24 for t in ts for Q in qs} | {(L ** b * t + 1) // 24 for t in ts}, reverse=True)
                n_args += len(args)
                val = dict(pool.imap_unordered(ST.part_mod, [(n, P) for n in args], chunksize=2))
                c = lambda Mm: val[(Mm + 1) // 24]
                mat = [[(c(Q * Q * L ** b * t) + ST.jacobi(L ** b * t, Q) * ST.jacobi(-3 % Q, Q) * pow(Q * Q, -1, P) * c(L ** b * t)) % P
                        for Q in qs] for t in ts]
                me = ST.smith_vals(mat, L, M); pr = predicted(sp, b, S, M); res.append((b, me == pr, me, pr))
            ok = all(x[1] for x in res); good += ok; tot += 1
            if not ok: misses.append(L)
            say('ℓ = %3d: spectrum (slope, r, m) %s …' % (L, [(s, r, m) for s, r, m, kind in sp[:S]]))
            for b, okb, me, pr in res:
                say('   level %d: %s%s' % (b, me, '' if okb else '   PREDICTED %s' % pr))
            say('   -> %s  (%d partition numbers, %.0f s)' % ('as predicted' if ok else 'NOT AS PREDICTED', n_args, time.time() - t1))
    say(); say('as predicted at %d of %d primes; not: %s  (%.0f s)' % (good, tot, misses or 'none', time.time() - t0))


if __name__ == '__main__':
    main()
