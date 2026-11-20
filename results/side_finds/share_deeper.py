# share_deeper.py, 14 Nov 2026 (README "The share condition, tested on partition numbers"; predicted there first).
# share_test.py's matrices one and two levels deeper: M_b[t, Q] = c(Q²·ℓᵇt) + (ℓᵇt|Q)(−3|Q)Q⁻²·c(ℓᵇt), c(M) = p((M+1)/24),
# on the same grid (S = R + K + 2), for b = 4 and 5, over ℤ/ℓ⁶. Predicted: level 4 has R units, K factors ℓ², the rest
# divisible by ℓ⁴; level 5 has R units, K factors ℓ³, the rest divisible by ℓ⁵. Partition numbers exact (FLINT), workers at
# idle priority. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/share_deeper.py LMIN LMAX [--procs=N]
#          [--plan] [--keep]
import os, sys, time, json
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import share_test as ST
M6 = 6
LEVELS = (4, 5)


def grid_b(L, S, b):
    cls = (23 * pow(L ** b, -1, 24)) % 24; ts = []; t = cls
    while len(ts) < S:
        if t % L and ST.squarefree(t): ts.append(t)
        t += 24
    return ts


def main():
    nums = [int(a) for a in sys.argv[1:] if a.isdigit()]
    opt = {a.split('=')[0]: a.split('=')[1] for a in sys.argv[1:] if a.startswith('--') and '=' in a}
    PLAN, KEEP = '--plan' in sys.argv, '--keep' in sys.argv; procs = int(opt.get('--procs', 4))
    primes = [p for p in range(nums[0], nums[1] + 1) if ST.is_prime(p) and p >= 13]
    out_txt = os.path.join(HERE, 'share_deeper_%d-%d.txt' % (nums[0], nums[1])); out_js = out_txt[:-4] + '.jsonl'
    if KEEP and not PLAN:
        for f in (out_txt, out_js):
            if os.path.exists(f): sys.exit('stopped before computing anything: %s exists' % f)
    fh = open(out_txt, 'x', encoding='utf-8') if KEEP and not PLAN else None
    fj = open(out_js, 'x', encoding='utf-8') if KEEP and not PLAN else None

    def say(*a):
        line = ' '.join(str(x) for x in a); print(line, flush=True)
        if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())

    jobs = []
    for L in primes:
        R, K = (L + 12) // 24, L // 24; S = R + K + 2
        tq0, qs = ST.grid(L, S); tq = {b: grid_b(L, S, b) for b in LEVELS}
        args = set()
        for b in LEVELS:
            for t in tq[b]:
                M0 = L ** b * t; args.add((M0 + 1) // 24)
                for Q in qs: args.add((Q * Q * M0 + 1) // 24)
        jobs.append((L, R, K, S, tq, qs, sorted(args, reverse=True)))
    if PLAN:
        tot = 0
        for L, R, K, S, tq, qs, args in jobs:
            e = sum(ST.est_seconds(n) for n in args); tot += e
            print('ℓ = %3d: grid %2d × %2d, %4d partition numbers, the largest p(%.2e); about %.0f s on one core' % (L, S, S, len(args), args[0], e))
        print('total about %.0f s on one core, %.0f s on %d' % (tot, tot / procs, procs)); return

    say('THE SHARE TEST, LEVELS 4 AND 5: invariant factors over ℤ/ℓ⁶ (valuations, smallest first; 6 = a factor ≡ 0), grid S = R + K + 2')
    say('   Predicted: level 4: R units, K factors ℓ², the rest ≥ ℓ⁴; level 5: R units, K factors ℓ³, the rest ≥ ℓ⁵.')
    t0 = time.time(); good = 0
    with Pool(procs, initializer=ST.init_worker) as pool:
        for L, R, K, S, tq, qs, args in jobs:
            t1 = time.time(); P = L ** M6
            val = dict(pool.imap_unordered(ST.part_mod, [(n, P) for n in args], chunksize=1))
            c = lambda M: val[(M + 1) // 24]
            vals = {}
            for b in LEVELS:
                rows = []
                for t in tq[b]:
                    M0 = L ** b * t
                    rows.append([(c(Q * Q * M0) + ST.jacobi(M0, Q) * ST.jacobi(-3 % Q, Q) * pow(Q * Q, -1, P) * c(M0)) % P for Q in qs])
                vals[b] = ST.smith_vals(rows, L, M6)
            Rexp = R - 1 if L == 607 else R
            ok4 = vals[4].count(0) == Rexp and vals[4].count(2) == K and all(v >= 4 for v in vals[4] if v not in (0, 2)) and vals[4].count(1) == 0 and vals[4].count(3) == 0
            ok5 = vals[5].count(0) == Rexp and vals[5].count(3) == K and all(v >= 5 for v in vals[5] if v not in (0, 3)) and all(vals[5].count(j) == 0 for j in (1, 2, 4))
            good += ok4 and ok5
            say('ℓ = %3d: R = %d, K = %d, grid %d × %d; level 4: %s; level 5: %s  ->  %s  (%d partition numbers up to p(%.2e), %.0f s)' % (
                L, R, K, S, S, vals[4], vals[5], 'as predicted' if ok4 and ok5 else 'NOT AS PREDICTED', len(args), args[0], time.time() - t1))
            if fj:
                fj.write(json.dumps({'l': L, 'R': R, 'K': K, 'S': S, 't': tq, 'Q': qs, 'vals': vals, 'ok': [ok4, ok5],
                                     'seconds': round(time.time() - t1, 1)}) + '\n'); fj.flush(); os.fsync(fj.fileno())
    say(); say('as predicted at %d of %d primes  (%.0f s)' % (good, len(jobs), time.time() - t0))


if __name__ == '__main__':
    main()
