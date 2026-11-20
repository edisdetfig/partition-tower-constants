# share_test.py, 14 Nov 2026 (README "The share condition, tested on partition numbers"; predicted there first).
# Does the tower carry a share of the first-step (slope-1) families? No forms are used, only partition numbers.
#   c(M) = p((M + 1)/24). Rows: the first S squarefree t > 0, ℓ ∤ t, with ℓᵇt ≡ 23 (mod 24). Columns: the first S primes
#   Q ≥ 5, Q ≠ ℓ. Entry M_b[t, Q] = c(Q²·ℓᵇt) + (ℓᵇt|Q)·(−3|Q)·Q⁻²·c(ℓᵇt)  (the T(Q²) rule of partitions_split.py).
#   Each family i of the tower adds λ_i(Q)·c_i(ℓᵇt), so M_b = (shares by t) × (eigenvalues by Q). Its invariant factors over
#   ℤ/ℓ⁵ (Smith form; 5 printed for a factor ≡ 0) sort the families by the depth of their share at level b.
# Levels b = 1, 2, 3. Partition numbers exact (FLINT), each worker at idle priority with one thread. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/share_test.py LMIN LMAX
#          [--only=p,q,...] [--small-from=X (S = R + 3 for ℓ ≥ X)] [--procs=N] [--plan] [--keep]
import os, sys, time, json
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
M5 = 5                                                       # work modulo ℓ^5


def is_prime(n): return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def squarefree(n):
    d = 2
    while d * d <= n:
        if n % (d * d) == 0: return False
        d += 1
    return True


def jacobi(a, n):                                            # Jacobi symbol (a|n), n odd positive
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0


def vl(x, L, cap):
    if x == 0: return cap
    v = 0
    while x % L == 0 and v < cap: x //= L; v += 1
    return v


def smith_vals(A, L, m):
    """Valuations of the invariant factors of the integer matrix A over ℤ/ℓ^m (m for a factor ≡ 0), smallest first."""
    P = L ** m; A = [[x % P for x in row] for row in A]; nr, nc = len(A), len(A[0]); out = []
    for k in range(min(nr, nc)):
        best = None
        for i in range(k, nr):
            for j in range(k, nc):
                if A[i][j]:
                    v = vl(A[i][j], L, m)
                    if best is None or v < best[0]: best = (v, i, j)
                    if v == 0: break
            if best and best[0] == 0: break
        if best is None: out += [m] * (min(nr, nc) - k); break
        v, i, j = best
        A[k], A[i] = A[i], A[k]
        for row in A: row[k], row[j] = row[j], row[k]
        u = (A[k][k] // L ** v) % P; ui = pow(u, -1, P)
        for i2 in range(k + 1, nr):
            if A[i2][k]:
                f = (A[i2][k] // L ** v) * ui % P
                A[i2] = [(a - f * b) % P for a, b in zip(A[i2], A[k])]
        out.append(v)
    return out


def grid(L, S):
    tq = {}
    for b in (1, 2, 3):
        cls = (23 * pow(L ** b, -1, 24)) % 24; ts = []; t = cls
        while len(ts) < S:
            if t % L and squarefree(t): ts.append(t)
            t += 24
        tq[b] = ts
    qs = []; q = 5
    while len(qs) < S:
        if is_prime(q) and q != L: qs.append(q)
        q += 1
    return tq, qs


def needed(L, S):
    tq, qs = grid(L, S); args = set()
    for b in (1, 2, 3):
        for t in tq[b]:
            M0 = L ** b * t; args.add((M0 + 1) // 24)
            for Q in qs: args.add((Q * Q * M0 + 1) // 24)
    return tq, qs, sorted(args, reverse=True)


def est_seconds(n): return 0.02 if n < 10 ** 9 else 2.0 * (n / 1e12) ** 0.6


def init_worker():
    try:
        import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
    except Exception: pass
    import flint; flint.ctx.threads = 1


def part_mod(job):
    import flint
    n, P = job
    return n, int(flint.fmpz(n).partitions_p() % P)


def main():
    nums = [int(a) for a in sys.argv[1:] if a.isdigit()]
    opt = {a.split('=')[0]: a.split('=')[1] for a in sys.argv[1:] if a.startswith('--') and '=' in a}
    PLAN, KEEP = '--plan' in sys.argv, '--keep' in sys.argv
    if '--only' in opt: primes = [int(x) for x in opt['--only'].split(',')]; tag = 'only_' + '_'.join(map(str, primes))
    else: primes = [p for p in range(nums[0], nums[1] + 1) if is_prime(p) and p >= 13]; tag = '%d-%d' % (nums[0], nums[1])
    small_from = int(opt.get('--small-from', 10 ** 9)); procs = int(opt.get('--procs', 8))
    out_txt = os.path.join(HERE, 'share_test_%s.txt' % tag); out_js = os.path.join(HERE, 'share_test_%s.jsonl' % tag)
    if KEEP and not PLAN:
        for f in (out_txt, out_js):
            if os.path.exists(f): sys.exit('stopped before computing anything: %s exists' % f)
    fh = open(out_txt, 'x', encoding='utf-8') if KEEP and not PLAN else None
    fj = open(out_js, 'x', encoding='utf-8') if KEEP and not PLAN else None

    def say(*a):
        line = ' '.join(str(x) for x in a); print(line, flush=True)
        if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())

    def grid_size(L):
        R, K = (L + 12) // 24, L // 24
        return R, K, (R + 3 if L >= small_from else R + K + 2)

    if PLAN:
        tot = 0
        for L in primes:
            R, K, S = grid_size(L); tq, qs, args = needed(L, S); e = sum(est_seconds(n) for n in args); tot += e
            print('ℓ = %3d: R = %2d, K = %2d, grid %2d × %2d, %4d partition numbers, the largest p(%.2e); about %.0f s on one core'
                  % (L, R, K, S, S, len(args), args[0], e))
        print('total about %.0f s on one core, %.0f s on %d' % (tot, tot / procs, procs)); return

    say('THE SHARE TEST: invariant factors of M_b[t, Q] = c(Q²·ℓᵇt) + (ℓᵇt|Q)(−3|Q)Q⁻²·c(ℓᵇt) over ℤ/ℓ⁵, levels b = 1, 2, 3')
    say('   (valuations, smallest first; 5 = a factor ≡ 0 mod ℓ⁵). R = ⌊(ℓ+12)/24⌋ loud notes, K = ⌊ℓ/24⌋ first-step forms.')
    say('   Predicted: level 1 and level 2 have exactly R units (607: R − 1); level 3 has R units, then K factors ℓ² (S = R + K + 2)')
    say('   or three (S = R + 3), then factors divisible by ℓ³.')
    t0 = time.time(); summary = []
    with Pool(procs, initializer=init_worker) as pool:
        for L in primes:
            R, K, S = grid_size(L); tq, qs, args = needed(L, S); P = L ** M5; t1 = time.time()
            val = dict(pool.imap_unordered(part_mod, [(n, P) for n in args], chunksize=1))
            c = lambda M: val[(M + 1) // 24]
            vals = {}
            for b in (1, 2, 3):
                rows = []
                for t in tq[b]:
                    M0 = L ** b * t; row = []
                    for Q in qs:
                        e = jacobi(M0, Q) * jacobi(-3 % Q, Q)
                        row.append((c(Q * Q * M0) + e * pow(Q * Q, -1, P) * c(M0)) % P)
                    rows.append(row)
                vals[b] = smith_vals(rows, L, M5)
            u = {b: sum(1 for v in vals[b] if v == 0) for b in (1, 2, 3)}
            two = sum(1 for v in vals[3] if v == 2); one3 = sum(1 for v in vals[3] if v == 1)
            Rexp = R - 1 if L == 607 else R
            k_exp = min(K, S - Rexp)                         # K on the full grid; on S = R + 3, three (607: four)
            ok = (u[1] == Rexp, u[2] == Rexp, u[3] == Rexp and one3 == 0 and two == k_exp)
            dt = time.time() - t1
            say('ℓ = %3d: R = %d, K = %d, grid %d × %d (Q = %d … %d); level 1: %s; level 2: %s; level 3: %s  ->  units %d, %d, %d '
                '(predicted %d); level 3 factors ℓ²: %d (predicted %d)%s  (%d partition numbers up to p(%.2e), %.0f s)' % (
                    L, R, K, S, S, qs[0], qs[-1], vals[1], vals[2], vals[3], u[1], u[2], u[3], Rexp, two, k_exp,
                    '' if all(ok) else '  <-- NOT AS PREDICTED', len(args), args[0], dt))
            summary.append((L, ok, two, k_exp))
            if fj:
                fj.write(json.dumps({'l': L, 'R': R, 'K': K, 'S': S, 't': tq, 'Q': qs, 'vals': vals, 'units': u, 'l2_factors': two,
                                     'predicted_units': Rexp, 'predicted_l2': k_exp, 'seconds': round(dt, 1)}) + '\n')
                fj.flush(); os.fsync(fj.fileno())
    say(); say('as predicted at %d of %d primes; not: %s  (%.0f s)' % (
        sum(all(s[1]) for s in summary), len(summary), [s[0] for s in summary if not all(s[1])] or 'none', time.time() - t0))


if __name__ == '__main__':
    main()
