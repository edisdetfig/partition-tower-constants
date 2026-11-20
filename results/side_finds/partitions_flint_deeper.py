# A dig, 12 Nov 2026: DEEPER CHECKS WITH FLINT. FLINT computes p(n) exactly (partitions_flint_arbiter.txt) and fast:
# 3 s at 10¹², 39 s at 10¹⁴, 172 s at 10¹⁵ (measured). So the partition-number side of three findings is pushed deeper.
#   1  THE LAWS: 10 members each of 37's laws at Q = 6599, 7541, 9547 (Ahlgren–Allen–Tang; class (n|Q) as there) and its
#      Ono-type law at Q = 2999 (24N − 1 = 37·2999³·m, 2999 ∤ m), and of 43's laws at Q = 4049, 8941 (partitions_split.txt).
#      PREDICTED: every member is a multiple of ℓ.
#   2  THE TWINS ONE LEVEL DEEPER: H_(b+4) − T·H_(b+2) + N·H_b with H_b = p(δ_b), b up to 9, for 37, 41, 43 (T, N from the
#      twins' constants in partitions_twins.txt and partitions_twins_43.txt). PREDICTED: lane ℓ³ at levels 5, 7, 9 (each
#      two levels add a factor, as derived in partitions_twins.py), lane ℓ² at least at 4, 6, 8.
#   3  THE CONSTANTS OF 13 TO 31: the ratio H_(b+2)(0)/H_b(0) for every pair of levels with δ_(b+2) ≤ 10¹⁷, against the
#      forms' digits (no partition number used there): 13 digits for 13, 17, 19, 23 (partitions_hida_digits_x12.txt),
#      11 for 29, 31 (partitions_hida_digits.txt). PREDICTED: every digit the partition numbers settle agrees with the
#      forms' digit (for 29 that includes the third digit, 17, which our data never reached and which Boylan–Webb's
#      14706 modulo 29³ carries).
# FLINT runs on 8 threads (5 times faster at 10¹⁴, measured). A first single-thread run of this script was stopped and its
# partial file replaced; it had reached the same results through 13's levels (11, 13). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_deeper.py [--keep]
import os, sys, time, re
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
import flint
from sympy import jacobi_symbol
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_flint_deeper.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); PC = {}; LIMIT = 10 ** 17
flint.ctx.threads = 8
leg = lambda a, q: int(jacobi_symbol(a % q, q))
read = lambda name: open(os.path.join(HERE, name), encoding='utf-8').read()


def P(n):
    if n not in PC: PC[n] = int(flint.fmpz(n).partitions_p())
    return PC[n]


# ---- 1. the laws ----
say('1  THE LAWS (10 members each, FLINT)')
for L, Q, e in ((37, 6599, -1), (37, 7541, 1), (37, 9547, -1), (43, 4049, -1), (43, 8941, -1)):
    d1 = pow(24, -1, L); j = 0; got = []; t = time.time()
    while len(got) < 10:
        Ns = L * j + d1; M = 24 * Ns - 1; j += 1
        if M % Q == 0 or leg(L, Q) * leg(M, Q) != e: continue
        N = Q * Q * Ns - (Q * Q - 1) // 24; got.append((N, P(N) % L))
    say('   ℓ = %d, Q = %d, (n|Q) = %+d: %d of 10 members a multiple of %d (largest N = %d; %.0f s)' % (
        L, Q, e, sum(1 for N, v in got if v == 0), L, max(N for N, v in got), time.time() - t))
ms = [m for m in range(1, 1000) if (37 * 2999 ** 3 * m + 1) % 24 == 0 and m % 2999][:10]; t = time.time()
vals = [P((37 * 2999 ** 3 * m + 1) // 24) % 37 for m in ms]
say('   ℓ = 37, Q = 2999 (24N − 1 = 37·2999³·m, m = %s): %d of 10 a multiple of 37 (largest N = %d; %.0f s)' % (
    ms, sum(1 for v in vals if v == 0), (37 * 2999 ** 3 * ms[-1] + 1) // 24, time.time() - t))

# ---- 2. the twins one level deeper ----
say(); say('2  THE TWINS ONE LEVEL DEEPER (FLINT)')
for L, name in ((37, 'partitions_twins.txt'), (41, 'partitions_twins.txt'), (43, 'partitions_twins_43.txt')):
    text = read(name); block = text[text.index('ℓ = %d, twin 1' % L):]
    digs = [[int(d) for d in s.split(',')] for s in re.findall(r'digits lowest first: \[([^\]]*)\]', block)[:2]]
    K = min(len(d) for d in digs); MK = L ** K; t = time.time()
    a1, a2 = (sum(d[i] * L ** i for i in range(K)) for d in digs); T, N = (a1 + a2) % MK, a1 * a2 % MK
    H = {b: P(pow(24, -1, L ** b)) % MK for b in range(1, 10)}
    val = lambda x: next((e for e in range(K) if x % L ** (e + 1)), K)
    say('   ℓ = %d (constants to %d digits): divisible by ℓ^ %s at the levels %s  (%.0f s)' % (
        L, K, [val((H[b + 4] - T * H[b + 2] + N * H[b]) % MK) for b in range(1, 6)], ['(%d,%d,%d)' % (b, b + 2, b + 4) for b in range(1, 6)], time.time() - t))

# ---- 3. the constants of 13 to 31 ----
say(); say('3  THE CONSTANTS OF 13 TO 31 (FLINT; pairs of levels two apart, n = 0, δ ≤ 10¹⁷)')
for L in (13, 17, 19, 23, 29, 31):
    text = read('partitions_hida_digits_x12.txt' if L <= 23 else 'partitions_hida_digits.txt')
    block = text[text.index('ℓ = %d:' % L):]
    forms = [int(d) for d in re.search(r'u\(−2\), digits lowest first: \[([^\]]*)\]', block).group(1).split(',')]
    K = len(forms); M = L ** K; H = {}; prev = None; rows = []
    say('   ℓ = %d, the forms\' digits: %s' % (L, forms))
    b = 1
    while pow(24, -1, L ** b) <= LIMIT:
        t = time.time(); H[b] = P(pow(24, -1, L ** b)) % M; dt = time.time() - t
        if b - 2 in H:
            a = b - 2
            if H[a] % L == 0: say('      levels %d, %d: H_%d(0) divisible by ℓ, skipped' % (a, b, a))
            else:
                r = H[b] * pow(H[a], -1, M) % M; dg = [(r // L ** i) % L for i in range(K)]
                same = next((i for i in range(K) if dg[i] != forms[i]), K)
                self_ = next((i for i in range(K) if prev is not None and dg[i] != prev[i]), K) if prev is not None else None
                say('      levels %d, %d (δ = %d, %.0f s): digits %s;  agree with the forms in %d%s' % (
                    a, b, pow(24, -1, L ** b), dt, dg, same, '' if self_ is None else ';  with the previous pair in %d' % self_))
                prev = dg
        b += 1
say(); say('seconds: %.0f' % (time.time() - t0))
