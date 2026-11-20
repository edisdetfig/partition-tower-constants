# A dig, 13 Nov 2026: THE ONO-TYPE LAWS FOR 73 AND 61, CHECKED THROUGH THEIR STEP (README "Lucky primes for 73"; predicted
# there first). Their members need p(N) with N ≈ ℓ·Q³/8 (about 10¹⁹ to 10²⁰), out of reach. But an Ono-type law at Q means
# λ ≡ 0 for every form of the class, and then Statement B's step gives, for every M with ℓ | M, Q ∤ M, in both classes,
#     c(Q²M) ≡ ε·δ·c(M) = −(M|Q)·(−3|Q)·Q⁻²·c(M)  (mod ℓ),     c(M) = p((M + 1)/24),
# with N = (Q²M + 1)/24 of about 10¹⁴ to 10¹⁶. PREDICTED: all values as stated, for the two smallest Ono-type Q at 73 and
# at 61, two M of each class. FLINT's exact partition numbers, at idle priority. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_lucky_check_ono.py [--keep]
import os, sys, time
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
flint.ctx.threads = 4
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_flint_lucky_check_ono.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def jac(a, n):
    a %= n
    if a == 0: return 0
    return 1 if pow(a, (n - 1) // 2, n) == 1 else -1


def p_mod(N, L):
    return int(flint.fmpz(N).partitions_p()) % L


t0 = time.time(); right = total = 0
say('Ono-type laws (from partitions_flint_lucky_dual_73_61_max50000000.txt), checked through c(Q²M) ≡ −(M|Q)(−3|Q)Q⁻²·c(M)')
for L, Qs in ((73, (1109167, 7725269)), (61, (2144719, 2146489))):
    say(); say('ℓ = %d' % L)
    for Q in Qs:
        got = {1: [], -1: []}; tt = 1
        while len(got[1]) < 2 or len(got[-1]) < 2:
            M = L * tt
            if M % 24 == 23 and M % Q and len(got[jac(M, Q)]) < 2: got[jac(M, Q)].append(M)
            tt += 1
        for eps in (1, -1):
            for M in got[eps]:
                t = time.time(); N = (Q * Q * M + 1) // 24
                v = p_mod(N, L); c = p_mod((M + 1) // 24, L)
                pred = (-eps * jac(-3, Q) * pow(Q, -2, L) * c) % L; total += 1; right += v == pred
                say('   Q = %8d, M = %5d ((M|Q) = %+d): p(%d) ≡ %2d;  predicted −(M|Q)(−3|Q)Q⁻²·p(%d) = −(%+d)(%+d)·%d·%d = %2d  %s  (%.0f s)' % (
                    Q, M, eps, N, v, (M + 1) // 24, eps, jac(-3, Q), pow(Q, -2, L), c, pred, 'as PREDICTED' if v == pred else 'NOT as predicted', time.time() - t))
say(); say('%d of %d as predicted;  seconds: %.0f' % (right, total, time.time() - t0))
