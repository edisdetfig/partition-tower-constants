# A dig, 13 Nov 2026: THE NEW LAWS FOR 73 (and 61), CHECKED ON EXACT PARTITION NUMBERS (README "Lucky primes for 73";
# predicted there first). partitions_flint_lucky_dual_73_61_max50000000.txt lists, from the forms alone, every prime
# Q < 5·10⁷ where (3|Q)·T_Q is one number on the class modulo ℓ. A class law with ε says: p(N) ≡ 0 (mod ℓ) whenever
# 24N − 1 = Q²·M, ℓ | M, Q ∤ M, (M|Q) = ε.
# PREDICTED before the run, for the two smallest Q of each class kind:
#   - the two smallest members are ≡ 0 (mod ℓ);
#   - a control M′ of the other class at the same Q is not 0 but EXACTLY 2ε·(−3|Q)·Q⁻²·p((M′ + 1)/24) (mod ℓ): Statement B's
#     step c(Q²M′) = (λ + ε′δ)·c(M′) with λ = −εδ (the law's own number), ε′ = −ε, δ = −(−3|Q)·Q⁻² (derived).
# Members are computed only up to N = 3·10¹⁶ (FLINT's exact p(N), at idle priority). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_lucky_check_73.py [--keep]
import os, sys, re, time
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
flint.ctx.threads = 6
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRC = 'partitions_flint_lucky_dual_73_61_max50000000.txt'
fh = open(os.path.join(HERE, 'partitions_flint_lucky_check_73.txt'), 'x', encoding='utf-8') if KEEP else None
NMAX = 3 * 10 ** 16


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def jac(a, n):
    a %= n
    if a == 0: return 0
    return 1 if pow(a, (n - 1) // 2, n) == 1 else -1


def p_mod(N, L):
    return int(flint.fmpz(N).partitions_p()) % L


t0 = time.time()
say('First the control rule on laws known before (37, 43, 53): two controls each, the other class at the same Q')
for L, Q, eps in ((37, 6599, 1), (37, 7541, 1), (37, 9547, -1), (43, 4049, 1), (43, 8941, -1), (53, 467, -1), (53, 3797, -1)):
    tt = 1; got = []
    while len(got) < 2:
        M = L * tt
        if M % 24 == 23 and M % Q and jac(M, Q) == -eps: got.append(M)
        tt += 1
    for Mp in got:
        N = (Q * Q * Mp + 1) // 24; vc = p_mod(N, L)
        pc = (2 * eps * jac(-3, Q) * pow(Q, -2, L) * p_mod((Mp + 1) // 24, L)) % L
        say('   ℓ = %d, Q = %5d, ε = %+d: control p(%d) ≡ %2d, predicted 2ε(−3|Q)Q⁻²·p(%d) = %2d  %s' % (
            L, Q, eps, N, vc, (Mp + 1) // 24, pc, 'as PREDICTED' if vc == pc else 'NOT as predicted'))
text = open(os.path.join(HERE, SRC), encoding='utf-8').read()
say(); say('lists from %s; members computed up to N = %.0e' % (SRC, NMAX))
for blk in re.split(r'\n(?=ℓ = )', text):
    mm = re.match(r'ℓ = (\d+)', blk)
    if not mm: continue
    L = int(mm.group(1)); lines = blk.split('\n'); say(); say('ℓ = %d' % L)
    for i, line in enumerate(lines):
        for kind, eps in (('class ε=+1', 1), ('class ε=−1', -1)):
            if not (line.strip().startswith(kind) and 'expected' in line): continue
            alls = next((l for l in lines[i + 1:i + 3] if l.strip().startswith('all:')), '')
            for Q in [int(v) for v in re.findall(r'\d+', alls)][:2]:
                mem, ctl = [], None; tt = 1
                while len(mem) < 2 or ctl is None:              # M = ℓ·t ≡ 23 (mod 24), Q ∤ M
                    M = L * tt
                    if M % 24 == 23 and M % Q:
                        if jac(M, Q) == eps and len(mem) < 2: mem.append(M)
                        elif jac(M, Q) == -eps and ctl is None: ctl = M
                    tt += 1
                Ns = [(Q * Q * M + 1) // 24 for M in mem]; Nc = (Q * Q * ctl + 1) // 24
                if max(Ns + [Nc]) > NMAX:
                    say('   %-11s Q = %8d: smallest members p(%d), p(%d), control p(%d): past N = %.0e, not computed' % (kind, Q, Ns[0], Ns[1], Nc, NMAX))
                    continue
                t = time.time()
                vals = [p_mod(N, L) for N in Ns]; vc = p_mod(Nc, L)
                pc = (2 * eps * jac(-3, Q) * pow(Q, -2, L) * p_mod((ctl + 1) // 24, L)) % L
                ok = vals == [0, 0] and vc == pc
                say('   %-11s Q = %8d: members p(%d), p(%d) on lane ℓ: %s;  control p(%d): %d, predicted 2ε(−3|Q)Q⁻²·p(%d) = %d  %s  (%.0f s)' % (
                    kind, Q, Ns[0], Ns[1], vals, Nc, vc, (ctl + 1) // 24, pc, 'as PREDICTED' if ok else 'NOT as predicted', time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
