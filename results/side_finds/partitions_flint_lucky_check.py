# A dig, 12 Nov 2026: NEW LUCKY PRIMES, CHECKED ON EXACT PARTITION NUMBERS. partitions_flint_lucky_*.txt lists, from the
# forms alone, every prime Q < 10⁶ where all forms of ℓ's class take a law's number on lane ℓ (STATEMENT.md, E1). Each such Q
# is a congruence (AAT's Lemma 5.1): for a class law with ε, p(N) ≡ 0 (mod ℓ) whenever 24N − 1 = Q²·M, ℓ | M, Q ∤ M,
# (M|Q) = ε; for Ono's kind, p(N) ≡ 0 whenever 24N − 1 = ℓ·Q³·m, Q ∤ m.
# PREDICTED before the run: for the two smallest new Q of each kind for 37 to 59 (2999, 4049, 6599, 7541, 8941, 9547 were
# checked before), the two smallest members are multiples of ℓ; a control (the other class, or for Ono's kind a member of
# Q²·M's shape) is not, as a rule (it is 0 with chance 1/ℓ). FLINT's exact partition numbers. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_lucky_check.py [--keep]
import os, sys, re, time
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
flint.ctx.threads = 4
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_flint_lucky_check.txt'), 'x', encoding='utf-8') if KEEP else None
SRC = 'partitions_flint_lucky_37_43_41_47_53_59_61_13_17_19_23_29_31.txt'
DONE = {2999, 4049, 6599, 7541, 8941, 9547}


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def jac(a, n):                                                 # Jacobi symbol (a|n), n odd prime here
    a %= n
    if a == 0: return 0
    return 1 if pow(a, (n - 1) // 2, n) == 1 else -1


def p_mod(N, L):
    return int(flint.fmpz(N).partitions_p()) % L


t0 = time.time(); text = open(os.path.join(HERE, SRC), encoding='utf-8').read()
say('lists from %s' % SRC)
for blk in re.split(r'\n(?=ℓ = )', text):
    mm = re.match(r'ℓ = (\d+)', blk)
    if not mm or not 37 <= int(mm.group(1)) <= 59: continue
    L = int(mm.group(1)); lines = blk.split('\n'); say(); say('ℓ = %d' % L)
    for i, line in enumerate(lines):
        for kind, eps in (('class ε=+1', 1), ('class ε=−1', -1), ('Ono', 0)):
            if not (line.strip().startswith(kind) and 'expected' in line): continue
            alls = next((l for l in lines[i + 1:i + 3] if l.strip().startswith('all:')), '')
            Qs = [q for q in (int(v) for v in re.findall(r'\d+', alls)) if q not in DONE][:2]
            for Q in Qs:
                t = time.time(); out = []
                if eps:                                        # 24N − 1 = Q²·M, M = ℓ·t ≡ 23 (mod 24), (M|Q) = ±ε
                    mem, ctl = [], None; tt = 1
                    while len(mem) < 2 or ctl is None:
                        M = L * tt
                        if M % 24 == 23 and M % Q:
                            if jac(M, Q) == eps and len(mem) < 2: mem.append(M)
                            elif jac(M, Q) == -eps and ctl is None: ctl = M
                        tt += 1
                    Ns = [(Q * Q * M + 1) // 24 for M in mem]; Nc = (Q * Q * ctl + 1) // 24
                else:                                          # 24N − 1 = ℓ·Q³·m; control 24N − 1 = Q²·M, ℓ | M, Q ∤ M
                    mem = []; mm_ = 1
                    while len(mem) < 2:
                        if (L * Q ** 3 * mm_) % 24 == 23 and mm_ % Q: mem.append(mm_)
                        mm_ += 1
                    Ns = [(L * Q ** 3 * x + 1) // 24 for x in mem]; tt = 1
                    while True:
                        M = L * tt
                        if M % 24 == 23 and M % Q: break
                        tt += 1
                    Nc = (Q * Q * M + 1) // 24
                vals = [p_mod(N, L) for N in Ns]; vc = p_mod(Nc, L)
                say('   %-11s Q = %6d: members p(%d), p(%d) on lane ℓ: %s;  control p(%d): %d  %s  (%.0f s)' % (
                    kind, Q, Ns[0], Ns[1], vals, Nc, vc, 'as PREDICTED' if vals == [0, 0] else 'NOT as predicted', time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
