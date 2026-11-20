# A dig, 13 Nov 2026: THE LAWS FOR 73 AND 61 WRITTEN AS CONGRUENCES (README "Lucky primes for 73"). DERIVED from the lists
# of partitions_flint_lucky_dual_73_61_max50000000.txt, by arithmetic alone:
#   class law (Q, ε): p(N) ≡ 0 (mod ℓ) whenever 24N − 1 = Q²·M, ℓ | M, Q ∤ M, (M|Q) = ε. For M ≡ M₀ (mod 24ℓQ) this is the
#     progression N = ℓ·Q³·n + β with β = (Q²·M₀ + 1)/24; there are (Q − 1)/2 such β modulo ℓ·Q³.
#   Ono's kind (Q): p(N) ≡ 0 whenever 24N − 1 = ℓ·Q³·m, Q ∤ m: N = ℓ·Q⁴·n + β, β = (ℓ·Q³·m₀ + 1)/24, Q − 1 values of β.
# Printed: the smallest β of each law, and a check that 24β − 1 has the stated shape. DERIVED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_lucky_progressions.py [--keep]
import os, sys, re
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRC = 'partitions_flint_lucky_dual_73_61_max50000000.txt'
fh = open(os.path.join(HERE, 'partitions_lucky_progressions.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def jac(a, n):
    a %= n
    if a == 0: return 0
    return 1 if pow(a, (n - 1) // 2, n) == 1 else -1


text = open(os.path.join(HERE, SRC), encoding='utf-8').read()
say('the laws of %s, as congruences (smallest β of each; DERIVED)' % SRC)
for blk in re.split(r'\n(?=ℓ = )', text):
    mm = re.match(r'ℓ = (\d+)', blk)
    if not mm: continue
    L = int(mm.group(1)); lines = blk.split('\n'); say(); say('ℓ = %d' % L)
    for i, line in enumerate(lines):
        for kind, eps in (('class ε=+1', 1), ('class ε=−1', -1), ('Ono', 0)):
            if not (line.strip().startswith(kind) and 'expected' in line): continue
            alls = next((l for l in lines[i + 1:i + 3] if l.strip().startswith('all:')), '')
            for Q in [int(v) for v in re.findall(r'\d+', alls)]:
                if eps:
                    t = 1
                    while True:
                        M = L * t
                        if M % 24 == 23 and M % Q and jac(M, Q) == eps: break
                        t += 1
                    beta = (Q * Q * M + 1) // 24; mod = L * Q ** 3
                    ok = (24 * beta - 1) % (Q * Q) == 0 and ((24 * beta - 1) // (Q * Q)) % L == 0 and jac((24 * beta - 1) // (Q * Q), Q) == eps
                    say('   %-11s Q = %8d: p(%d·n + %d) ≡ 0 (mod %d)   [modulus ℓ·Q³; %d values of β in all; 24β − 1 = Q²·%d%s]' % (
                        kind, Q, mod, beta, L, (Q - 1) // 2, M, '' if ok else ' — SHAPE WRONG'))
                else:
                    m0 = 1
                    while (L * Q ** 3 * m0) % 24 != 23 or m0 % Q == 0: m0 += 1
                    beta = (L * Q ** 3 * m0 + 1) // 24; mod = L * Q ** 4
                    ok = (24 * beta - 1) % (L * Q ** 3) == 0 and ((24 * beta - 1) // (L * Q ** 3)) % Q != 0
                    say('   %-11s Q = %8d: p(%d·n + %d) ≡ 0 (mod %d)   [modulus ℓ·Q⁴; %d values of β in all; 24β − 1 = ℓ·Q³·%d%s]' % (
                        kind, Q, mod, beta, L, Q - 1, m0, '' if ok else ' — SHAPE WRONG'))
