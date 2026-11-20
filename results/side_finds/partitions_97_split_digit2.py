# A dig, 13 Nov 2026: THE SECOND DIGIT OF 97'S CONSTANT ON PARTITION NUMBERS (README "97's tower is two towers", "Next: the
# second digit"; predicted there first). partitions_97_split.py cancelled the three conjugate forms of 97's class with the
# step at Q′ = 5 (κ from the cubic factor h of λ(5)'s polynomial) and saw G(97²M) ≡ 28·G(M) (mod 97). Here κ is taken to
# two and three digits, from the weight −2 member's λ(5) (partitions_flint_family_steps_97.txt: its root λ_g(5) and the
# rest h on lane 97⁴), and the multiplier G(97²M)/G(M) is read off exact partition numbers modulo 97² and 97³:
#   levels 3 → 5 (M = 97³·(24n + 23)): PREDICTED 2550 = 28 + 26·97 modulo 97², for n = 0 … 9 (5 ∤ M);
#   levels 1 → 3 (M = 97·(24n + 23)): PREDICTED 28 modulo 97, but not 2550 modulo 97² for most n;
#   modulo 97³ at levels 3 → 5: measured, not predicted (the constant's third digit is 81).
# FLINT's exact partition numbers, at idle priority. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_97_split_digit2.py [--keep]
import os, sys, re, time
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
flint.ctx.threads = 6
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_97_split_digit2.txt'), 'x', encoding='utf-8') if KEEP else None
L = 97; QP = 5


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


txt = open(os.path.join(HERE, 'partitions_flint_family_steps_97.txt'), encoding='utf-8').read()
lam_d = [int(x) for x in re.search(r'a root on lane ℓ, digits lowest first: \[([^\]]*)\]', txt).group(1).split(',')]
rest = [int(x) for x in re.search(r'the rest \(lowest first\) on lane ℓ\^4: \[([^\]]*)\]', txt).group(1).split(',')]
al_d = [int(x) for x in re.search(r'a constant on lane ℓ, digits lowest first: \[([^\]]*)\]', txt).group(1).split(',')]
num = lambda dg, e: sum(x * L ** i for i, x in enumerate(dg[:e]))
say('from partitions_flint_family_steps_97.txt: λ_g(5) digits %s, the constant α digits %s, h on lane 97⁴ %s' % (lam_d, al_d, rest))
say('   h on lane 97: %s (partitions_97_split.txt had [22, 13, 13, 1])' % [x % L for x in rest])


def kappa(e5, mod):
    kap = e5 * leg(-3, QP) * pow(QP, -2, mod) % mod; q3 = pow(QP, -3, mod); h = [x % mod for x in rest]
    out = []
    for j in range(4):
        y = [1 if jj == j else 0 for jj in range(4)]; tot = 0
        for i, p in enumerate(h):
            tot = (tot + p * y[0]) % mod
            if i < len(h) - 1: y = [(y[1] + kap * y[0]) % mod] + [(y[t + 1] + q3 * y[t - 1]) % mod for t in range(1, len(y) - 1)]
        out.append(tot)
    return out


PC = {}


def c(M, mod):
    N = (M + 1) // 24
    if N not in PC: PC[N] = int(flint.fmpz(N).partitions_p())
    return PC[N] % mod


def G(M, mod):
    K = kappa(leg(M, QP), mod)
    return sum(K[j] * c(QP ** (2 * j) * M, mod) for j in range(4)) % mod


def digits(x, e): return [(x // L ** i) % L for i in range(e)]


t0 = time.time()
for e in (2, 3):
    say('   κ on lane 97^%d: (M|5) = +1: %s;  (M|5) = −1: %s' % (e, kappa(1, L ** e), kappa(-1, L ** e)))
for (lo, hi), ns in (((3, 5), range(14)), ((1, 3), range(14))):
    say(); say('Levels %d → %d: M = 97^%d·(24n + 23); the multiplier G(97²M)/G(M) read off exact partition numbers' % (lo, hi, lo))
    t = time.time(); rows = []; hit2 = hit1 = cnt = 0; hit3 = 0
    for n in ns:
        M = L ** lo * (24 * n + 23)
        if M % QP == 0: continue
        if cnt == 10: break
        g0, g1 = G(M, L ** 3), G(L * L * M, L ** 3); cnt += 1
        if g0 % L == 0: rows.append('n = %d: G(M) ≡ 0 (mod 97)' % n); continue
        r = g1 * pow(g0, -1, L ** 3) % L ** 3; dg = digits(r, 3)
        hit1 += dg[0] == al_d[0]; hit2 += dg[:2] == al_d[:2]; hit3 += dg == al_d[:3]
        rows.append('n = %d: %s' % (n, dg))
    say('   multiplier digits (lowest first), mod 97³: %s' % ';  '.join(rows))
    say('   first digit 28: %d of %d;  first two digits 28, 26 (2550 mod 97²): %d of %d;  three digits 28, 26, 81: %d of %d  (%.0f s)' % (
        hit1, cnt, hit2, cnt, hit3, cnt, time.time() - t))
say(); say('the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
