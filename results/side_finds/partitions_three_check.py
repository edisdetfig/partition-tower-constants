# A dig, 12 Nov 2026: THREE FORMS FROM 61 ON, THE PARTITION NUMBERS' SIDE. partitions_three_forms.txt found the tower's class
# one cubic Galois orbit for 61 (two forms for 53, the control) and gives, for each Q, the polynomial P_Q of the class's
# numbers λ = (3|Q)·a_Q on lane ℓ. PREDICTED before the run, by the same reasoning as for the twins (partitions_split.py):
#   ALONG ℓ: Σ p_i·H_(1+2i)(n) ≡ 0 (mod ℓ) with the coefficients p_i of P_ℓ (61: H₇ with H₅, H₃, H₁; 53: H₅ with H₃, H₁),
#            for n = 0 to 9; and no law with fewer terms fits (no T, N with H₅ ≡ T·H₃ − N·H₁ for 61).
#   ALONG Q: P_Q(T)·c ≡ 0 at every chain M, Q²M, Q⁴M, Q⁶M (ℓ | M, Q ∤ M; c(M) = p(N), M = 24N − 1), where T is the step rule
#            T·c(Q^(2j)M) = c(Q^(2j+2)M) + [j = 0]·ε(−3|Q)Q⁻²·c(M) + [j ≥ 1]·Q⁻³·c(Q^(2j−2)M), ε = (M|Q); 12 chains for each
#            Q = 5, 7, 11, 13; and no law of lower degree fits all of them.
# Exact partition numbers from FLINT (python-flint 0.9.0, 8 threads). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_three_check.py [--keep]
import os, sys, time, re
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
import flint
from sympy import jacobi_symbol
flint.ctx.threads = 8
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_three_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); PC = {}
leg = lambda a, q: int(jacobi_symbol(a % q, q))


def P(n):
    if n not in PC: PC[n] = int(flint.fmpz(n).partitions_p())
    return PC[n]


POLY = {}
for L, Q, co in re.findall(r'ℓ = (\d+), Q = (\d+): \[([^\]]*)\]', open(os.path.join(HERE, 'partitions_three_forms.txt'), encoding='utf-8').read()):
    POLY[int(L), int(Q)] = [int(c) for c in co.split(',')]


def apply_T(x, kap, q3, L):                                     # the step rule on a chain x_j = c(Q^(2j)·M)
    return [(x[1] + kap * x[0]) % L] + [(x[j + 1] + q3 * x[j - 1]) % L for j in range(1, len(x) - 1)]


def law_at_0(poly, x, kap, q3, L):                              # (P(T)·c)(M) for P = Σ poly[i]·T^i
    tot, y = 0, list(x)
    for i, p in enumerate(poly):
        tot = (tot + p * y[0]) % L
        if i < len(poly) - 1: y = apply_T(y, kap, q3, L)
    return tot


for L in (53, 61):
    t = time.time(); d = len(POLY[L, L]) - 1; say(); say('ℓ = %d: %d forms in the class (partitions_three_forms.txt)' % (L, d))
    # ---- along ℓ ----
    H = {b: [P(L ** b * n + pow(24, -1, L ** b)) % L for n in range(10)] for b in range(1, 2 * d + 2, 2)}
    pl = POLY[L, L]; ok = sum(1 for n in range(10) if sum(pl[i] * H[1 + 2 * i][n] for i in range(d + 1)) % L == 0)
    say('   along ℓ: the %d-term law from P_ℓ = %s holds for %d of 10 values of n  (%.0f s)' % (d + 1, pl, ok, time.time() - t))
    if d == 3:
        two = [(T, N) for T in range(L) for N in range(L) if all((H[5][n] - T * H[3][n] + N * H[1][n]) % L == 0 for n in range(10))]
        one = [A for A in range(L) if all((H[3][n] - A * H[1][n]) % L == 0 for n in range(10))]
        say('      fewer terms: (T, N) with H₅ ≡ T·H₃ − N·H₁: %s;  one constant with H₃ ≡ A·H₁: %s' % (two or 'none', one or 'none'))
    # ---- along Q ----
    d1 = pow(24, -1, L); rows = []; lower_total = 0
    for Q in (5, 7, 11, 13):
        q3 = pow(Q, -3, L); chains = []; j = 0
        while len(chains) < 12:
            Ns = L * j + d1; M = 24 * Ns - 1; j += 1
            if M % Q == 0: continue
            kap = leg(M, Q) * leg(-3, Q) * pow(Q, -2, L) % L
            x = [P(Q ** (2 * i) * Ns - (Q ** (2 * i) - 1) // 24) % L for i in range(d + 1)]
            chains.append((x, kap))
        ok = sum(1 for x, kap in chains if law_at_0(POLY[L, Q], x, kap, q3, L) == 0)
        lower = sum(1 for a in range(L) for b in range(L) if d == 3 and all(law_at_0([b, a, 1], x, kap, q3, L) == 0 for x, kap in chains))
        lower_total += lower
        rows.append('Q = %d: %d of 12%s' % (Q, ok, '' if d < 3 else ' (laws of degree 2 that fit all 12: %d)' % lower))
    say('   along Q: the %d-step law holds at %s  (%.0f s)' % (d, ';  '.join(rows), time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
