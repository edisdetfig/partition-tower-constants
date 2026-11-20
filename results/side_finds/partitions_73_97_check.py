# A dig, 13 Nov 2026: THE TOWERS OF 73 AND 97 ON EXACT PARTITION NUMBERS (README "The towers of 73 and 97"; predicted
# there first). The class S_ℓ and its polynomials come from the forms alone (partitions_flint_family_73.txt, _97.txt):
# for each Q the polynomial of λ = (3|Q)·a_Q on lane ℓ, and P_ℓ, the constants' polynomial on lane ℓ.
#   ALONG ℓ (73 only; 97 would need p at about 10¹⁸):
#     Σ [xⁱ]P_ℓ·H_(1+2i)(n) ≡ 0 (mod ℓ) for n = 0 … 9, and no law with fewer terms.
#   ALONG Q (both): P_Q(T)·c ≡ 0 on 12 chains M, Q²M, …, Q^(2R)M at level 1 (ℓ | M, Q ∤ M; c(M) = p(N), M = 24N − 1), with
#     the step rule of partitions_three_check.py, for Q = 5, 7, 11; and no law of degree R − 1 fits all 12.
#   ONE LEVEL AT A TIME (the user's question about the one-level zoom at ℓ ≡ 1 mod 24): the laws in single-level steps,
#     Σ κᵢ·H_(1+i)(n) ≡ 0 (mod ℓ), i = 0 … d. For each d, the dimension of the space of such laws fitting every n checked.
#     73: d = 1 … 6 (levels 1 to 7, n = 0 … 9); 97: d = 1 … 4 (levels 1 to 5, n = 0 … 12); controls 71, 79: d = 1 … 3.
# Exact partition numbers from FLINT (python-flint 0.9.0, 8 threads). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_73_97_check.py [--keep]
import os, sys, time, re
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
sys.set_int_max_str_digits(0)
import flint
from sympy import jacobi_symbol
flint.ctx.threads = 8
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_73_97_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


t0 = time.time(); PC = {}
leg = lambda a, q: int(jacobi_symbol(a % q, q))


def P(n):
    if n not in PC: PC[n] = int(flint.fmpz(n).partitions_p())
    return PC[n]


def H(L, b, n): return P(L ** b * n + pow(24, -1, L ** b))


def apply_T(x, kap, q3, L):
    return [(x[1] + kap * x[0]) % L] + [(x[j + 1] + q3 * x[j - 1]) % L for j in range(1, len(x) - 1)]


def law_at_0(poly, x, kap, q3, L):
    tot, y = 0, list(x)
    for i, p in enumerate(poly):
        tot = (tot + p * y[0]) % L
        if i < len(poly) - 1: y = apply_T(y, kap, q3, L)
    return tot


def nullity(rows, L):
    """Dimension of {κ : Σ κᵢ·row[i] ≡ 0 (mod L) for every row} (Gaussian elimination over the field of L elements)."""
    A = [[x % L for x in r] for r in rows]; m = len(A[0]); rank = 0
    for c in range(m):
        k = next((i for i in range(rank, len(A)) if A[i][c]), None)
        if k is None: continue
        A[rank], A[k] = A[k], A[rank]; inv = pow(A[rank][c], -1, L); A[rank] = [x * inv % L for x in A[rank]]
        for i in range(len(A)):
            if i != rank and A[i][c]: f = A[i][c]; A[i] = [(x - f * y) % L for x, y in zip(A[i], A[rank])]
        rank += 1
    return m - rank, A[:rank]


def read_family(L):
    txt = open(os.path.join(HERE, 'partitions_flint_family_%d.txt' % L), encoding='utf-8').read()
    polyQ = {int(q): [int(c) for c in co.split(',')] for q, co in re.findall(r'\(3\|Q\)·a_Q on lane ℓ, Q = (\d+): polynomial \[([^\]]*)\]', txt)}
    pl = [int(c) for c in re.search(r"the constants' polynomial \(lowest first\) on lane ℓ: \[([^\]]*)\]", txt).group(1).split(',')]
    dim = int(re.search(r'the class at weight \d+[^:]*: dimension (\d+)', txt).group(1))
    return dim, pl, polyQ


for L in (73, 97):
    t = time.time(); R, pl, polyQ = read_family(L)
    say(); say('ℓ = %d (≡ %d mod 24): %d forms in the class; P_ℓ on lane ℓ (lowest first) = %s; (3|Q)·a_Q polynomials for Q = %s' % (L, L % 24, R, pl, [q for q in sorted(polyQ) if q != L]))
    if L == 73:
        Hs = {b: [H(L, b, n) % L for n in range(10)] for b in range(1, 2 * R + 2)}
        ok = sum(1 for n in range(10) if sum(pl[i] * Hs[1 + 2 * i][n] for i in range(R + 1)) % L == 0)
        fewer = []
        for d in range(1, R):
            k, _ = nullity([[Hs[1 + 2 * i][n] for i in range(d + 1)] for n in range(10)], L)
            fewer.append('%d-term: %d' % (d + 1, k))
        say('   along ℓ: the %d-term law from P_ℓ holds for %d of 10 values of n; laws with fewer terms two levels apart (dimension of the fits): %s  (%.0f s)' % (
            R + 1, ok, ', '.join(fewer), time.time() - t))
    d1 = pow(24, -1, L); rows = []
    for Q in [q for q in sorted(polyQ) if q != L]:
        q3 = pow(Q, -3, L); chains = []; j = 0
        while len(chains) < 12:
            Ns = L * j + d1; M = 24 * Ns - 1; j += 1
            if M % Q == 0: continue
            kap = leg(M, Q) * leg(-3, Q) * pow(Q, -2, L) % L
            x = [P(Q ** (2 * i) * Ns - (Q ** (2 * i) - 1) // 24) % L for i in range(R + 1)]
            chains.append((x, kap))
        ok = sum(1 for x, kap in chains if law_at_0(polyQ[Q], x, kap, q3, L) == 0)
        # laws of degree R − 1 that fit all 12 chains: the linear conditions on the coefficients (monic), counted by elimination
        cond = []
        for x, kap in chains:
            y = list(x); vals = []
            for i in range(R):
                vals.append(y[0]);
                if i < R - 1: y = apply_T(y, kap, q3, L)
            cond.append(vals)
        k_lower, _ = nullity(cond, L)
        rows.append('Q = %d: %d of 12 (laws of degree %d fitting all 12: dimension %d)' % (Q, ok, R - 1, k_lower))
    say('   along Q: the %d-step law from the forms holds at %s  (%.0f s)' % (R, ';  '.join(rows), time.time() - t))

say(); say('One level at a time: the space of laws Σ κᵢ·H_(1+i)(n) ≡ 0 (mod ℓ), i = 0 … d, fitting every n checked')
for L, dmax, nn in ((73, 6, 10), (97, 4, 13), (71, 3, 12), (79, 3, 12)):
    t = time.time()
    Hs = {b: [H(L, b, n) % L for n in range(nn)] for b in range(1, dmax + 2)}
    res = []
    for d in range(1, dmax + 1):
        k, _ = nullity([[Hs[1 + i][n] for i in range(d + 1)] for n in range(nn)], L)
        res.append('d = %d: %d' % (d, k))
    say('   ℓ = %d (≡ %2d mod 24), levels 1 to %d, n = 0 … %d: %s  (%.0f s)' % (L, L % 24, dmax + 1, nn - 1, ';  '.join(res), time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
