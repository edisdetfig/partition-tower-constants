# A dig, 13 Nov 2026: THE TWINS TAKEN APART, EACH CONSTANT ALONE, DIGIT BY DIGIT (README "The twins taken apart"; predicted
# there first). At 37, 41, 43, 53 the class is two forms f₁, f₂ with ℓ-adic constants α₁, α₂ and step numbers at 5,
# λ₁(5), λ₂(5) (partitions_flint_family_steps_ℓ.txt, 8 digits; paired with each other through the forms' eigenvectors on
# lane ℓ). One step at 5 cancels f_j:
#     G_i(M) = c(25M) + (κ − λ_j(5))·c(M) = (λ_i(5) − λ_j(5))·c_i(M),   κ = (M|5)·(−3|5)·5⁻²,   c(M) = p((M + 1)/24),
# so G_i(ℓ²M)/G_i(M) should be α_i. Read off exact partition numbers modulo ℓ⁸ at levels b → b + 2 (M = ℓ^b·(24n + k),
# k ≡ −ℓ^(−b) mod 24), b = 1, 3, 5, 7, while p(N) has N ≤ 1.2·10¹⁶; and, for contrast, the whole tower's c(ℓ²M)/c(M).
# FLINT's exact partition numbers, at idle priority. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_twins_split.py [--keep]
import os, sys, re, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat
import forms6_flint as F
flint.ctx.threads = 6
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_twins_split.txt'), 'x', encoding='utf-8') if KEEP else None
D = 8; NMAX = 12 * 10 ** 15


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


PC = {}


def pn(N):
    if N not in PC: PC[N] = int(flint.fmpz(N).partitions_p())
    return PC[N]


def agree(x, digs, L):
    """How many lowest base-L digits of x agree with digs."""
    n = 0
    while n < len(digs) and (x // L ** n) % L == digs[n]: n += 1
    return n


t0 = time.time()
for L in (37, 41, 43, 53):
    txt = open(os.path.join(HERE, 'partitions_flint_family_steps_%d.txt' % L), encoding='utf-8').read()
    alphas = [[int(x) for x in g.split(',')] for g in re.findall(r'a constant on lane ℓ, digits lowest first: \[([^\]]*)\]', txt)]
    lams = [[int(x) for x in g.split(',')] for g in re.findall(r'a root on lane ℓ, digits lowest first: \[([^\]]*)\]', txt)]
    # pair them on lane ℓ: the joint eigenvectors of T_5 and T_ℓ on the class
    k = L - 3; n0 = k + 1
    Kz, d = F.class_basis_exact(k, leg(2, L), leg(3, L))
    Xc, Rc = F.normalized_basis(nmod_mat(n0, d, [int(Kz[i, j]) % L for i in range(n0) for j in range(d)], L), d, L)
    H = F.hecke_matrices(k, L, [('T', 5), ('T', L)])
    T5 = F.restrict(H[('T', 5)], Xc, Rc, d); TL = F.restrict(H[('T', L)], Xc, Rc, d)
    pairs = []
    for lam in lams:
        r = lam[0]; Kv, nul = (T5 * leg(3, 5) - F.identity(d, L) * r).nullspace(); assert nul == 1
        v = [int(Kv[i, 0]) for i in range(d)]; i0 = next(i for i in range(d) if v[i])
        e = sum(int(TL[i0, j]) * v[j] for j in range(d)) * pow(v[i0], -1, L) % L
        al = next(a for a in alphas if a[0] == leg(3, L) * e % L)
        pairs.append((lam, al))
    say(); say('ℓ = %d: twin 1: λ(5) %s…, constant %s…;  twin 2: λ(5) %s…, constant %s…  (paired on lane ℓ by the eigenvectors)' % (
        L, pairs[0][0][:4], pairs[0][1][:5], pairs[1][0][:4], pairs[1][1][:5]))
    m = L ** D; num = lambda dg: sum(x * L ** i for i, x in enumerate(dg)) % m
    lamv = [num(p[0]) for p in pairs]
    kk = None
    for b in (1, 3, 5, 7):
        kk = (-pow(L ** b, -1, 24)) % 24
        t = time.time(); res = {0: [], 1: []}; whole = []; skipped = 0; n = 0; used = 0
        while used < 6:
            M = L ** b * (24 * n + kk); n += 1
            if M % 5 == 0: continue
            if (25 * L * L * M + 1) // 24 > NMAX: break
            used += 1
            kap = leg(M, 5) * leg(-3, 5) * pow(5, -2, m) % m
            cM, c25M = pn((M + 1) // 24) % m, pn((25 * M + 1) // 24) % m
            cL, c25L = pn((L * L * M + 1) // 24) % m, pn((25 * L * L * M + 1) // 24) % m
            for i in (0, 1):
                lj = lamv[1 - i]
                g0 = (c25M + (kap - lj) * cM) % m; g1 = (c25L + (kap - lj) * cL) % m
                if g0 % L == 0: skipped += 1; continue
                ratio = g1 * pow(g0, -1, m) % m
                res[i].append((agree(ratio, pairs[i][1], L), agree(ratio, pairs[1 - i][1], L)))
            if cM % L: whole.append(max(agree(cL * pow(cM, -1, m) % m, pairs[i][1], L) for i in (0, 1)))
        if not used: say('   levels %d → %d: past N = %.1e, not computed' % (b, b + 2, NMAX)); continue
        say('   levels %d → %d (%d values of M): digits agreeing with its own constant / with the other twin\'s:  twin 1: %s;  twin 2: %s;  '
            'the whole tower c(ℓ²M)/c(M), best of the two constants: %s%s  (%.0f s)' % (
                b, b + 2, used, ' '.join('%d/%d' % x for x in res[0]), ' '.join('%d/%d' % x for x in res[1]), whole,
                ('; %d skipped (G(M) ≡ 0 on lane ℓ)' % skipped) if skipped else '', time.time() - t))
say(); say('the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
