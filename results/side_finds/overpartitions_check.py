# A dig, 12 Nov 2026: OVERPARTITIONS AGAINST THE FORMS. p̄(n) counts overpartitions: Σ p̄(n)qⁿ = ∏(1 + qⁿ)/(1 − qⁿ) = 1/θ₄,
# θ₄ = Σ(−1)ⁿq^(n²) = η(z)²/η(2z), weight −1/2. The tower is H_b(n) = p̄(ℓᵇn); two levels are the step U(ℓ²).
# The data-only laws (group_members_….txt) are, in y = the two-level step:
#   - y − 1 at 5, 7 and 11;
#   - at 13, y² + 5y + 7 = (y − 1)(y − 7) on lane 13.
# The newform of weight 10 = 13 − 3 on Γ0(2) has a₁₃ ≡ 7 (mod 13), and the weight-8 form h₈ (opposite sign at 2) shows no
# trace at 11.
# PREDICTED (derived from those four primes, as the recipe at weight −1/2 → −2 reads): the constants are the Eisenstein 1
# and ψ(ℓ)·a_ℓ(f) for the newforms f of weight ℓ − 3 on Γ0(2) in one Atkin–Lehner class at 2. Which class (w₂ = −1 always,
# or carried by (2|ℓ)) and which twist ψ (1, (−1|·), (2|·), (−2|·)) the data must decide. Tested here at 13 … 31 on lane
# ℓ: every candidate law, from every level, every n. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/overpartitions_check.py N ℓ … [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly, nmod_mat
flint.ctx.threads = 8
import forms2_flint as F2
import forms6_flint as F
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; N, PS = nums[0], nums[1:]; KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'overpartitions_check_N%d_%s.txt' % (N, '_'.join(map(str, PS)))), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def euler(M, step, n):
    P = nmod_poly([], M); terms = [(0, 1)]; j = 1
    while step * j * (3 * j - 1) // 2 < n:
        terms.append((step * j * (3 * j - 1) // 2, (-1) ** j))
        if step * j * (3 * j + 1) // 2 < n: terms.append((step * j * (3 * j + 1) // 2, (-1) ** j))
        j += 1
    for e, c in sorted(terms, reverse=True): P[e] = c % M
    return P


def class_roots(L, k, s2):
    """T_ℓ's eigenvalues on lane ℓ (as the characteristic polynomial mod ℓ) on the newforms of weight k on Γ0(2) with
    U₂ = s2·2^((k−2)/2) (exact class, forms2_flint)."""
    K, d = F2.class_basis_exact2(k, s2)
    if d == 0: return [1], 0
    n = k // 4 + 1
    H = F2.hecke_matrices2(k, L, [('T', L)])[('T', L)]
    X = nmod_mat(n, d, [int(K[i, j]) % L for i in range(n) for j in range(d)], L)
    if len(F.pivot_rows(X, L)) < d: return None, d
    Xn, Rn = F.normalized_basis(X, d, L)
    T = F.restrict(H, Xn, Rn, d)
    return [int(c) for c in T.charpoly().coeffs()], d


def val(x, L, K):
    x %= L ** K
    if x == 0: return K
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


t0 = time.time()
for L in PS:
    t = time.time(); K = 1
    while L ** (K + 1) < 2 ** 30: K += 1
    M = L ** K; NN = N
    E1, E2 = euler(M, 1, NN), euler(M, 2, NN)
    Pb = E2.mul_low(E1.pow_trunc(2, NN).inverse_series_trunc(NN), NN); del E1, E2
    lv = {}
    for b in range(0, 40):
        if L ** b >= N: break
        lv[b] = [int(Pb[L ** b * n]) % L for n in range(12) if L ** b * n < N]
    del Pb
    top = max(lv); k = L - 3
    say(); say('ℓ = %d (%.0f s): tower levels 0..%d; (2|ℓ) = %+d, (−1|ℓ) = %+d' % (L, time.time() - t, top, leg(2, L), leg(-1, L)))
    roots = {}
    for s2 in (1, -1):
        cp, d = class_roots(L, k, s2)
        roots[s2] = (cp, d)
        say('   newforms of weight %d on Γ0(2) with U₂ = %+d·2^%d (w₂ = %+d): %d; T_ℓ charpoly on lane ℓ %s' % (
            k, s2, (k - 2) // 2, -s2, d, cp))
    hyps = {'w₂ = −1': 1, 'w₂ = +1': -1, 'w₂ = (2|ℓ)': -leg(2, L), 'w₂ = −(2|ℓ)': leg(2, L)}
    tw = {'ψ = 1': 1, 'ψ = (−1|ℓ)': leg(-1, L), 'ψ = (2|ℓ)': leg(2, L), 'ψ = (−2|ℓ)': leg(-2, L)}
    for hn, s2 in hyps.items():
        cp, d = roots[s2]
        if cp is None: say('   %s: basis degenerate on lane ℓ' % hn); continue
        for tn, ps in tw.items():
            # law in y: (y − 1)·∏(y − ψ·a): if χ(x) = ∏(x − a) then ∏(y − ψa) = ψ^d·χ(ψ·y)
            chi = nmod_poly(cp, L); co = [int(c) for c in chi.coeffs()]
            tw_co = [(co[i] * pow(ps, i, L) * pow(ps, len(co) - 1, L)) % L for i in range(len(co))]   # χ(ψy)·ψ^d
            Q = (nmod_poly([L - 1, 1], L) * nmod_poly(tw_co, L)).coeffs(); Q = [int(c) for c in Q]; D = len(Q) - 1
            ok = True; tested = 0; worst = None
            for b in sorted(lv):
                if b < 1 or b + 2 * D > top: continue
                nn = min(len(lv[b + 2 * i]) for i in range(D + 1))
                for n in range(nn):
                    s = sum(Q[i] * lv[b + 2 * i][n] for i in range(D + 1)) % L; tested += 1
                    if s: ok = False; worst = (b, n); break
                if not ok: break
            say('   %s, %s: law in y %s (lowest first): %s on %d checks%s' % (
                hn, tn, Q, 'HOLDS' if ok and tested else ('no checks' if not tested else 'fails'), tested,
                '' if ok else ' (first failure at level %d, n = %d)' % worst))
say(); say('seconds: %.0f' % (time.time() - t0))
