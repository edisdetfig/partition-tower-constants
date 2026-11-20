# A dig, 12 Nov 2026: LINK 1 OF THE PROOF, CHECKED ON NUMBERS — Chen–Du–Hou–Sun, Thm 3.1, with its exact weights.
# Claim (CDHS, read 12 Nov): for a prime ℓ ≥ 5 and b ≥ 1, F_b(z) = Σ_{ℓᵇn ≡ −1 (24)} p((ℓᵇn + 1)/24)·qⁿ is congruent modulo
# ℓᵇ to η(24z)^γ·φ(24z) with φ ∈ M_λ(SL₂(ℤ)), γ = (24β − 1)/ℓᵇ (β ≡ 1/24 mod ℓᵇ, 0 ≤ β < ℓᵇ), and
# λ = (ℓᵇ − ℓᵇ⁻¹)/2 − (γ + 1)/2 for odd b, λ = (ℓᵇ − ℓᵇ⁻¹) − (γ + 1)/2 for even b. The weight on Yang's side is then
# K_b = γ + 2λ − 1 = (ℓ − 1)ℓᵇ⁻¹ − 2 (odd b) or 2(ℓ − 1)ℓᵇ⁻¹ − 2 (even b).
# Method: write n = 24j + γ; then φ(Q) = Σ_j p((ℓᵇ(24j + γ) + 1)/24)·Q^j · ∏(1 − Qⁱ)^(−γ) modulo ℓᵇ. Fit φ to the basis
# E₄^(3(d−1−i))·Δⁱ·E_w of M_λ (d = dim M_λ; unitriangular) from its first d coefficients, and compare every further
# coefficient up to J. PREDICTED: all agree (and if λ is odd or negative, φ ≡ 0). COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/cdhs_check.py ℓ:b … [--J=400] [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly
flint.ctx.threads = 8
import forms6_flint as F6
PAIRS = [tuple(int(x) for x in a.split(':')) for a in sys.argv[1:] if ':' in a]
J = next((int(a[4:]) for a in sys.argv[1:] if a.startswith('--J=')), 400); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'cdhs_check_%s_J%d.txt' % ('_'.join('%d-%d' % p for p in PAIRS), J)), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def sigma_list(L, e):
    s = [0] * L
    for d in range(1, L):
        for n in range(d, L, d): s[n] += d ** e
    return s


def level1_basis(lam, L, M):
    """d series (nmod_poly mod M, length L): E₄^(3(d−1−i))·Δⁱ·E_w, i = 0..d−1, starting at qⁱ."""
    d = lam // 12 + (0 if lam % 12 == 2 else 1); w = lam - 12 * (d - 1)
    s3, s5 = sigma_list(L, 3), sigma_list(L, 5)
    E4 = nmod_poly([1] + [240 * s3[n] % M for n in range(1, L)], M); E6 = nmod_poly([1] + [(-504 * s5[n]) % M for n in range(1, L)], M)
    mul = lambda a, b: a.mul_low(b, L)
    Ew = {0: nmod_poly([1], M), 4: E4, 6: E6, 8: mul(E4, E4), 10: mul(E4, E6), 14: mul(mul(E4, E4), E6)}[w]
    Dq = F6.eta_series(1, L, M).pow_trunc(24, L)                    # Δ/q
    E43 = E4.pow_trunc(3, L); out = []
    for i in range(d):
        g = mul(mul(E43.pow_trunc(d - 1 - i, L), Dq.pow_trunc(i, L)), Ew)
        out.append(mul(g, nmod_poly([0] * i + [1], M)))            # times qⁱ
    return out, d


t0 = time.time()
for L_, b in PAIRS:
    t = time.time(); M = L_ ** b
    beta = pow(24, -1, M); gam = (24 * beta - 1) // M
    lam = ((M - M // L_) // 2 if b % 2 else (M - M // L_)) - (gam + 1) // 2
    K = gam + 2 * lam - 1
    nmax = (M * (24 * (J - 1) + gam) + 1) // 24 + 1
    P = F6.eta_series(1, nmax + 1, M).inverse_series_trunc(nmax + 1)              # p(n) modulo M
    C = nmod_poly([int(P[(M * (24 * j + gam) + 1) // 24]) for j in range(J)], M)
    phi = C.mul_low(F6.eta_series(1, J, M).pow_trunc(gam, J).inverse_series_trunc(J), J)
    head = 'ℓ = %d, b = %d (lane ℓ^%d): γ = %d, λ = %d, K_b = γ + 2λ − 1 = %d (≡ −2 mod (ℓ−1)ℓ^(b−1): %s)' % (
        L_, b, b, gam, lam, K, (K + 2) % ((L_ - 1) * L_ ** (b - 1)) == 0)
    if lam < 0 or lam % 2:
        z = all(int(phi[j]) == 0 for j in range(J))
        say(head + ';  M_λ = 0, so φ ≡ 0 predicted: %s (%.0f s)' % (z, time.time() - t)); continue
    B, d = level1_basis(lam, J, M)
    co = []                                                          # unitriangular solve from φ's first d coefficients
    for i in range(d):
        cur = sum(co[j] * int(B[j][i]) for j in range(i)) % M
        co.append((int(phi[i]) - cur) * pow(int(B[i][i]), -1, M) % M)
    fit = nmod_poly([0], M)
    for i in range(d): fit += B[i] * co[i]
    bad = [j for j in range(d, J) if int(fit[j]) != int(phi[j])]
    say(head + ';  dim M_λ = %d; coefficients beyond the fit that agree: %d of %d%s (%.0f s)' % (
        d, J - d - len(bad), J - d, '' if not bad else ', first disagreement at Q^%d' % bad[0], time.time() - t))
say('seconds: %.0f' % (time.time() - t0))
