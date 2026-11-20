# A dig, 12 Nov 2026: THE "SECOND KIND" OF EXCEPTION — A CLASS FORM WHOSE NUMBER AT ℓ IS ITSELF DIVISIBLE BY ℓ.
# The reasoning (derived, mod-ℓ theory of theta cycles, Jochnowitz; Serre weights): if f has weight k = ℓ − 3 on Γ0(6) and
# a_ℓ(f) ≡ 0 (mod ℓ), then g = θ⁴f is a mod-ℓ eigenform of weight ℓ + 3 − k = 6 on Γ0(6), again with a_ℓ(g) ≡ 0, and
# U₂(g) ≡ 2⁴·a₂(f), U₃(g) ≡ 3⁴·a₃(f). With a₂(f) = −w₂·2^((k−2)/2) and 2^((ℓ−1)/2) ≡ (2|ℓ): U₂(g) ≡ −w₂(f)·(2|ℓ)·4, and
# U₃(g) ≡ −w₃(f)·(3|ℓ)·9. By Deligne–Serre g lifts to an eigenform of weight 6 on Γ0(6). Those are: Eisenstein series
# (ordinary, excluded); the level-3 newform η(z)⁶η(3z)⁶ stabilized at 2 (U₂ a root of x² + 6x + 32, U₃ = +9); and the one
# newform h₆ of level 6, weight 6, with U₂ = +4, U₃ = −9 (class_basis_exact). So:
#   - the partition class (w₂, w₃) = (−(2|ℓ), −(3|ℓ)) needs U₂ ≡ 4, U₃ ≡ 9: no lift for any ℓ ≥ 5 (h₆ has U₃ = −9; the
#     level-3 form needs ℓ | 72). PREDICTED: every partition-class form is ordinary at every ℓ ≥ 7;
#   - the class (−(2|ℓ), +(3|ℓ)) needs U₂ ≡ 4, U₃ ≡ −9: exactly h₆. PREDICTED: it has a non-ordinary form exactly when
#     ℓ | a_ℓ(h₆);
#   - the classes (+(2|ℓ), ±(3|ℓ)) need U₂ ≡ −4: no lift (ℓ | 24 for the level-3 form). PREDICTED: always ordinary.
# Part A: h₆'s q-expansion (FLINT, modulo a 61-bit prime; |a_p| ≤ 2p^(5/2) is recovered exactly for p ≤ 10⁷) and every
# prime p with p | a_p(h₆). Part B: for each prime ℓ from 13 to LMAX and each of the four sign classes at weight ℓ − 3 (exact
# over ℤ, class_basis_exact), whether T_ℓ has the eigenvalue 0 on lane ℓ. COUNTED.
# A first run (n < 20,000) found h₆ non-ordinary at exactly one prime, 859. PREDICTED before Part B is run there: at
# ℓ = 859 the class (−(2|ℓ), +(3|ℓ)) has a form with T_ℓ-eigenvalue 0 on lane ℓ, and the other three classes (the partition
# class among them) do not; the neighbours 853, 857, 863 (controls) have none in any class.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/second_kind.py NA LMAX [--at=ℓ,ℓ,…] [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_poly, nmod_mat, fmpz_mat
import forms6_flint as F
flint.ctx.threads = 8
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; NA, LMAX = nums[0], nums[1]; KEEP = '--keep' in sys.argv
AT = [int(x) for x in next((a[5:] for a in sys.argv[1:] if a.startswith('--at=')), '').split(',') if x]
fh = open(os.path.join(HERE, 'second_kind_N%d_L%d%s.txt' % (NA, LMAX, ('_at' + '_'.join(map(str, AT))) if AT else '')), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if s[i]]


t0 = time.time()
# ---- Part A: h₆ = the newform of weight 6 on Γ0(6) (U₂ = +4, U₃ = −9), to NA terms -----------------------------------
K, d = F.class_basis_exact(6, 1, -1); assert d == 1
first = [int(K[a, 0]) for a in range(7)]; g0 = first[1]; first = [c // g0 for c in first]
say('h₆ (weight 6, level 6, U₂ = +4, U₃ = −9, so w₂ = −1, w₃ = +1): first coefficients %s' % first)
P = 2 ** 61 - 1; L = NA + 1; k = 6
E = {dd: F.eta_series(dd, L, P) for dd in (1, 2, 3, 6)}
inv = lambda f: f.inverse_series_trunc(L)
mul = lambda a, b: a.mul_low(b, L)
u = mul(mul(E[1].pow_trunc(5, L), E[3]), inv(mul(E[2], E[6].pow_trunc(5, L))))          # t·q = u: η⁵(z)η(3z)/(η(2z)η⁵(6z))·q
phi = mul(mul(E[1], E[6].pow_trunc(6, L)), inv(mul(E[2].pow_trunc(2, L), E[3].pow_trunc(3, L)))).pow_trunc(k, L)  # G_k/q^k
H = []                                                                                      # H_m = q^m·φ·u^(k−m)
upow = [nmod_poly([1], P)]
for j in range(k): upow.append(mul(upow[-1], u))
for m in range(k + 1):
    H.append(mul(mul(phi, upow[k - m]), nmod_poly([0] * m + [1], P)))
c = [0] * (k + 1); target = [x % P for x in first]                                         # solve Σ c_m H_m = h₆ (unitriangular)
for m in range(k + 1):
    cur = sum(c[j] * int(H[j][m]) for j in range(m)) % P
    c[m] = (target[m] - cur) * pow(int(H[m][m]), -1, P) % P
h = nmod_poly([], P)
for m in range(k + 1):
    if c[m]: h += H[m] * c[m]
coef = lambda n: (lambda v: v - P if v > P // 2 else v)(int(h[n]))
assert [coef(n) for n in range(7)] == first
ps = primes_upto(NA - 1); nonord = [p for p in ps if p >= 5 and coef(p) % p == 0]
say('   a_p for primes p < %d computed (%.0f s); checks: a₅ = %d, a₇ = %d, a₁₁ = %d, a₁₃ = %d' % (
    NA, time.time() - t0, coef(5), coef(7), coef(11), coef(13)))
say('   primes 5 ≤ p < %d with p | a_p(h₆) (h₆ non-ordinary at p): %s' % (NA, nonord))
# ---- Part B: the four sign classes at weight ℓ − 3, for ℓ = 13 … LMAX -------------------------------------------------
say(); say('Part B: does T_ℓ have the eigenvalue 0 on lane ℓ, class by class (weight ℓ − 3, Γ0(6), classes over the ℓ-adic '
          'integers, class_padic)?')
say('   classes named by (w₂, w₃) relative to ((2|ℓ), (3|ℓ)): "−,−" is the partition class; "−,+" is the one h₆ predicts')
hits = {}
for ell in [p for p in ps if 13 <= p <= LMAX] + AT:
    t = time.time(); k = ell - 3; row = []; N3 = 3
    Hp = F.hecke_matrices(k, ell ** N3, [('U', 2), ('U', 3), ('T', ell)])
    for name, (a2, a3) in (('−,−', (-1, -1)), ('−,+', (-1, 1)), ('+,−', (1, -1)), ('+,+', (1, 1))):
        w2, w3 = a2 * leg(2, ell), a3 * leg(3, ell)                    # the actual Atkin–Lehner signs
        res = F.class_padic(Hp, k, ell, N3, -w2, -w3)                 # U₂ = −w₂·2^((k−2)/2), U₃ likewise
        if res is None or res[5] == 0: row.append('%s: empty' % name); continue
        W, RW, r, Wn, Rn, dd, Np = res
        if Wn is None: row.append('%s: %d forms, precision too low' % (name, dd)); continue
        mp = ell ** Np
        TL = F.restrict(F.mat_mod(k + 1, k + 1, [int(v) % mp for v in Hp[('T', ell)].entries()], mp), Wn, Rn, dd)
        det0 = F.reduce_mat(TL, ell).det() == 0
        row.append('%s: %d forms%s' % (name, dd, ', EIGENVALUE 0' if det0 else ''))
        if det0: hits.setdefault(name, []).append(ell)
    say('   ℓ = %d (%.0f s): %s' % (ell, time.time() - t, ';  '.join(row)))
say(); say('non-ordinary forms found, by class: %s' % hits)
say("h₆'s non-ordinary primes among the ℓ tested: %s" % [p for p in nonord if 13 <= p <= LMAX or p in AT])
say('seconds: %.0f' % (time.time() - t0))
