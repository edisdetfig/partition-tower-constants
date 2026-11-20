# A dig, 13 Nov 2026 (version 2: the slope-S form's Hecke numbers on lane ℓ, --skip, Lagrange; README "The fading
# form: which family"). WHERE THE FADING NOTES LIVE (README "What fades"; a copy of partitions_flint_slopes.py, 12 Nov, aimed at
# the class forms of 13 to 31 and at a chosen slope S; predicted in the README first).
# The original header follows.
# A dig, 12 Nov 2026: THE CONSTANTS OF 5, 7 AND 11 ON FLINT. For these primes the tower's class is empty and the tower
# follows the smallest positive slope: s = 2, 1, 2 (the laws of Watson and Atkin), with α_ℓ = (3|ℓ)·u(−2)/ℓ^s, u the U_ℓ-
# number of the family of slope s in the tower's system on lane ℓ (STATEMENT.md, C; PARI runs: partitions_hida_laws.txt,
# partitions_hida_seven_14.txt). Here on FLINT alone (forms6_flint.py):
#   the system on lane ℓ: 5: a_Q ≡ Q²(1 + Q); 7: the weight-16 form, 11: the weight-28 form of the class (both the only form
#   of their class, found exactly);
#   at each weight k = −2 + (ℓ − 1)·ℓ·(x + 1), modulo ℓ^N′: the piece with the system's T_q numbers (q = 5, 7, 11, 13, not ℓ),
#   each lifted over the ℓ-adic integers (Fitting's splitting), inside it the forms new at 2 and 3 with the class's
#   U₂, U₃ values (signs carried along the weights), T_ℓ on them, its characteristic polynomial (Berkowitz), its Newton
#   polygon, and the one root of valuation s (Hensel); u from it; Newton's series at x = −1.
# PREDICTED before the run: at every weight exactly one root of slope s; the digits of (3|ℓ)·u(−2)/ℓ^s agree with the
# partition numbers' (5: 33 settled; 7: 14; 11: 9, 3, 2, 0) and with PARI's (5: 19; 7: 13 of 14; 11: 9, 3, 2, 0, 5, 8)
# as far as the series reaches. No partition number is used. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_slopes.py ℓ X N [--extra=E]
#          [--keep] [--threads=t]          (N′ = N + E; E defaults to 30)
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat
import forms6_flint as F
flint.ctx.threads = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--threads=')), 8)
args = [int(a) for a in sys.argv[1:] if a.isdigit()]; L, X, N = args[0], args[1], args[2]
EXTRA = next((int(a[8:]) for a in sys.argv[1:] if a.startswith('--extra=')), 30); NP = N + EXTRA
KEEP = '--keep' in sys.argv
S_ARG = next((int(a[8:]) for a in sys.argv[1:] if a.startswith('--slope=')), 2)
ANY = '--any' in sys.argv
SKIP = '--skip' in sys.argv
SYS = next(({int(kv.split(':')[0]): int(kv.split(':')[1]) for kv in a[6:].split(',')} for a in sys.argv[1:] if a.startswith('--sys=')), None)
OUT = os.path.join(HERE, 'fading_forms2_%d_s%d%s.txt' % (L, S_ARG, '_any' if ANY else ('_sys' if SYS else '')))
if KEEP and os.path.exists(OUT): OUT = OUT[:-4] + '_N%d.txt' % N
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None
SETTLED = {13: [10, 11, 6, 2, 10, 10], 17: [3, 9, 1, 13, 4, 3], 29: [23, 3, 10], 31: [2, 28, 13], 19: [], 23: []}   # tower_fading.txt: the fading notes' unit parts


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


S = S_ARG; QS = [q for q in (5, 7, 11, 13) if q != L]; t0 = time.time()
sign = lambda p, k: leg(p, L) ** (1 + (k - (L - 3)) // (L - 1))            # the class's sign of U_p at weight k
say('ℓ = %d, slope %d, weights −2 + %d·(x + 1) for x = 0..%d, modulo ℓ^%d (N = %d and %d spare); FLINT %s' % (
    L, S, (L - 1) * L, X, NP, N, EXTRA, flint.__version__))

# ---- the system on lane ℓ ----
if True:
    kr = L - 3; Kz, dr = F.class_basis_exact(kr, sign(2, kr), sign(3, kr))
    Xr, Rr = F.normalized_basis(nmod_mat(kr + 1, dr, [int(Kz[i, j]) % L for i in range(kr + 1) for j in range(dr)], L), dr, L)
    Hr = F.hecke_matrices(kr, L, [('T', q) for q in QS])
    if dr != 1: say('   the class at weight %d has %d forms; stopped' % (kr, dr)); sys.exit()
    ref = {q: int(F.restrict(Hr[('T', q)], Xr, Rr, 1)[0, 0]) for q in QS}
    how = 'the weight-%d form (U₂ = %+d·2^%d, U₃ = %+d·3^%d), the only one of its class' % (kr, sign(2, kr), (kr - 2) // 2, sign(3, kr), (kr - 2) // 2)
say('the system on lane ℓ: %s; T_q numbers %s' % (how, ref))
if SYS:
    say('   --sys: the piece is taken instead with the Hecke numbers %s (the slope-%d form found by --any)' % (SYS, S))
    class_ref = dict(ref); ref = dict(SYS)
else:
    class_ref = dict(ref)


def compose(Y, R, Z, Rz):
    return Y * Z, [R[j] for j in Rz]


# ---- the members ----
ws = []; precs = []; xs = []
for x in range(X + 1):
    k = -2 + (L - 1) * L * (x + 1); t = time.time(); m = L ** NP; n = k + 1
    ops = [('T', q) for q in QS] + [('T', L), ('U', 2), ('U', 3)]
    H = F.hecke_matrices(k, m, ops)
    Y, R, r = None, None, n                                    # the piece: the system's numbers, one T_q after another
    dims = []
    for q in ([] if ANY else QS):                              # --any: every system on lane ℓ, only the class's signs below
        Tq = H[('T', q)] if Y is None else F.restrict(H[('T', q)], Y, R, r)
        Z, Rz, r2 = F.kernel_part(Tq - F.identity(r, m) * ref[q], L, NP)
        if r2 == 0: Y = None; r = 0; break
        if Y is None: Y, R = Z, Rz
        else: Y, R = compose(Y, R, Z, Rz)
        r = r2; dims.append(r)
    if r == 0: say('   weight %d: no form with the system\'s numbers; stopped' % k); break
    c2 = sign(2, k) * 2 ** ((k - 2) // 2) % m; c3 = sign(3, k) * 3 ** ((k - 2) // 2) % m
    if Y is None: U2 = H[('U', 2)]; Z, Rz, r2 = F.kernel_part(U2 - F.identity(r, m) * c2, L, NP); Y, R = Z, Rz; r = r2; Z, Rz = F.identity(r, m), list(range(r))
    else: U2 = F.restrict(H[('U', 2)], Y, R, r); Z, Rz, r2 = F.kernel_part(U2 - F.identity(r, m) * c2, L, NP)
    if r2 == 0: say('   weight %d: no form with the class\'s U₂ value; stopped' % k); break
    Y, R = compose(Y, R, Z, Rz); r = r2
    U3 = F.restrict(H[('U', 3)], Y, R, r); Z, Rz, r2 = F.kernel_part(U3 - F.identity(r, m) * c3, L, NP)
    if r2 == 0: say('   weight %d: no form with the class\'s U₃ value; stopped' % k); break
    Y, R = compose(Y, R, Z, Rz); r = r2
    U2W = F.restrict(H[('U', 2)], Y, R, r); U3W = F.restrict(H[('U', 3)], Y, R, r); TLW = F.restrict(H[('T', L)], Y, R, r)
    Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)] + \
         [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(r)] for i in range(r)]
    ker, N1 = F.padic_kernel(Bm, 2 * r, r, L, NP); d = len(ker); m1 = L ** N1
    Xn = [[ker[c][i] % m1 for c in range(d)] for i in range(r)]  # the new forms, in the piece's coordinates (r × d)
    Rn = F.pivot_rows(nmod_mat(r, d, [Xn[i][c] % L for i in range(r) for c in range(d)], L), L)
    XR = F.inverse_mod(F.mat_mod(d, d, [Xn[Rn[a]][b] for a in range(d) for b in range(d)], m1), L, N1)
    TX = F.mat_mod(r, r, [int(v) % m1 for v in TLW.entries()], m1) * F.mat_mod(r, d, [Xn[i][c] for i in range(r) for c in range(d)], m1)
    Mt = XR * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m1)
    chi = F.berkowitz([[int(Mt[i, j]) for j in range(d)] for i in range(d)], m1)
    lam, pr, slopes = F.slope_root(chi, L, S, N1)
    poly = ', '.join('%s×%d' % (sl, ln) for sl, ln in slopes)
    if lam is None:
        say('   weight %d: the piece %s, new at 2 and 3: %d (to ℓ^%d); slopes %s: no single root of slope %d; %s' % (k, dims + [r], d, N1, poly, S, 'skipped' if SKIP else 'stopped'))
        if SKIP: continue
        break
    mp = L ** (pr - S); w = (lam // L ** S) % mp; e = k - 1 - 2 * S          # u = ℓ^S·w, w² − (λ/ℓ^S)·w + ℓ^e = 0
    for _ in range(pr + 2): w = ((lam // L ** S) - (L ** e % mp) * pow(w, -1, mp)) % mp if e < pr - S else (lam // L ** S) % mp
    ws.append(w); precs.append(pr - S); xs.append(x)
    # the slope-S form's Hecke numbers on lane ℓ: its eigenvector in the new forms, then T_q on it
    mpr = L ** min(pr, N1); Bk = [[(int(Mt[i, j]) - (lam if i == j else 0)) % mpr for j in range(d)] for i in range(d)]
    kv, Nk = F.padic_kernel(Bk, d, d, L, min(pr, N1))
    sysq = {}
    if len(kv) == 1:
        wv = kv[0]; i0 = next(i for i in range(d) if wv[i] % L)
        for q in QS:
            TQW = F.restrict(H[('T', q)], Y, R, r); TQp = F.mat_mod(r, r, [int(v) % m1 for v in TQW.entries()], m1)
            TXq = TQp * F.mat_mod(r, d, [Xn[i][c] for i in range(r) for c in range(d)], m1)
            Mq = XR * F.mat_mod(d, d, [int(TXq[Rn[a], b]) for a in range(d) for b in range(d)], m1)
            vq = sum(int(Mq[i0, j]) * wv[j] for j in range(d)) % L
            sysq[q] = vq * pow(wv[i0] % L, -1, L) % L
        twists = [j for j in range(L - 1) if all(sysq[q] == pow(q, j, L) * class_ref[q] % L for q in QS)]
        say("      its Hecke numbers on lane ℓ: %s; the class's: %s; θ-twists j with a_q ≡ q^j·a_q(class) for every q: %s" % (sysq, class_ref, twists or 'none'))
    else:
        say('      (its eigenvector was not found alone: %d kernel vectors)' % len(kv))
    say('   weight %d (dimension %d): the piece %s, new at 2 and 3 with the class\'s signs: %d (to ℓ^%d); slopes of T_ℓ there: %s; '
        'u/ℓ^%d on lane ℓ: %d, known to ℓ^%d; (3|ℓ)·u/ℓ^%d, digits lowest first: %s  (%.0f s)' % (k, n, dims + [r], d, N1, poly, S, w % L, pr - S, S,
        [((leg(3, L) * w % L ** (pr - S)) // L ** i) % L for i in range(min(8, pr - S))], time.time() - t))

# ---- weight −2 ----
if len(ws) >= 2:
    P = min(precs); mP = L ** P; v = 0
    for i, xi in enumerate(xs):
        num, den = 1, 1
        for j, xj in enumerate(xs):
            if j != i: num *= (-1 - xj); den *= (xi - xj)
        v = (v + ws[i] % mP * num * pow(den % mP, -1, mP)) % mP
    val = 'Lagrange through x = %s' % xs
    a = leg(3, L) * v % L ** P; K = min(P, len(xs))
    say('the differences are divisible by ℓ to %s' % val)
    dg = [(a // L ** j) % L for j in range(P)]
    same = next((i for i in range(min(len(dg), len(SETTLED[L]))) if dg[i] != SETTLED[L][i]), min(len(dg), len(SETTLED[L])))
    say('(3|ℓ)·u(−2)/ℓ^%d, digits lowest first (the first %d are the reach of %d weights): %s' % (S, K, len(xs), dg[:max(K, len(SETTLED[L]))]))
    say('the fading note measured from partition numbers (tower_fading.txt): %s  →  the first %d digits agree' % (SETTLED[L], same))
say(); say('seconds: %.0f' % (time.time() - t0))
