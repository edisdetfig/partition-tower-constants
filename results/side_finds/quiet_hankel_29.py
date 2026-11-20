# A dig, 13 Nov 2026: 29'S QUIET NOTE BY THE SHORTCUT (README "29's quiet note by the shortcut, beside the forms run";
# predicted there first). As quiet_hankel_check.py, at 29 and deeper: H_b(n) = p(29^b·n + δ_b) for b = 1, 3, …, 11 (n = 0 …
# 7, as far as N ≤ 4·10¹⁶), exact partition numbers from FLINT kept modulo 29^30.
#   - the shortcut: D_b = H_(b+4)·H_b − H_(b+2)²; the ratio D_(b+2)/D_b, its valuation and the digits of its unit part that
#     all n share; divided by α (partitions_flint_family_29.txt, 30 digits), β/29;
#   - the cross-check: E_b = H_(b+2) − α·H_b and the ratio E_(b+2)/E_b, the same digits.
# COUNTED. Runs at idle priority beside the forms run.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/quiet_hankel_29.py [--keep] [--threads=t]
import os, sys, re, time
HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
flint.ctx.threads = next((int(a[10:]) for a in sys.argv[1:] if a.startswith('--threads=')), 4)
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'quiet_hankel_29.txt'), 'x', encoding='utf-8') if KEEP else None
D, NMAX, L = 30, 4 * 10 ** 16, 29


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


PC = {}
def c(M):
    N = (M + 1) // 24
    if N not in PC: PC[N] = int(flint.fmpz(N).partitions_p())
    return PC[N]


def val(x, cap):
    x %= L ** cap
    if x == 0: return cap
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


def shared(rats):
    """rats: (valuation, unit, modulus) per n. The digits of the unit part that all n share."""
    k = 0; mpmin = min(r[2] for r in rats)
    while L ** (k + 1) <= mpmin and len({r[1] % L ** (k + 1) for r in rats}) == 1: k += 1
    return k, [(rats[0][1] // L ** i) % L for i in range(k)]


def ratios(X, b, cap):
    out = []
    for n in X[b + 2]:
        if n not in X[b]: continue
        x0, x2 = X[b][n], X[b + 2][n]; v0, v2 = val(x0, cap), val(x2, cap)
        mp = L ** max(1, min(cap - v0, cap - v2) - v0)
        out.append((v2 - v0, (x2 // L ** v2) * pow((x0 // L ** v0) % mp, -1, mp) % mp, mp))
    return out


t0 = time.time(); m = L ** D
al = [int(x) for x in re.search(r'a constant on lane ℓ, digits lowest first: \[([^\]]*)\]', open(os.path.join(HERE, 'partitions_flint_family_29.txt'), encoding='utf-8').read()).group(1).split(',')]
alpha = sum(x * L ** i for i, x in enumerate(al[:D])) % m
H = {}
for b in range(1, 12, 2):
    kk = (-pow(L ** b, -1, 24)) % 24; row = {}
    for n in range(8):
        M = L ** b * (24 * n + kk)
        if (M + 1) // 24 > NMAX: break
        row[n] = c(M) % m
    H[b] = row
    say('level %2d: %d values of n (%.0f s)' % (b, len(row), time.time() - t0))
Dv = {b: {n: (H[b + 4][n] * H[b][n] - H[b + 2][n] ** 2) % m for n in H[b + 4] if n in H[b] and n in H[b + 2]} for b in range(1, 8, 2) if b + 4 in H}
Ev = {b: {n: (H[b + 2][n] - alpha * H[b][n]) % m for n in H[b + 2] if n in H[b]} for b in range(1, 10, 2) if b + 2 in H}
say()
for b in sorted(Dv):
    if b + 2 not in Dv: continue
    rats = ratios(Dv, b, D)
    if len(rats) < 2: continue
    k, dg = shared(rats); u = sum(d * L ** i for i, d in enumerate(dg)); mk = L ** max(k, 1)
    beta = [((u * pow(alpha % mk, -1, mk)) % mk // L ** i) % L for i in range(k)] if k else []
    say('the shortcut, levels %d … %d (%d values of n): D(next)/D has valuation %s; α·β/29 shared by all n: %s (%d digits); so β/29 = %s' % (
        b, b + 6, len(rats), sorted({r[0] for r in rats}), dg, k, beta))
for b in sorted(Ev):
    if b + 2 not in Ev: continue
    rats = ratios(Ev, b, D)
    if len(rats) < 2: continue
    k, dg = shared(rats)
    say('the leftover E, levels %d → %d (%d values of n): E(next)/E has valuation %s; β/29 shared by all n: %s (%d digits)' % (
        b, b + 2, len(rats), sorted({r[0] for r in rats}), dg, k))
say(); say('the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
