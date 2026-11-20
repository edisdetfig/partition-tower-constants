# A dig, 14 Nov 2026 (README "What sets the settling rate", the check at 7 and 11): tower_fading.py for 7 and 11, where
# Watson's and Atkin's laws first divide the tower by a power of ℓ. H(M) = p((M + 1)/24) / ℓ^e(b) with M = ℓ^b·(24n + k)
# at odd levels b, e(b) = ⌊b/2⌋ + 1 at 7 (Watson) and b at 11 (Atkin). E_b(n) = H(ℓ²M) − α·H(M), α from
# partitions_flint_slopes_ℓ.txt, modulo ℓ^D with D the digits α is known to (20 at 7, 12 at 11). Printed for each level:
# the valuations v(E_b(n)) over n = 0 … 7, and the ratio E_(b+2)(n)/E_b(n): its valuation and its unit part's digits
# shared by all n. FLINT's exact partition numbers, at idle priority. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/tower_fading_law.py [--keep]
import os, sys, re, time
HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
flint.ctx.threads = 6
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'tower_fading_law.txt'), 'x', encoding='utf-8') if KEEP else None
NMAX = 10 ** 14


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


PC = {}


def p(N):
    if N not in PC: PC[N] = int(flint.fmpz(N).partitions_p())
    return PC[N]


def val(x, L, cap):
    if x == 0: return cap
    v = 0
    while x % L == 0: x //= L; v += 1
    return v


t0 = time.time()
for L, law in ((7, lambda b: b // 2 + 1), (11, lambda b: b)):
    txt = open(os.path.join(HERE, 'partitions_flint_slopes_%d.txt' % L), encoding='utf-8').read()
    al = [int(x) for x in re.search(r'u\(−2\)/ℓ\^\d, digits lowest first[^:]*: \[([^\]]*)\]', txt).group(1).split(',')]
    D = len(al); m = L ** D; aInt = sum(x * L ** i for i, x in enumerate(al))
    say(); say('ℓ = %d: α to %d digits (%s…); the tower divided by ℓ^e(b), e(b) = %s' % (
        L, D, al[:6], '⌊b/2⌋ + 1 (Watson)' if L == 7 else 'b (Atkin)'))
    E = {}
    for b in range(1, 40, 2):
        kk = (-pow(L ** b, -1, 24)) % 24; row = {}
        for n in range(8):
            M = L ** b * (24 * n + kk)
            if (L * L * M + 1) // 24 > NMAX: break
            N0, N2 = (M + 1) // 24, (L * L * M + 1) // 24
            P0, P2 = p(N0), p(N2)
            assert P0 % L ** law(b) == 0 and P2 % L ** law(b + 2) == 0, 'the law fails at ℓ = %d, b = %d, n = %d' % (L, b, n)
            row[n] = (P2 // L ** law(b + 2) - aInt * (P0 // L ** law(b))) % m
        if not row: break
        E[b] = row
    for b in sorted(E):
        vs = [val(E[b][n], L, D) for n in sorted(E[b])]
        line = '   level %2d: v(E) over n = %s' % (b, vs)
        if b + 2 in E:
            rats = []
            for n in [n for n in E[b + 2] if n in E[b]]:
                e0, e2 = E[b][n], E[b + 2][n]; v0, v2 = val(e0, L, D), val(e2, L, D)
                if v0 >= D or v2 >= D: continue
                u0, u2 = e0 // L ** v0, e2 // L ** v2
                mp = L ** max(1, min(D - v0, D - v2))
                rats.append((v2 - v0, u2 * pow(u0 % mp, -1, mp) % mp, mp))
            if rats:
                vr = sorted({r[0] for r in rats}); k = 0; mpmin = min(r[2] for r in rats)
                while L ** (k + 1) <= mpmin and len({r[1] % L ** (k + 1) for r in rats}) == 1: k += 1
                dg = [(rats[0][1] // L ** i) % L for i in range(min(8, k + 1))]
                line += ';  E(b + 2)/E(b): valuation %s, unit part digits %s… shared by all %d values of n: %d digit%s' % (
                    vr, dg, len(rats), k, '' if k == 1 else 's')
        say(line)
say(); say('the largest partition number used: p(%d);  seconds: %.0f' % (max(PC), time.time() - t0))
