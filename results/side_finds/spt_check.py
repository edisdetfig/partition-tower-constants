# A dig, 12 Nov 2026: THE spt TOWERS, NORMALIZED. spt_tower_N16000000_*.txt: for every prime from 11 to 43 except 13, two
# levels multiply the tower by (3|ℓ) on ever deeper lanes. At 5, 7, 13 the values carry a law (Garvan), so there the test
# needs H′_b = H_b / ℓ^e(b) first, e(b) = ⌊(b+1)/2⌋ (from memory of Garvan's theorem; measured below).
# PREDICTED before the run (derived): the law primes follow the Eisenstein family's OTHER root (U_ℓ = ℓ at weight 2, slope 1)
# because Andrews' congruences kill the ordinary part there; divided by the law, the constant is again (3|ℓ)·ℓ/ℓ = (3|ℓ):
# H′_(b+2) ≡ (3|ℓ)·H′_b on ever deeper lanes, at 5, 7 and 13 as at every other prime. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/spt_check.py [--keep]
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'spt_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


d = json.load(open(os.path.join(HERE, 'spt_tower_N16000000_11_17_19_23_29_31_37_41_43_5_7_13.json'), encoding='utf-8'))
for Ls in ('5', '7', '13'):
    L = int(Ls); rec = d['primes'][Ls]; K = rec['K']; lv = {int(b): v for b, v in rec['levels'].items()}; s3 = leg(3, L)
    def val(x):
        if x % L ** K == 0: return K
        v = 0
        while x % L == 0: x //= L; v += 1
        return v
    e = {b: min(val(x) for x in v) for b, v in lv.items()}
    say('ℓ = %d, (3|ℓ) = %+d: the smallest valuation at each level %s (⌊(b+1)/2⌋: %s)' % (L, s3, e, {b: (b + 1) // 2 for b in lv}))
    for b in sorted(lv):
        if b + 2 not in lv: continue
        eb, ec = (b + 1) // 2, (b + 3) // 2
        P = K - ec; M = L ** P
        if P <= 0: continue
        nn = min(len(lv[b]), len(lv[b + 2]))
        deps = []
        for i in range(nn):
            x = (lv[b][i] // L ** eb) % M; y = (lv[b + 2][i] // L ** ec) % M
            w = (y - s3 * x) % M
            deps.append(next((t for t in range(P) if w % L ** (t + 1)), P))
        say('   H′_%d − (3|ℓ)·H′_%d (divided by ℓ^%d and ℓ^%d): ℓ-adic depth %s (of %d)' % (b + 2, b, ec, eb, deps, P))
