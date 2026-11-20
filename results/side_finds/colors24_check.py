# A dig, 12 Nov 2026: 24 COLOURS (1/Δ), THE FORMS' CONSTANTS AGAINST THE EXACT TOWERS. The tower is
# H_b(n) = p₂₄(ℓᵇn + 1) = c(ℓᵇn), c the coefficients of 1/Δ, one level per step (1/Δ | U_ℓᵇ; no Shimura step at integral
# weight). forms1_family_w-12_ℓ.txt gives, from level-1 forms alone, the polynomial P(x) of the constants at weight −12:
# the Eisenstein family's 1 times the ordinary cusp forms' families.
# PREDICTED before the run (derived, PROOF.md's classical route at integral weight): Σ_i [xⁱ]P · H_(b+i)(n) ≡ 0 on lane ℓ
# for every n at the first levels, and on ever deeper lanes as b grows; at 17, where only the Eisenstein family is
# ordinary (M₄ = ⟨E₄⟩), simply H_(b+1)(n) ≡ H_b(n): the 24-colour tower at 17 stops changing. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/colors24_check.py <tower json> … [--keep]
import os, sys, re, json
HERE = os.path.dirname(os.path.abspath(__file__)); KEEP = '--keep' in sys.argv
SRCS = [a for a in sys.argv[1:] if a.endswith('.json')]
fh = open(os.path.join(HERE, 'colors24_check.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


for src in SRCS:
    data = json.load(open(os.path.join(HERE, src), encoding='utf-8'))
    say('towers: %s (n < %d)' % (src, data['N']))
    for Ls, rec in sorted(data['primes'].items(), key=lambda kv: int(kv[0])):
        L = int(Ls); f = os.path.join(HERE, 'forms1_family_w-12_%d.txt' % L)
        if not os.path.exists(f): say('   ℓ = %d: no family file' % L); continue
        mt = re.search(r'on lane ℓ\^(\d+): \[([^\]]*)\];  on lane ℓ:', open(f, encoding='utf-8').read())
        KF = int(mt.group(1)); P = [int(v) for v in mt.group(2).split(',')]; D = len(P) - 1
        K = min(KF, rec['K']); M = L ** K; lv = {int(b): [x % M for x in v] for b, v in rec['levels'].items()}
        rows = []
        for b in sorted(lv):
            if b + D not in lv: continue
            nn = min(len(lv[b + i]) for i in range(D + 1))
            deps = []
            for n in range(nn):
                w = sum(P[i] * lv[b + i][n] for i in range(D + 1)) % M
                deps.append(next((e for e in range(K) if w % L ** (e + 1)), K))
            rows.append('from level %d: %s' % (b, deps))
        say('   ℓ = %d, %d constant%s (P on lane ℓ %s): the law\'s ℓ-adic depth for n = 0, 1, …  %s' % (
            L, D, '' if D == 1 else 's', [c % L for c in P], ';  '.join(rows) or 'no level span below N'))
