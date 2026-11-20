# A check, 14 Nov 2026 (README "The quiet notes of 7 and 11, digit by digit"): carry 7's and 11's quiet forms to weight −2
# and compare with the fading part measured on partition numbers (tower_fading_law.txt).
# Rows from fading_forms5.py's .jsonl files. At 7 the quiet form (system (5, 6, 0)) has slope 4 at weights 40, 46, 52
# (x = 6, 7, 8) and slope 5 at 58 … 76 (x = 9 … 12; its zeros on the first digit). Each row's root r = u/7^S is put on one
# scale, u/7^4 = 7^(S − 4)·r. At 11 the quiet form (the system present at x = 4) has slope 3 at x = 4 … 11 (at 12 it is
# not of slope 3 in its system: a zero on the first digit, as are 13 and 14).
# The weights used differ by less than ℓ, so the curve's denominators are units (no division by ℓ).
# Plain Lagrange at x = −1 (these forms are not 0 at the row's start), then (3|ℓ)·u/ℓ^s, digits lowest first.
# Printed with the agreement of the carries from every one weight fewer. DERIVED from COUNTED forms.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/fading_quiet_fit_7_11.py
import os, json, re, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
fh = open(os.path.join(HERE, 'fading_quiet_fit_7_11.txt'), 'x', encoding='utf-8')


def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); fh.write(s + '\n'); fh.flush()


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def rows(fn): return [json.loads(l) for l in open(os.path.join(HERE, fn), encoding='utf-8') if l.strip()]


def lag(xs, f, m):
    v = 0
    for a in xs:
        nu = de = 1
        for b in xs:
            if b != a: nu *= (-1 - b); de *= (a - b)
        v = (v + f[a] * nu * pow(de % m, -1, m)) % m
    return v


def dig(v, L, n): return [(v // L ** i) % L for i in range(n)]


def agree(u, w, L, n):
    k = 0
    while k < n and (u // L ** k) % L == (w // L ** k) % L: k += 1
    return k


tf = open(os.path.join(HERE, 'tower_fading_law.txt'), encoding='utf-8').read()
cases = [(7, 4, ['fading_forms5_7_s4_x6-17.jsonl', 'fading_forms5_7_s5_x9-12.jsonl'], list(range(6, 13))),
         (11, 3, ['fading_forms5_11_s3_x4-15.jsonl'], list(range(4, 12)))]
for L, SQ, files, xs in cases:
    allr = [r for fn in files for r in rows(fn)]
    sysl = next(r['system'] for r in allr if r['x'] == xs[0] and r['slope'] == SQ)
    f, prec = {}, {}
    for r in allr:
        if r['system'] != sysl or r['deg'] != 1 or r['x'] not in xs: continue
        root = (-r['G'][0]) % L ** r['prec']                 # u/ℓ^S, known to ℓ^prec
        f[r['x']] = root * L ** (r['slope'] - SQ); prec[r['x']] = r['prec'] + r['slope'] - SQ
    assert sorted(f) == xs, 'missing weights at ℓ = %d: %s' % (L, sorted(set(xs) - set(f)))
    P = min(prec.values()); m = L ** P; s3 = leg(3, L)
    full = s3 * lag(xs, f, m) % m
    ag = min(agree(s3 * lag(list(s), f, m) % m, full, L, P) for s in itertools.combinations(xs, len(xs) - 1))
    say(); say('ℓ = %d: the quiet form in the system %s, slope %d, weights x = %s (slopes there: %s), modulo ℓ^%d' % (
        L, sysl, SQ, xs, [r['slope'] for r in sorted([r for r in allr if r['system'] == sysl and r['x'] in xs], key=lambda r: r['x'])], P))
    say('   (3|ℓ)·u(−2)/ℓ^%d from the forms, digits lowest first: %s' % (SQ, dig(full, L, min(P, ag + 1))))
    say('   the carries from every %d of the %d weights agree with it on the first %d digits (digit %d likely)' % (
        len(xs) - 1, len(xs), ag, ag + 1))
    blk = tf.split('ℓ = %d:' % L)[1].split('ℓ = ')[0]
    best = max(((int(k), [int(x) for x in d.split(',')]) for d, k in re.findall(r'unit part digits \[([^\]]*)\]… shared by all \d+ values of n: (\d+) digit', blk)), key=lambda t: t[0])
    k, dg = best
    say('   the fading part on partition numbers (tower_fading_law.txt), its best level: %d digits shared by all n: %s' % (k, dg[:k]))
    fd = dig(full, L, k)
    n_ag = 0
    while n_ag < k and fd[n_ag] == dg[n_ag]: n_ag += 1
    say('   agreement: %d of %d digits' % (n_ag, k))
fh.close()
