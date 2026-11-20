# serre_weights.py, 14 Nov 2026 (README "Why level 6ℓ gains one or two forms per step").
# The counts of level6l_counts.py and level6l_counts2.py, predicted from small weights at level 6 alone.
#   1. GL2(F_ℓ). The composition factors of the principal series Ind(a^s·d^t) and of Sym^n (Glover's recursion: Sym^n has
#      Sym^(n−ℓ−1) ⊗ det inside, with the principal series Ind(a^n) as quotient). Each is checked by Brauer characters on
#      every semisimple class, split and elliptic. Also checked: Ind_B^G Sym^n = Σ_i Ind(a^(n−i)·d^i), by characters.
#   2. The seeds. Take every Serre weight σ = Sym^r ⊗ det^m with r + 2m ≡ ℓ − 5 (mod ℓ − 1), the row's central character.
#      For each, count the forms of weight r + 2 at level 6 in the class's system twisted by ω^(−m)
#      (U_2 → 2^(−m)·U_2, U_3 → 3^(−m)·U_3, T_q → q^(−m)·T_q). Two counts:
#        new: new at 2 and 3 with those signs exactly (the Shimura curve's forms);
#        all: the whole generalized eigenspace mod ℓ in M_(r+2)(Γ0(6)).
#   3. The predictions. Level 6ℓ is the sum of the seeds over the factors of Ind_B^G Sym^(k−2); level 6 is the sum over the
#      factors of Sym^(k−2). They are compared with every system at every weight in the .jsonl files given.
#      One correction, derived after the trial at 7 (README): on the compact Shimura curve the forms counted are minus the
#      Euler characteristic, which adds up exactly over the factors. Only the factors det^m′ (r = 0) have invariants. Each
#      det^m′ whose numbers (T_q ≡ q^m′·(1 + q), U_2 ≡ 2^m′, U_3 ≡ 3^m′) are the class's takes one form away
#      (h⁰ = h² = 1, h¹ = 0). This applies to the class counts. The whole-eigenspace counts (Γ0(6), with Eisenstein series)
#      are predicted only for systems that are not of Eisenstein type (T_q ≡ q^a + q^b).
# COUNTED (the seeds), DERIVED (the predictions), compared with COUNTED (the .jsonl files).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/serre_weights.py ℓ a.jsonl b.jsonl …
#              [--keep]
import os, sys, time, json, ast
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat
import forms6_flint as F

flint.ctx.threads = 2
args = [a for a in sys.argv[1:] if not a.startswith('--')]
L = int(args[0]); FILES = args[1:]; KEEP = '--keep' in sys.argv
QS = [q for q in (5, 7, 11, 13) if q != L]
LM = L - 1; NPW = 10
OUT = os.path.join(HERE, 'serre_weights_%d.txt' % L)
if KEEP and os.path.exists(OUT): sys.exit('stopped before computing anything: %s already exists' % OUT)
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def sign(p, k): return leg(p, L) ** (1 + (k - (L - 3)) // (L - 1))      # the class's sign of U_p (fading_forms5.py)


# ---- the counts to predict ----
rows = {}
for fn in FILES:
    p = fn if os.path.isabs(fn) else os.path.join(HERE, fn)
    for line in open(p, encoding='utf-8'):
        r = json.loads(line)
        assert r['l'] == L
        old = rows.get(r['x'])
        if old:                                             # the same weight in two files: the counts must agree
            for s, v in r['systems'].items():
                o = old['systems'].get(s, {'d6': 0, 'd6l': 0})
                assert (o['d6'], o['d6l']) == (v['d6'], v['d6l']), 'two runs disagree at x = %d, %s' % (r['x'], s)
            if 'm6' in next(iter(r['systems'].values()), {}): rows[r['x']] = r
        else: rows[r['x']] = r
XS = sorted(rows); KS = {x: rows[x]['k'] for x in XS}
say('ℓ = %d: the counts of %s, at x = %s (weights %s)' % (L, ', '.join(os.path.basename(f) for f in FILES), XS, [KS[x] for x in XS]))

# ---- 1. GL2(F_ℓ): composition factors, checked by Brauer characters ----
def jh_ind(s, t):
    s, t = s % LM, t % LM
    if s == t: return [(0, t), (LM, t)]
    r = (s - t) % LM
    return [(r, t), (LM - r, s)]


def jh_sym(n):
    if n <= LM: return [(n, 0)]
    out = jh_ind(n, 0)
    if n >= L + 1: out += [(r, (m + 1) % LM) for r, m in jh_sym(n - L - 1)]
    return out


LL = L * L - 1
UU, VV = [a.ravel() for a in np.meshgrid(np.arange(LM), np.arange(LM), indexing='ij')]
WW = np.array([w for w in range(LL) if w % (L + 1)])
z = lambda e, N: np.exp(2j * np.pi * (np.asarray(e) % N) / N)
def ch_serre(r, m):
    return np.concatenate([sum(z((m + i) * UU + (m + r - i) * VV, LM) for i in range(r + 1)),
                           sum(z(WW * ((L + 1) * m + i + L * (r - i)), LL) for i in range(r + 1))])
def ch_sym(n):
    return np.concatenate([sum(z(i * UU + (n - i) * VV, LM) for i in range(n + 1)),
                           sum(z(WW * (i + L * (n - i)), LL) for i in range(n + 1))])
def ch_ind(s, t):
    sp = np.where(UU != VV, z(s * UU + t * VV, LM) + z(s * VV + t * UU, LM), (L + 1) * z((s + t) * UU, LM))
    return np.concatenate([sp, np.zeros(len(WW))])
def same(a, b): return float(np.max(np.abs(a - b))) < 1e-6


t0 = time.time(); nmax = max(KS.values()) - 2
bad_ind = [(s, t) for s in range(LM) for t in range(LM) if not same(ch_ind(s, t), sum(ch_serre(*g) for g in jh_ind(s, t)))]
bad_sym = [n for n in range(nmax + 1) if not same(ch_sym(n), sum(ch_serre(*g) for g in jh_sym(n)))]
bad_big = [k for k in KS.values() if not same(sum(ch_ind(k - 2 - i, i) for i in range(k - 1)), ch_sym(k - 2) * ch_ind(0, 0))]
say('GL2(F_ℓ), by Brauer characters on %d split and %d elliptic classes (%.0f s):' % (LM * LM, len(WW), time.time() - t0))
say('   the factors of Ind(a^s·d^t), all %d characters: %s' % (LM * LM, 'right' if not bad_ind else 'WRONG at %s' % bad_ind[:5]))
say('   the factors of Sym^n, n = 0 … %d (Glover\'s recursion): %s' % (nmax, 'right' if not bad_sym else 'WRONG at n = %s' % bad_sym[:5]))
say('   Ind_B^G Sym^(k−2) = Σ_i Ind(a^(k−2−i)·d^i) = Sym^(k−2) ⊗ Ind(1), at every weight: %s' % ('right' if not bad_big else 'WRONG at %s' % bad_big))
if bad_ind or bad_sym or bad_big: say('STOP'); sys.exit(1)
cnt = {}
for i in range(LM):
    for g in jh_ind(LM - 4 - i, i): cnt[g] = cnt.get(g, 0) + 1
say('   one step adds the %d principal series Ind(a^(k−2−i)·d^i), i mod %d; in them each Serre weight with 1 ≤ r ≤ ℓ − 2 '
    'appears %s times, with r = 0 or ℓ − 1 %s times' % (LM, LM, sorted({v for (r, m), v in cnt.items() if 0 < r < LM}),
                                                         sorted({v for (r, m), v in cnt.items() if r in (0, LM)})))


# ---- 2. the seeds ----
def split_systems(Mq, d, Nn):
    m1 = L ** Nn; pieces = [((), F.identity(d, m1), list(range(d)), d)]
    for q in QS:
        newp = []
        for lab, Yp, Rp, rp in pieces:
            Mqr = F.restrict(Mq[q], Yp, Rp, rp)
            for f, ee in F.reduce_mat(Mqr, L).charpoly().factor()[1]:
                co = [int(c) for c in f.coeffs()]
                Bq = Mqr - F.identity(rp, m1) * ((-co[0]) % L) if len(co) == 2 else F.mat_poly(co, Mqr)
                Zm, Rz, r2 = F.kernel_part(Bq, L, Nn)
                if r2: newp.append((lab + (((-co[0]) % L) if len(co) == 2 else tuple(co),), Yp * Zm, [Rp[j] for j in Rz], r2))
        pieces = newp
    return {lab: rp for lab, Yp, Rp, rp in pieces}


def class_and_eigenspace(H, w, s2, s3, N):
    """At weight w on Γ0(6): ({system: forms new at 2 and 3 with U_2 = s2·2^((w−2)/2), U_3 = s3·3^((w−2)/2) exactly},
    {system: forms in the whole generalized eigenspace of those numbers mod ℓ})."""
    n = w + 1; m = L ** N; c2 = s2 * 2 ** ((w - 2) // 2) % m; c3 = s3 * 3 ** ((w - 2) // 2) % m
    Y, R, r = F.kernel_part(H[('U', 2)] - F.identity(n, m) * c2, L, N)
    if r == 0: return {}, {}
    Zm, Rz, r2 = F.kernel_part(F.restrict(H[('U', 3)], Y, R, r) - F.identity(r, m) * c3, L, N)
    if r2 == 0: return {}, {}
    Y, R = Y * Zm, [R[j] for j in Rz]; r = r2
    gen = split_systems({q: F.restrict(H[('T', q)], Y, R, r) for q in QS}, r, N)
    U2W = F.restrict(H[('U', 2)], Y, R, r); U3W = F.restrict(H[('U', 3)], Y, R, r)
    Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)] + \
         [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(r)] for i in range(r)]
    ker, N1 = F.padic_kernel(Bm, 2 * r, r, L, N); d = len(ker)
    if d == 0: return {}, gen
    assert N1 >= 1, 'the class is not exact at this precision'
    m1 = L ** N1
    Xm = F.mat_mod(r, d, [ker[c][i] % m1 for i in range(r) for c in range(d)], m1)
    Rn = F.pivot_rows(nmod_mat(r, d, [ker[c][i] % L for i in range(r) for c in range(d)], L), L)
    XR = F.inverse_mod(F.mat_mod(d, d, [ker[b][Rn[a]] % m1 for a in range(d) for b in range(d)], m1), L, N1)
    def on_new(T):
        TX = F.mat_mod(r, r, [int(v) % m1 for v in T.entries()], m1) * Xm
        return XR * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m1)
    return split_systems({q: on_new(F.restrict(H[('T', q)], Y, R, r)) for q in QS}, d, N1), gen


def twist(lab, m):
    """The system with every T_q multiplied by q^m (a factor's roots multiplied by q^m)."""
    out = []
    for q, v in zip(QS, lab):
        f = pow(q, m, L)
        out.append(v * f % L if isinstance(v, int) else tuple(v[i] * pow(f, len(v) - 1 - i, L) % L for i in range(len(v))))
    return tuple(out)


c2 = sign(2, L - 3) * pow(2, (L - 5) // 2, L) % L; c3 = sign(3, L - 3) * pow(3, (L - 5) // 2, L) % L
assert all((sign(2, k) * pow(2, (k - 2) // 2, L) - c2) % L == 0 and (sign(3, k) * pow(3, (k - 2) // 2, L) - c3) % L == 0
           for k in KS.values()), 'the class\'s U_2, U_3 are not constant mod ℓ along the row'
say(); say('the class mod ℓ, the same at every step of the row: U_2 ≡ %d, U_3 ≡ %d' % (c2, c3))
say('the seeds: for each Serre weight Sym^r ⊗ det^m (r + 2m ≡ ℓ − 5 mod ℓ − 1), the forms of weight r + 2 at level 6 with '
    'U_2 ≡ 2^(−m)·%d, U_3 ≡ 3^(−m)·%d, by system twisted back by q^m (new = new at 2 and 3; all = with the old forms)' % (c2, c3))
mu_new, mu_all = {}, {}; t0 = time.time()
for r in range(2, L, 2):
    w = r + 2
    H = F.hecke_matrices(w, L ** NPW, [('T', q) for q in QS] + [('U', 2), ('U', 3)])
    for m in [m for m in range(LM) if (r + 2 * m - (L - 5)) % LM == 0]:
        u2 = pow(2, -m, L) * c2 % L; u3 = pow(3, -m, L) * c3 % L
        s2 = u2 * pow(2, -(r // 2), L) % L; s3 = u3 * pow(3, -(r // 2), L) % L
        assert s2 in (1, L - 1) and s3 in (1, L - 1)
        s2, s3 = (1 if s2 == 1 else -1), (1 if s3 == 1 else -1)
        ex, gen = class_and_eigenspace(H, w, s2, s3, NPW)
        for lab, c in ex.items(): mu_new.setdefault(twist(lab, m), {})[(r, m)] = c
        for lab, c in gen.items(): mu_all.setdefault(twist(lab, m), {})[(r, m)] = c
        say('   Sym^%d ⊗ det^%d (weight %d, signs %s%s): new %s; all %s' % (
            r, m, w, '+' if s2 > 0 else '−', '+' if s3 > 0 else '−',
            {repr(twist(l, m)): c for l, c in ex.items()} or '—', {repr(twist(l, m)): c for l, c in gen.items()} or '—'))
say('   (%.0f s)' % (time.time() - t0))


# ---- 3. the predictions, against every count ----
def eis_twists(lab):
    """The m′ for which det^m′ carries the class's numbers: T_q ≡ q^m′·(1 + q) for q in QS, U_2 ≡ 2^m′, U_3 ≡ 3^m′."""
    return {mp for mp in range(LM) if all(isinstance(v, int) and v == pow(q, mp, L) * (1 + q) % L for q, v in zip(QS, lab))
            and pow(2, mp, L) == c2 and pow(3, mp, L) == c3}


def eis_type(lab):
    """Is the system ≡ q^a + q^b on lane ℓ for every q in QS (a reducible representation mod ℓ)?"""
    return any(all(isinstance(v, int) and v == (pow(q, a, L) + pow(q, b, L)) % L for q, v in zip(QS, lab))
               for a in range(LM) for b in range(a, LM))


def predict(mu, lab, k, correct):
    s = mu.get(lab, {}); E = eis_twists(lab) if correct else set()
    val = lambda g: s.get(g, 0) - (1 if g[0] == 0 and g[1] in E else 0)
    return sum(val(g) for i in range(k - 1) for g in jh_ind(k - 2 - i, i)), sum(val(g) for g in jh_sym(k - 2))


def gain(mu, lab, correct):
    s = mu.get(lab, {}); E = eis_twists(lab) if correct else set()
    return sum(s.get(g, 0) - (1 if g[0] == 0 and g[1] in E else 0) for i in range(LM) for g in jh_ind(LM - 4 - i, i))


found = {}
for x in XS:
    for s, v in rows[x]['systems'].items(): found.setdefault(ast.literal_eval(s), {})[x] = v
labs = sorted(set(found) | set(mu_new) | set(mu_all), key=repr)
has_m = all('m6' in v for x in XS for v in rows[x]['systems'].values())
say(); say('THE PREDICTIONS against the counts (level 6ℓ, then level 6; "pred/counted"):')
tally = {'class': [0, 0], 'eigenspace': [0, 0]}
for lab in labs:
    seeds = sorted(mu_new.get(lab, {}).items()); seeds_all = sorted(mu_all.get(lab, {}).items())
    et = eis_type(lab); E = sorted(eis_twists(lab))
    say('   system %s%s: seeds new %s, all %s; gain per step predicted: class %d%s' % (
        repr(lab), ' (Eisenstein type%s)' % ('; det^%s carries its numbers' % E if E else '') if et else '',
        ', '.join('weight %d·ω^%d ×%d' % (r + 2, m, c) for (r, m), c in seeds) or '—',
        ', '.join('weight %d·ω^%d ×%d' % (r + 2, m, c) for (r, m), c in seeds_all) or '—', gain(mu_new, lab, True),
        '' if et else ', eigenspace %d' % gain(mu_all, lab, False)))
    for kind, mu, key6l, key6 in (('class', mu_new, 'd6l', 'd6'), ('eigenspace', mu_all, 'm6l', 'm6')):
        if kind == 'eigenspace' and (not has_m or et): continue
        parts = []
        for x in XS:
            D, d = predict(mu, lab, KS[x], kind == 'class'); v = found.get(lab, {}).get(x, {})
            cD, cd = v.get(key6l, 0), v.get(key6, 0)
            ok = (D == cD and d == cd); tally[kind][0] += ok; tally[kind][1] += 1
            parts.append('x=%d %d/%d, %d/%d%s' % (x, D, cD, d, cd, '' if ok else ' DIFFERENT'))
        say('      %-10s: %s' % (kind, '; '.join(parts)))
say(); say('agreeing (system, weight) pairs, both levels at once: the class %d of %d; the whole eigenspace %d of %d' % (
    tally['class'][0], tally['class'][1], tally['eigenspace'][0], tally['eigenspace'][1]))
