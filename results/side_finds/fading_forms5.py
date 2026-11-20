# A dig, 14 Nov 2026 (README "607's quiet braids"). The same forms, weights and per-weight computation as fading_forms4.py
# (which is kept unchanged), with what fading_forms4.py lacked:
#   1. Each weight's G, the factor of T_ℓ's polynomial whose roots are root/ℓ^S for a group's roots of slope S, is printed
#      in full and written to a .jsonl file as soon as it is found. A stopped run loses nothing, and the weights can run as
#      separate processes (--from=x with X = x), joined afterwards with --fit.
#   2. The carry to weight −2 is done two ways, and both are printed:
#      plain : Lagrange of each coefficient at x = −1 (fading_forms4.py's way);
#      scaled: Lagrange of (coefficient of y^i)/(x − x0)^j, j = d − i, then times (−1 − x0)^j.
#      A quiet root is 0 at x0 = 0, the weight ℓ − 3, where the class has no quiet form (at 29 and 31 its first digit is a
#      constant times x at all 7 weights). So a coefficient that multiplies j roots grows like x^j, which a curve through a
#      few weights cannot follow for j ≥ 3. The scaled carry is the answer for slope S ≥ 1 on the small row.
#   3. Checks printed with them: the first digit of every scaled coefficient at every weight (the same at all weights if
#      the picture holds); the carries from every subset of one weight fewer, and the digits on which they agree with the
#      full carry; the factors of each carry on lane ℓ.
#   4. The answer split into its strands on the first digit, in its own field F_ℓ[t]/(f) for each factor f: t, t^ℓ,
#      t^(ℓ²), …, each checked to be a root, all distinct, the ℓ-th power stepping through them; their product and sum.
# No partition number is used. COUNTED (the forms), DERIVED (the carry to −2).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/fading_forms5.py ℓ X N [--slope=S]
#              [--from=x] [--extra=E] [--minp=P] [--disc] [--x0=x0] [--keep] [--dir=folder] [--threads=t]
#          .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/fading_forms5.py ℓ --fit=a.jsonl,b.jsonl
#              [--minp=P] [--x0=x0] [--keep] [--dir=folder]
import os, sys, time, json, itertools
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat, nmod_poly
import forms6_flint as F


def opt(name, default):
    return next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--' + name + '=')), default)


flint.ctx.threads = int(opt('threads', 8))
nums = [int(a) for a in sys.argv[1:] if a.isdigit()]; L = nums[0]
FIT = opt('fit', None); MINP = int(opt('minp', 4)); X0F = int(opt('x0', 0))
KEEP = '--keep' in sys.argv; DIR = opt('dir', HERE)
rows = []                                   # one per group per weight: what the .jsonl holds
readlog = []
if FIT:
    for p in FIT.split(','):
        p = p if os.path.isabs(p) else os.path.join(DIR, p)
        with open(p, encoding='utf-8') as f: got = [json.loads(ln) for ln in f if ln.strip()]
        rows += got; readlog.append('   read %s: %d group%s at x = %s' % (os.path.basename(p), len(got), '' if len(got) == 1 else 's', sorted({r['x'] for r in got})))
    assert rows and all(r['l'] == L for r in rows), 'the .jsonl files must hold ℓ = %d' % L
    S, WSTEP = rows[0]['slope'], rows[0]['wstep']
    assert all(r['slope'] == S and r['wstep'] == WSTEP for r in rows), 'the files mix slopes or weight steps'
    DISC = WSTEP != L - 1
else:
    X, N = nums[1], nums[2]; EXTRA = int(opt('extra', 20)); NP = N + EXTRA; XS = int(opt('from', 1))
    S = int(opt('slope', 1)); DISC = '--disc' in sys.argv; WSTEP = (L - 1) * L if DISC else L - 1
SCALED = S >= 1 and not DISC and '--plain' not in sys.argv   # --plain (14 Nov): for a group that is not 0 at x0 (7, 11)
base = os.path.join(DIR, 'fading_forms5_%d_s%d%s_%s' % (L, S, '_disc' if DISC else '', 'fit' if FIT else 'x%d-%d' % (XS, X)))
OUT, JS = base + '.txt', base + '.jsonl'
if KEEP and (os.path.exists(OUT) or (not FIT and os.path.exists(JS))):   # never rename: the join must read what it names
    sys.exit('stopped before computing anything: %s already exists, and results are never overwritten or renamed. Use another '
             '--dir, or check whether that earlier run is the one meant.' % (OUT if os.path.exists(OUT) else JS))
fh = open(OUT, 'x', encoding='utf-8') if KEEP else None
jh = open(JS, 'x', encoding='utf-8') if KEEP and not FIT else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def dig(v, n): return [(v // L ** i) % L for i in range(n)]


def agree(u, w, n):                         # how many lowest digits of u and w are the same, up to n
    k = 0
    while k < n and (u // L ** k) % L == (w // L ** k) % L: k += 1
    return k


# ---- the computation at each weight: fading_forms4.py's, unchanged, plus the printing and the .jsonl line ----
def sign(p, k): return leg(p, L) ** (1 + (k - (L - 3)) // (L - 1))      # the class's sign of U_p at weight k


def compose(Y, R, Z, Rz): return Y * Z, [R[j] for j in Rz]


def segment_factor(chi, s, Nn):
    """The Newton polygon of chi (lowest first, mod ℓ^Nn); for its segment of slope s: (G monic lowest first, whose roots
    are root/ℓ^s for the roots of valuation s, mod ℓ^prec; prec; exact?; the polygon's slopes)."""
    d = len(chi) - 1; m = L ** Nn
    def val(x):
        x %= m
        if x == 0: return Nn
        v = 0
        while x % L == 0: x //= L; v += 1
        return v
    vs = [val(c) for c in chi]; hull = [0]
    for j in range(1, d + 1):
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            if (vs[b] - vs[a]) * (j - a) >= (vs[j] - vs[a]) * (b - a): hull.pop()
            else: break
        hull.append(j)
    slopes = [(Fraction(vs[a] - vs[b], b - a), b - a) for a, b in zip(hull, hull[1:])]
    seg = [(a, b) for a, b in zip(hull, hull[1:]) if Fraction(vs[a] - vs[b], b - a) == s]
    if not seg: return None, 0, False, slopes
    j1, j2 = seg[0]; exact = vs[j1] < Nn and vs[j2] < Nn; c = vs[j1] + s * j1; prec = Nn - c
    if not exact or prec <= 0: return None, 0, exact, slopes
    mp = L ** prec
    h = [((chi[j] % m) * L ** (s * j - c)) % mp if s * j - c >= 0 else ((chi[j] % m) // L ** (c - s * j)) % mp for j in range(d + 1)]
    lc = h[j2] % L; inv = pow(lc, -1, L)
    g = [h[j] % L * inv % L for j in range(j1, j2 + 1)]
    G, _ = F.hensel_factor(h, g, [0] * j1 + [lc], L, prec)
    return [x % mp for x in G], prec, True, slopes


def compute():
    QS = [q for q in (5, 7, 11, 13) if q != L]
    for x in range(XS, X + 1):
        k = -2 + WSTEP * (x + 1); t = time.time(); m = L ** NP; n = k + 1
        H = F.hecke_matrices(k, m, [('T', q) for q in QS] + [('T', L), ('U', 2), ('U', 3)])
        c2 = sign(2, k) * 2 ** ((k - 2) // 2) % m; c3 = sign(3, k) * 3 ** ((k - 2) // 2) % m
        Y, R, r = F.kernel_part(H[('U', 2)] - F.identity(n, m) * c2, L, NP)
        if r == 0:                          # added 14 Nov for 5, 7, 11: an empty space (no effect where forms exist)
            say('   weight %d (dimension %d): no forms with the class\'s sign at 2  (%.0f s)' % (k, n, time.time() - t)); continue
        Z, Rz, r2 = F.kernel_part(F.restrict(H[('U', 3)], Y, R, r) - F.identity(r, m) * c3, L, NP)
        if r2 == 0:
            say('   weight %d (dimension %d): no forms with the class\'s signs at 2 and 3  (%.0f s)' % (k, n, time.time() - t)); continue
        Y, R = compose(Y, R, Z, Rz); r = r2
        U2W = F.restrict(H[('U', 2)], Y, R, r); U3W = F.restrict(H[('U', 3)], Y, R, r)
        Bm = [[(int(U2W[i, j]) - (c2 if i == j else 0)) % m for j in range(r)] for i in range(r)] + \
             [[(int(U3W[i, j]) - (c3 if i == j else 0)) % m for j in range(r)] for i in range(r)]
        ker, N1 = F.padic_kernel(Bm, 2 * r, r, L, NP); d = len(ker); m1 = L ** N1
        Xn = [[ker[c][i] % m1 for c in range(d)] for i in range(r)]
        Xm = F.mat_mod(r, d, [Xn[i][c] for i in range(r) for c in range(d)], m1)
        Rn = F.pivot_rows(nmod_mat(r, d, [Xn[i][c] % L for i in range(r) for c in range(d)], L), L)
        XR = F.inverse_mod(F.mat_mod(d, d, [Xn[Rn[a]][b] for a in range(d) for b in range(d)], m1), L, N1)
        def on_new(T):
            TX = F.mat_mod(r, r, [int(v) % m1 for v in T.entries()], m1) * Xm
            return XR * F.mat_mod(d, d, [int(TX[Rn[a], b]) for a in range(d) for b in range(d)], m1)
        Mt = on_new(F.restrict(H[('T', L)], Y, R, r)); Mq = {q: on_new(F.restrict(H[('T', q)], Y, R, r)) for q in QS}
        pieces = [((), F.identity(d, m1), list(range(d)), d)]
        for q in QS:
            newp = []
            for lab, Yp, Rp, rp in pieces:
                Mqr = F.restrict(Mq[q], Yp, Rp, rp)
                for f, ee in F.reduce_mat(Mqr, L).charpoly().factor()[1]:
                    co = [int(c) for c in f.coeffs()]
                    B = Mqr - F.identity(rp, m1) * ((-co[0]) % L) if len(co) == 2 else F.mat_poly(co, Mqr)
                    Z, Rz, r2 = F.kernel_part(B, L, N1)
                    if r2: Y2, R2 = compose(Yp, Rp, Z, Rz); newp.append((lab + (((-co[0]) % L) if len(co) == 2 else tuple(co),), Y2, R2, r2))
            pieces = newp
        say('   weight %d (dimension %d): %d forms new at 2 and 3 with the class\'s signs (to ℓ^%d), in %d systems  (%.0f s)' % (k, n, d, N1, len(pieces), time.time() - t))
        for lab, Yp, Rp, rp in pieces:
            Mtp = F.restrict(Mt, Yp, Rp, rp)
            chip = F.berkowitz([[int(Mtp[i, j]) for j in range(rp)] for i in range(rp)], m1)
            G, prec, exact, sl = segment_factor(chip, S, N1)
            degs = [1 if not isinstance(v, tuple) else len(v) - 1 for v in lab]
            line = '      a system of degree %s on lane ℓ: %d form%s; slopes %s' % (degs[0], rp, '' if rp == 1 else 's', ', '.join('%s×%d' % (a_, b_) for a_, b_ in sl))
            if '--labels' in sys.argv: line += '  [system %s]' % (lab,)            # added 14 Nov: name each system
            if G is not None:
                Gb = nmod_poly([c % L for c in G], L); fdeg = sorted(ff.degree() for ff, e in Gb.factor()[1] for _ in range(e))
                line += ';  its %d root%s of slope %d: G known to ℓ^%d, factors on lane ℓ of degrees %s' % (len(G) - 1, '' if len(G) == 2 else 's', S, prec, fdeg)
            elif exact is False and any(a_ == S for a_, b_ in sl): line += ';  its slope-%d segment NOT exact at this precision' % S
            say(line)
            if G is not None:
                say('         G (roots root/ℓ^%d), coefficients lowest first, each as digits lowest first: %s' % (S, [dig(c, prec) for c in G]))
                row = {'l': L, 'slope': S, 'wstep': WSTEP, 'k': k, 'x': x, 'system': repr(lab), 'deg': len(G) - 1,
                       'prec': prec, 'G': [int(c) for c in G]}
                rows.append(row)
                if jh: jh.write(json.dumps(row) + '\n'); jh.flush(); os.fsync(jh.fileno())


# ---- the carry to weight −2, the checks, and the split into strands ----
def lag(pts, vals, m):                      # Lagrange through the points (x, vals[x]), read at x = −1, modulo m
    v = 0
    for a in pts:
        nu = de = 1
        for b in pts:
            if b != a: nu *= (-1 - b); de *= (a - b)
        v = (v + vals[a] * nu * pow(de % m, -1, m)) % m
    return v


def split(co):
    Pq = lambda c: nmod_poly([int(v) % L for v in c], L)
    def powmod(a, e, f):
        r, b = Pq([1]), a % f
        while e:
            if e & 1: r = (r * b) % f
            b = (b * b) % f; e >>= 1
        return r
    def at(f, s):                           # f(s), computed in F_ℓ[t]/(f)
        acc = Pq([0])
        for g in reversed([int(v) for v in f.coeffs()]): acc = (acc * s + Pq([g])) % f
        return acc
    def vec(a, k):
        c = [int(v) for v in a.coeffs()]; return c + [0] * (k - len(c))
    say('   the answer split into its strands, on the first digit:')
    for f, e in sorted(Pq(co).factor()[1], key=lambda fe: fe[0].degree()):
        k = f.degree(); fl = [int(v) for v in f.coeffs()]; rep = '' if e == 1 else ' (repeated %d times)' % e
        if k == 1:
            say('      a lone root, an ordinary number on lane ℓ: %d%s' % ((-fl[0]) % L, rep)); continue
        s = [Pq([0, 1]) % f]
        for i in range(k): s.append(powmod(s[-1], L, f))
        prod, sm = Pq([1]), Pq([0])
        for x in s[:k]: prod = (prod * x) % f; sm = (sm + x) % f
        say('      a braid of %d strands%s, in the field with ℓ^%d elements; t is a root of %s (lowest first)' % (k, rep, k, fl))
        say('         every strand a root: %s; all %d distinct: %s; the ℓ-th power steps through them and returns: %s' % (
            all(at(f, x).degree() < 0 for x in s[:k]), k, len({tuple(vec(x, k)) for x in s[:k]}) == k, s[k] == s[0]))
        say('         the strands multiply to %s and add up to %s' % (vec(prod, k)[0] if prod.degree() <= 0 else vec(prod, k),
                                                                     vec(sm, k)[0] if sm.degree() <= 0 else vec(sm, k)))
        for i, x in enumerate(s[:k]): say('         strand %d = t^(ℓ^%d) = %s' % (i + 1, i, vec(x, k)))


def carry_all():
    s3 = leg(3, L); groups = {}
    for r in rows: groups.setdefault(r['system'], []).append(r)
    for sysl, allr in groups.items():
        rs = sorted([r for r in allr if r['prec'] >= MINP], key=lambda r: r['x']); xs = [r['x'] for r in rs]
        if len(set(xs)) != len(xs):
            say(); say('a group in the system %s: the same weight appears twice, x = %s (two files of one run?): NOT carried' % (sysl, xs)); continue
        if len(rs) < 2 or len({r['deg'] for r in rs}) != 1:
            say(); say('a group in the system %s, at x = %s: not carried (it needs two or more weights of one degree, each known to ℓ^%d)' % (
                sysl, sorted(r['x'] for r in allr), MINP)); continue
        d = rs[0]['deg']; P = min(r['prec'] for r in rs); m = L ** P
        G = {r['x']: [c % m for c in r['G']] for r in rs}
        say(); say('a group of %d root%s of slope %d, weights x = %s (k = %s), carried to weight −2 modulo ℓ^%d' % (
            d, '' if d == 1 else 's', S, xs, [r['k'] for r in rs], P))
        say('   its system on lane ℓ: %s' % sysl)
        bad = []
        if SCALED:
            for i in range(d):
                j = d - i; f = [G[x][i] * pow((x - X0F) % L, -j, L) % L for x in xs]
                say('   coefficient of y^%d divided by (x − %d)^%d, first digit at x = %s: %s' % (i, X0F, j, xs, f))
                if len(set(f)) > 1: bad.append(i)
            say('   so the same at every weight for %d of %d coefficients%s' % (d - len(bad), d, '' if not bad else '; NOT for i = %s' % bad))
        def carry(pts, scaled):
            out = []
            for i in range(d + 1):
                j = d - i
                if scaled:
                    v = lag(pts, {x: G[x][i] * pow((x - X0F) % m, -j, m) % m for x in pts}, m) * pow((-1 - X0F) % m, j, m) % m
                else:
                    v = lag(pts, {x: G[x][i] for x in pts}, m)
                out.append(v * pow(s3 % m, j, m) % m)          # the polynomial of (3|ℓ)·y, y = root/ℓ^S
            return out
        answer = None
        for name, sc in (('plain (fading_forms4.py\'s way)', False), ('scaled', True)):
            if sc and not SCALED: continue
            co = carry(xs, sc); answer = (co, name)
            # how many digits the weights support: those on which every carry from one weight fewer agrees with this one
            # (confirmed), plus one (likely); nothing deeper is printed, since nothing deeper is supported
            if len(xs) >= 3:
                ag = min(agree(a, b, P) for s in itertools.combinations(xs, len(xs) - 1) for a, b in zip(carry(list(s), sc), co))
                shown = min(P, ag + 1)
                note = 'the first %d digit%s confirmed (the carries from every %d of the %d weights agree with it), digit %d likely' % (
                    ag, '' if ag == 1 else 's', len(xs) - 1, len(xs), ag + 1)
            else:
                shown = min(P, len(xs)); note = 'only %d weights, so no cross-check: at most %d digits' % (len(xs), len(xs))
            say('   %s carry: the polynomial of the constants (3|ℓ)·u/ℓ^%d at weight −2, coefficients lowest first, each as its first %d digit%s (lowest first):' % (
                name, S, shown, '' if shown == 1 else 's'))
            say('      %s' % [dig(c, shown) for c in co])
            say('      %s' % note)
            fb = nmod_poly([c % L for c in co], L)
            say('      its factors on lane ℓ: degrees %s%s' % (sorted(f.degree() for f, e in fb.factor()[1] for _ in range(e)),
                                                         ';  the root: %s' % dig((-co[0]) % m, shown) if d == 1 else ''))
        if bad:
            say('   NO ANSWER for this group: the first-digit check failed at i = %s, so a coefficient does not grow like x^j and' % bad)
            say('   neither carry above is to be used. Not split.'); continue
        say('   the answer is the %s carry' % answer[1].split(' ')[0])
        split(answer[0])


t0 = time.time()
if FIT:
    say('ℓ = %d, slope %d: the groups kept in these files, carried to weight −2' % (L, S))
    for ln in readlog: say(ln)
else:
    say('ℓ = %d, slope %d, weights −2 + %d·(x + 1) for x = %d..%d, modulo ℓ^%d; FLINT %s; each weight\'s G also goes to %s' % (
        L, S, WSTEP, XS, X, NP, flint.__version__, os.path.basename(JS) if jh else '(nothing: run without --keep)'))
    compute()
carry_all()
say(); say('seconds: %.0f' % (time.time() - t0))
