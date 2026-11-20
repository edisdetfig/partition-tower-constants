# A dig, 13 Nov 2026: 607'S PAIR ONE DIGIT DEEPER IN THE TOWER ITSELF (README "607's pair, one digit deeper in the tower
# itself"; predicted there first; prompted by a review of the page "The Exception 607"). At the weights k = 604 + 606·j
# (j = 0, 1, 2), modulo 607³: the class (forms6_flint.class_padic, its signs carried along the family), T₆₀₇ on it, the
# piece where T₆₀₇ ≡ 367 (Fitting), C = (T₆₀₇ − 367)/607 there and T₅ there. Form A is followed by its T₅ number on lane
# 607, the root t of x² + 561x + 66 in F_607[t]/(t² + 561t + 66) (as in chords_page_data607.py), and its δ is C's number on
# T₅'s eigenvector. Then δ(−2) = 2·δ(604) − δ(1210) (Hida: along j the n-th differences are divisible by 607ⁿ), checked by
# the second difference through 1816. The tower's constants are (3|607)·u(−2) = −(367 + 607·δ(−2)) = 240 + 607·(−1 − δ(−2))
# modulo 607². Their second digits are also given in size and direction in the fixed model of chords_page_data607b.py
# (F_607(√3), t ↦ τ = 23 + 137·√3, directions from G = 1 + 2·√3). Writes partitions_flint_607_deeper.json.
# COUNTED (the forms' numbers), DERIVED (the value at weight −2).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_607_deeper.py [--keep]
import os, sys, time, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
import forms6_flint as F
flint.ctx.threads = 8
KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_flint_607_deeper.txt'), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


L = 607; N = 3; m = L ** N; A0 = 367; C0, C1 = 66, 561                  # t² + 561t + 66: T₅ on the pair, on lane 607


def mul(x, y):                                                            # F_607[t]/(t² + C1·t + C0), x = (a, b) = a + b·t
    a, b = x; c, d = y; bd = b * d
    return ((a * c - bd * C0) % L, (a * d + b * c - bd * C1) % L)


def pw(x, n):
    r = (1, 0)
    while n:
        if n & 1: r = mul(r, x)
        x = mul(x, x); n >>= 1
    return r


inv = lambda x: pw(x, L * L - 2)
sub = lambda x, y: ((x[0] - y[0]) % L, (x[1] - y[1]) % L)
scal = lambda c, x: (c * x[0] % L, c * x[1] % L)
t0 = time.time(); delta = {}
for j in range(3):
    k = 604 + 606 * j; t = time.time(); e = 1 + j
    s2, s3 = leg(2, L) ** e, leg(3, L) ** e
    H = F.hecke_matrices(k, m, [('U', 2), ('U', 3), ('T', L), ('T', 5)])
    W, RW, r, Wn, Rn, d, Np = F.class_padic(H, k, L, N, s2, s3)
    mp = L ** Np
    TL = F.restrict(F.mat_mod(k + 1, k + 1, [int(v) % mp for v in H[('T', L)].entries()], mp), Wn, Rn, d)
    T5 = F.restrict(F.mat_mod(k + 1, k + 1, [int(v) % mp for v in H[('T', 5)].entries()], mp), Wn, Rn, d)
    co = [int(c) for c in F.reduce_mat(TL, L).charpoly().coeffs()]; zero = next(i for i, c in enumerate(co) if c)
    Y, R, s = F.kernel_part(TL - F.identity(d, mp) * A0, L, Np)
    B = F.restrict(TL, Y, R, s); B5 = F.restrict(T5, Y, R, s)
    Bl = [[int(B[i, jj]) for jj in range(s)] for i in range(s)]
    semis = all((Bl[i][jj] - (A0 if i == jj else 0)) % L == 0 for i in range(s) for jj in range(s))
    say('weight %d (signs U₂ %+d, U₃ %+d): the class has %d forms (to 607^%d), %d of them ordinary; T₆₀₇ ≡ 367 on a piece of '
        'dimension %d; T₆₀₇ = 367·I there on lane 607: %s  (%.0f s)' % (k, s2, s3, d, Np, d - zero, s, semis, time.time() - t))
    if s != 2 or not semis or Np < 2: say('   stopped: not the expected pair'); break
    C = [[((Bl[i][jj] - (A0 if i == jj else 0)) // L) % L for jj in range(2)] for i in range(2)]
    b = [[int(B5[i, jj]) % L for jj in range(2)] for i in range(2)]
    chi5 = [(b[0][0] * b[1][1] - b[0][1] * b[1][0]) % L, (-(b[0][0] + b[1][1])) % L, 1]
    chiC = [(C[0][0] * C[1][1] - C[0][1] * C[1][0]) % L, (-(C[0][0] + C[1][1])) % L, 1]
    say('   T₅ on the pair, on lane 607: %s (predicted [66, 561, 1]); C = (T₆₀₇ − 367)/607: %s' % (chi5, chiC))
    if chi5 != [C0, C1, 1]: say('   stopped: T₅\'s polynomial changed'); break
    tt = (0, 1)
    v = ((b[0][1], 0), sub(tt, (b[0][0], 0))) if b[0][1] else (sub(tt, (b[1][1], 0)), (b[1][0], 0))
    Cv = [((C[i][0] * v[0][0] + C[i][1] * v[1][0]) % L, (C[i][0] * v[0][1] + C[i][1] * v[1][1]) % L) for i in range(2)]
    i0 = 0 if v[0] != (0, 0) else 1
    dA = mul(Cv[i0], inv(v[i0]))
    assert mul(dA, v[1 - i0]) == Cv[1 - i0]                                                     # v is C's eigenvector too
    assert (lambda z: ((mul(z, z)[0] + chiC[1] * z[0] + chiC[0]) % L, (mul(z, z)[1] + chiC[1] * z[1]) % L))(dA) == (0, 0)
    delta[k] = dA
    say('   form A (T₅ ≡ t): δ = %d + %d·t; form B: δ = %d + %d·t' % (dA[0], dA[1], *pw(dA, L)))
if len(delta) == 3:
    d0, d1, d2 = delta[604], delta[1210], delta[1816]
    sec = sub(sub(d2, scal(2, d1)), scal(-1, d0))
    say(); say('second difference δ(1816) − 2·δ(1210) + δ(604) = %d + %d·t (predicted 0)' % sec)
    dm2 = sub(scal(2, d0), d1)
    say('δ(−2) = 2·δ(604) − δ(1210) = %d + %d·t for form A; %s' % (dm2[0], dm2[1], 'in F₆₀₇: the two forms AGREE on the second digit too'
                                                                     if dm2[1] == 0 else 'not in F₆₀₇: the two forms PART on the second digit'))
    secA = sub((L - 1, 0), dm2); secB = pw(secA, L)
    say('the tower\'s constants modulo 607²: 240 + 607·(%d + %d·t) and 240 + 607·(%d + %d·t)' % (*secA, *secB))
    # size and direction in the fixed model of chords_page_data607b.py
    tau = (23, 137); m3 = lambda x, y: ((x[0] * y[0] + 3 * x[1] * y[1]) % L, (x[0] * y[1] + x[1] * y[0]) % L)
    def p3(x, n):
        rr = (1, 0)
        while n:
            if n & 1: rr = m3(rr, x)
            x = m3(x, x); n >>= 1
        return rr
    assert m3(tau, tau) == ((-C1 * tau[0] - C0) % L, (-C1 * tau[1]) % L)
    emb = lambda ab: ((ab[0] + ab[1] * tau[0]) % L, (ab[1] * tau[1]) % L)
    E = p3((1, 2), L - 1); dtab = {}; cur = (1, 0)
    for i in range(L + 1): dtab[cur] = i; cur = m3(cur, E)
    sd = lambda u: [p3(u, L + 1)[0], dtab[p3(u, L - 1)]]
    pairT = {'A': sd(emb(secA)), 'B': sd(emb(secB))}
    pairF = {'A': sd(emb(d0)), 'B': sd(emb(pw(d0, L)))}
    say('in size and direction (F₆₀₇(√3), from G = 1 + 2·√3): the tower\'s second digits: form A %s, form B %s; '
        'the forms\' C-numbers at weight 604: form A %s, form B %s' % (pairT['A'], pairT['B'], pairF['A'], pairF['B']))
    if KEEP: json.dump({'tower_second_digit': pairT, 'forms604': pairF, 'delta': {str(kk): list(vv) for kk, vv in delta.items()},
                        'delta_m2': list(dm2), 'parts': dm2[1] != 0, 'second_difference': list(sec)},
                       open(os.path.join(HERE, 'partitions_flint_607_deeper.json'), 'x'), separators=(',', ':'))
say(); say('seconds: %.0f' % (time.time() - t0))
