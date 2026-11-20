# A dig, 12 Nov 2026: THE EXCEPTION 607. Boylan–Webb (Trans. AMS 2013; reports/Partition tower literature check.md)
# computed the rank of the partition tower's stable space to be ⌊(ℓ + 12)/24⌋ for every prime up to 1297 except 607.
# PREDICTED before the run (my reading): at weight ℓ − 3 = 604 the tower's class still has all ⌊(607 + 12)/24⌋ = 25 forms,
# and exactly one of them is non-ordinary at 607 (T_607 ≡ 0 on it, on lane 607), so one form fewer survives the levels.
# Controls, predicted all ordinary: the neighbouring primes 599, 601, 613.
# Method (forms6_flint.py, FLINT only, everything modulo ℓ): M_{ℓ−3}(Γ0(6)) over the field with ℓ elements; the class as
# the joint generalized eigenspace of U₂ and U₃ at (2|ℓ)·2^((k−2)/2) and (3|ℓ)·3^((k−2)/2); T_ℓ on the class; its
# characteristic polynomial, factored: the multiplicity of the root 0 counts the non-ordinary forms. COUNTED.
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/partitions_flint_607.py [ℓ …] [--keep]
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
try:
    import psutil; psutil.Process().nice(psutil.IDLE_PRIORITY_CLASS if os.name == 'nt' else 19)
except Exception: pass
import flint
from flint import nmod_mat
import forms6_flint as F
flint.ctx.threads = 4
NAMED = [int(a) for a in sys.argv[1:] if a.isdigit()]; KEEP = '--keep' in sys.argv
fh = open(os.path.join(HERE, 'partitions_flint_607%s.txt' % ''.join('_%d' % x for x in NAMED)), 'x', encoding='utf-8') if KEEP else None


def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True)
    if fh: fh.write(line + '\n'); fh.flush(); os.fsync(fh.fileno())


def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


t0 = time.time()
for L in NAMED or [599, 601, 607, 613]:
    t = time.time(); k = L - 3; dexp = (L + 12) // 24
    H = F.hecke_matrices(k, L, [('U', 2), ('U', 3), ('T', L)]); n = k + 1
    c2 = leg(2, L) * pow(2, (k - 2) // 2, L) % L; c3 = leg(3, L) * pow(3, (k - 2) // 2, L) % L
    U2, U3, TL = H[('U', 2)], H[('U', 3)], H[('T', L)]
    K1, r1 = F.generalized_kernel(U2 - F.identity(n, L) * c2); X1, R1 = F.normalized_basis(K1, r1, L)
    M3 = F.restrict(U3, X1, R1, r1); K2, r0 = F.generalized_kernel(M3 - F.identity(r1, L) * c3)
    Xc0 = X1 * nmod_mat(r1, r0, [int(K2[i, j]) for i in range(r1) for j in range(r0)], L)
    Xc, Rc = F.normalized_basis(Xc0, r0, L)
    chi = F.restrict(TL, Xc, Rc, r0).charpoly(); fac = chi.factor()[1]
    zero = sum(e for f, e in fac if f.degree() == 1 and int(f.coeffs()[0]) == 0)
    degs = sorted(f.degree() for f, e in fac for _ in range(e))
    say('ℓ = %d (weight %d, dimension %d): the class has %d forms (⌊(ℓ + 12)/24⌋ = %d); T_ℓ on it: %d with the root 0 on lane ℓ '
        '(non-ordinary), the rest in irreducible pieces of degrees %s  (%.0f s)' % (
            L, k, n, r0, dexp, zero, [g for g in degs if not (g == 1 and zero and False)], time.time() - t))
say(); say('seconds: %.0f' % (time.time() - t0))
