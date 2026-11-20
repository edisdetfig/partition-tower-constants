# A helper, 12 Nov 2026: the Appell–Lerch numerator of Ramanujan's f(q), exactly, for the WSL towers (WSL has FLINT for
# n < 10⁸, but no numpy). A(q) = 1 + 4·Σ_(n≥1) (−1)ⁿ q^(n(3n+1)/2)/(1 + qⁿ), so f(q) = A(q)/(q;q)_∞
# (mock_theta_tower.py checks the identity). The coefficients are small integers (|A_m| ≤ 4·#{n : n(3n+1)/2 ≤ m}), written
# as raw little-endian int32 to the given file. COUNTED (exact).
#   usage: .venv/Scripts/python.exe -I -B -u -X utf8 -W ignore results/side_finds/mock_theta_appell.py N <out file>
import sys, time
import numpy as np
N = int(sys.argv[1]); OUT = sys.argv[2]
t = time.time()
a = np.zeros(N, dtype=np.int64); a[0] = 1; n = 1
while n * (3 * n + 1) // 2 < N:
    s = n * (3 * n + 1) // 2; cnt = len(range(s, N, n))
    pat = np.empty(cnt, dtype=np.int64); pat[0::2] = 1; pat[1::2] = -1
    a[s::n] += 4 * (-1) ** n * pat; n += 1
assert np.abs(a).max() < 2 ** 31
a.astype('<i4').tofile(OUT)
print('A(q) to %d terms, %d values of n, largest |coefficient| %d; first values %s; written to %s (%.0f s)' % (
    N, n - 1, int(np.abs(a).max()), a[:12].tolist(), OUT, time.time() - t))
