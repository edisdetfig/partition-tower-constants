# Minimum ell-adic valuation of p(ell^b n + delta_b), 24 delta_b = 1 mod ell^b, for ell = 5, 7, 11.
# Checks whether the classical exponents (5^b, 7^(floor(b/2)+1), 11^b) are attained (sharp).
import time
from sympy import npartitions

def v(x, l):
    if x == 0:
        return 10**9
    c = 0
    while x % l == 0:
        x //= l; c += 1
    return c

for l, bmax in ((5, 6), (7, 8), (11, 4)):
    for b in range(1, bmax + 1):
        M = l**b
        d = pow(24, -1, M)          # delta_b
        ns = range(0, 6) if M * 6 < 3_000_000 else range(0, 2)
        t = time.time()
        vals = [v(npartitions(M * n + d), l) for n in ns]
        print(f"ell={l} b={b} delta={d} n={list(ns)} v={vals} min={min(vals)}  ({time.time()-t:.1f}s)", flush=True)
