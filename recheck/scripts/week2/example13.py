# Worked example at 13: one level of the tower against the level two above it.
# Needs python-flint (pip install python-flint). Runs in about a second.
from flint import fmpz

l = 13
def delta(b):                       # 24*delta = 1 (mod l^b), 0 <= delta < l^b
    return pow(24, -1, l**b)
def H(b, n):                        # exact partition number p(l^b n + delta_b)
    return int(fmpz(l**b * n + delta(b)).partitions_p())

print("delta_1, delta_3, delta_5 =", delta(1), delta(3), delta(5))
h1 = [H(1, n) for n in range(8)]; h3 = [H(3, n) for n in range(8)]; h5 = [H(5, n) for n in range(8)]
print("level 1 mod 13:", [x % 13 for x in h1])
print("level 3 mod 13:", [x % 13 for x in h3])
print("6 x level 1   :", [6 * x % 13 for x in h1])
# second digit: the constant is 6 + 3*13 = 45 (mod 169); the law mod 169 starts at level 3
print("level 5 mod 169 :", [x % 169 for x in h5])
print("45 x level 3    :", [45 * x % 169 for x in h3])
print("first digit holds :", all((a - 6 * b) % 13 == 0 for a, b in zip(h3, h1)))
print("two digits hold   :", all((a - 45 * b) % 169 == 0 for a, b in zip(h5, h3)))
print("45 x level 1 = level 3 mod 169?", all((a - 45 * b) % 169 == 0 for a, b in zip(h3, h1)))
