"""Certify that 701 cannot be the least split prime in the exponent-seven gap."""

from math import isqrt
import numpy as np

p = 701
lower = 31 * 10**19
constant = 4 * p**7
assert constant < 4 * lower  # q>lower implies y=1
limit = isqrt(constant - lower - 1)
primes = [n for n in range(3, p, 2) if all(n % d for d in range(3, isqrt(n)+1, 2))]


def good_table(ell):
    return np.array([pow((a*a-constant) % ell, (ell-1)//2, ell) == ell-1
                     for a in range(ell)], dtype=bool)


modulus = 2
residues = np.array([1], dtype=np.int64)
for ell in primes:
    if ell > 23:
        break
    candidates = np.concatenate([residues+k*modulus for k in range(ell)])
    residues = candidates[good_table(ell)[candidates % ell]]
    modulus *= ell
print('trace bound, wheel modulus, residue count', limit, modulus, len(residues))

tables = [(ell, good_table(ell)) for ell in primes if ell > 23]
counts = np.zeros(1+len(tables), dtype=np.int64)
for offset in range(0, limit+1, modulus):
    values = residues + offset
    values = values[values <= limit]
    counts[0] += values.size
    for i, (ell, table) in enumerate(tables, 1):
        values = values[table[values % ell]]
        counts[i] += values.size
        if not values.size:
            break
for ell, count in zip([23]+[ell for ell, table in tables], counts):
    print(ell, count)
    if count == 0:
        break
assert counts[-1] == 0
